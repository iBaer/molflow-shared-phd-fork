//
// Created by pascal on 5/17/21.
//

#include "SettingsIO.h"
#include <Helper/ConsoleLogger.h>
#include <Helper/StringHelper.h>
#include <filesystem>

// zip
#include <File.h>
#include <ziplib/ZipArchive.h>
#include <ziplib/ZipFile.h>

namespace SettingsIO {
    bool overwrite = false;
    bool isArchive = false;
    std::string workFile;
    std::string workPath;
    std::string inputFile;
    std::string inputPath;
    std::string outputFile;
    std::string outputPath;

    const std::string supportedFileFormats[]{".xml", ".zip"};

    int prepareIO() {
        if (initDirectories()) {
            Log::console_error("Error preparing folders\n");
            return 1;
        }

        if (initFromZip()) {
            Log::console_error("Error handling input from zip file\n");
            return 2;
        }
        return 0;
    }

//! Set all directory / output related variables and initialises directories depending on set parameters
    int initDirectories() {

        int err = 0;

        if (std::filesystem::path(SettingsIO::inputFile).has_parent_path()) {
            SettingsIO::inputPath =
                    std::filesystem::path(SettingsIO::inputFile).parent_path().string();
        }

        // Overwrite excludes outputpath/filename
        if (SettingsIO::overwrite) {
            SettingsIO::outputFile = SettingsIO::inputFile;
            std::string tmpFolder = MFMPI::world_size > 1 ? std::string("tmp{}/",MFMPI::world_rank) : "tmp/";
            SettingsIO::workPath = tmpFolder;
        } else if(!SettingsIO::outputPath.empty()){
            SettingsIO::workPath = SettingsIO::outputPath;
        }
        else if (SettingsIO::outputPath.empty() && std::filesystem::path(SettingsIO::outputFile).has_parent_path()) { // Use a default outputpath if unset
            SettingsIO::workPath =
                    std::filesystem::path(SettingsIO::outputFile).parent_path().string();
        } else if (SettingsIO::outputPath
                .empty()) { // Use a default outputpath if unset
            SettingsIO::outputPath = "Results_" + Util::getTimepointString();
            SettingsIO::workPath = SettingsIO::outputPath;
        }
        else if (std::filesystem::path(SettingsIO::outputFile).has_parent_path()) {
            Log::console_error(
                    "Output path was set to %s, but Output file also contains a parent "
                    "path %s\n"
                    "Output path will be appended!\n",
                    SettingsIO::outputPath.c_str(),
                    std::filesystem::path(SettingsIO::outputFile).parent_path().string().c_str());
        }

        // Use a default outputfile name if unset
        if (SettingsIO::outputFile.empty())
            SettingsIO::outputFile =
                    "out_" +
                    std::filesystem::path(SettingsIO::inputFile).filename().string();

        bool formatIsSupported = false;
        for (const auto &format : supportedFileFormats) {
            if (std::filesystem::path(SettingsIO::outputFile).extension().string() ==
                format) {
                formatIsSupported = true;
            }
        }
        if (!formatIsSupported) {
            Log::console_error("File format is not supported: %s\n",
                               std::filesystem::path(SettingsIO::outputFile).extension().string().c_str());
            return 1;
        }

        // Try to create directories
        // First for outputpath, with tmp/ and lastly ./ as fallback plans
        try {
            if (!std::filesystem::exists(SettingsIO::workPath))
                std::filesystem::create_directory(SettingsIO::workPath);
        } catch (const std::exception &e) {
            Log::console_error("Couldn't create directory [ %s ], falling back to "
                               "tmp folder for output files\n",
                               SettingsIO::workPath.c_str());
            ++err;

            // use fallback dir
            SettingsIO::workPath = MFMPI::world_size > 1 ? std::string("tmp{}/",MFMPI::world_rank) : "tmp/";
            SettingsIO::outputPath = SettingsIO::workPath;

            try {
                if (!std::filesystem::exists(SettingsIO::workPath))
                    std::filesystem::create_directory(SettingsIO::workPath);
            } catch (const std::exception &e) {
                SettingsIO::workPath = "./";
                Log::console_error("Couldn't create fallback directory [ tmp/ ], falling "
                                   "back to binary folder instead for output files\n");
                ++err;
            }
        }

        // Next check if outputfile name has parent path as name
        // Additional directory in outputpath
        if (std::filesystem::path(SettingsIO::outputFile).has_parent_path()) {
            std::string outputFilePath;
            if(SettingsIO::outputPath.empty())
                outputFilePath =
                        std::filesystem::path(SettingsIO::outputFile).parent_path().string();
            else
                outputFilePath =
                        std::filesystem::path(SettingsIO::outputPath)
                                .append(std::filesystem::path(SettingsIO::outputFile)
                                                .parent_path()
                                                .string())
                                .string();
            try {
                if (!std::filesystem::exists(outputFilePath))
                    std::filesystem::create_directory(outputFilePath);
            } catch (const std::exception &e) {
                Log::console_error(
                        "Couldn't create parent directory set by output filename [ %s ], "
                        "will only use default output path instead\n",
                        outputFilePath.c_str());
                ++err;
            }
        }

        return 0;
    }

    int initFromZip() {
        if (std::filesystem::path(SettingsIO::inputFile).extension() == ".zip") {
            SettingsIO::isArchive = true;

            // decompress file
            std::string parseFileName;
            Log::console_msg_master(2, "Decompressing zip file...\n");

            ZipArchive::Ptr zip = ZipFile::Open(SettingsIO::inputFile);
            if (zip == nullptr) {
                Log::console_error("Can't open ZIP file\n");
                return 1;
            }
            size_t numItems = zip->GetEntriesCount();
            bool notFoundYet = true;
            for (int i = 0; i < numItems && notFoundYet;
                 i++) { // extract first xml file found in ZIP archive
                auto zipItem = zip->GetEntry(i);
                std::string zipFileName = zipItem->GetName();

                if (std::filesystem::path(zipFileName).extension() ==
                    ".xml") { // if it's an .xml file
                    notFoundYet = false;

                    if (SettingsIO::outputPath != "tmp/")
                        FileUtils::CreateDir("tmp"); // If doesn't exist yet

                    parseFileName = "tmp/" + zipFileName;
                    ZipFile::ExtractFile(SettingsIO::inputFile, zipFileName, parseFileName);
                }
            }
            if (parseFileName.empty()) {
                Log::console_error("Zip file does not contain a valid geometry file!\n");
                return 1;
            }
            SettingsIO::workFile = parseFileName;
            Log::console_msg_master(2, "New input file: %s\n",
                                    SettingsIO::workFile.c_str());
        } else {
            SettingsIO::workFile = SettingsIO::inputFile;
        }
        return 0;
    }

// Out
    void cleanup_files() {
        // a) tmp folder if it is not our output folder
        if (std::filesystem::path(SettingsIO::outputPath)
                    .relative_path()
                    .compare(std::filesystem::path("tmp")) &&
            std::filesystem::path(SettingsIO::outputPath)
                    .parent_path()
                    .compare(std::filesystem::path("tmp"))) {
            // Settings::tmpfile_dir
            std::filesystem::remove_all("tmp");
        }
    }
} // namespace SettingsIO