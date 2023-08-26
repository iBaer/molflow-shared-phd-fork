/*
Program:     MolFlow+ / Synrad+
Description: Monte Carlo simulator for ultra-high vacuum and synchrotron radiation
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY / Pascal BAEHR
Copyright:   E.S.R.F / CERN
Website:     https://cern.ch/molflow

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Full license text: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
*/

#include "SettingsIO.h"
#include <Helper/ConsoleLogger.h>
#include <Helper/StringHelper.h>
#include <filesystem>

// zip
#include <File.h>
#include <ZipLib/ZipArchive.h>
#include <ZipLib/ZipFile.h>

// Input Output related settings and handy functions for the CLI application
namespace SettingsIO {
    // explanation in header
    bool overwrite = false;
    bool isArchive = false;
    bool outputFacetDetails = false;
    bool outputFacetQuantities = false;
    double autogenerateTest = 0.0;

    std::string workFile;
    std::string workPath;
    std::string inputFile;
    std::string inputPath;
    std::string outputFile;
    std::string outputPath;
    std::vector<std::string> extraFiles;
    std::map<std::string, std::vector<std::string>> cachedLines;
    std::vector<std::vector<std::string>> formulas;
    std::vector<SelectionGroup> selections;
    const std::string supportedFileFormats[]{".xml", ".zip", ".syn", ".syn7z"};

    //! prepares work/output directories and unpacks from archive
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
            std::string tmpFolder = MFMPI::world_size > 1 ? fmt::format("tmp{}/",MFMPI::world_rank) : "tmp/";
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
                    "Output path was set to {}, but Output file also contains a parent "
                    "path {}\n"
                    "Output path will be appended!\n",
                    SettingsIO::outputPath,
                    std::filesystem::path(SettingsIO::outputFile).parent_path().string());
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
            Log::console_error("File format is not supported: {}\n",
                               std::filesystem::path(SettingsIO::outputFile).extension().string());
            return 1;
        }

        // Try to create directories
        // First for outputpath, with tmp/ and lastly ./ as fallback plans
        try {
            if (!std::filesystem::exists(SettingsIO::workPath))
                std::filesystem::create_directory(SettingsIO::workPath);
        } catch (const std::exception &e) {
            Log::console_error("Couldn't create directory [ {} ], falling back to "
                               "tmp folder for output files\n",
                               SettingsIO::workPath);
            ++err;

            // use fallback dir
            SettingsIO::workPath = MFMPI::world_size > 1 ? fmt::format("tmp{}/",MFMPI::world_rank) : "tmp/";
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
                        "Couldn't create parent directory set by output filename [ {} ], "
                        "will only use default output path instead\n",
                        outputFilePath);
                ++err;
            }
        }

        return 0;
    }

    //! Unzip file and set correct variables (e.g. uncompressed work file)
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

                    std::string tmpFolder = MFMPI::world_size > 1 ? fmt::format("tmp{}/",MFMPI::world_rank) : "tmp/";
                    if (SettingsIO::outputPath != tmpFolder)
                        FileUtils::CreateDir(tmpFolder); // If doesn't exist yet

                    parseFileName = tmpFolder + zipFileName;
                    ZipFile::ExtractFile(SettingsIO::inputFile, zipFileName, parseFileName);
                }
            }
            if (parseFileName.empty()) {
                Log::console_error("Zip file does not contain a valid geometry file!\n");
                return 1;
            }
            SettingsIO::workFile = parseFileName;
            Log::console_msg_master(2, "New input file: {}\n",
                                    SettingsIO::workFile);
        } else {
            SettingsIO::workFile = SettingsIO::inputFile;
        }
        return 0;
    }

//! Cleanup tmp files
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

        if(!SettingsIO::workPath.empty() && is_empty(std::filesystem::path(SettingsIO::workPath))){
            std::filesystem::remove_all(SettingsIO::workPath);
        }
    }
} // namespace SettingsIO