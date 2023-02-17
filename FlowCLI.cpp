//
// Created by Pascal Baehr on 28.04.20.
//

#include "SimulationManager.h"
#include "SMP.h"
#include "Buffer_shared.h"
#include <fstream>
#include <filesystem>
#include <Parameter.h>
#include <IO/LoaderXML.h>
#include <IO/WriterXML.h>
#include "GeometrySimu.h"
#include "Initializer.h"
#include "Helper/MathTools.h"
#include <sstream>
#include <Helper/Chronometer.h>
#include <Helper/StringHelper.h>
#include <Helper/ConsoleLogger.h>
#include <ZipLib/ZipFile.h>

#include "FlowMPI.h"
#include "SettingsIO.h"
#include "File.h"

static constexpr const char* molflowCliLogo = R"(
  __  __     _  __ _
 |  \/  |___| |/ _| |_____ __ __
 | |\/| / _ \ |  _| / _ \ V  V /
 |_|  |_\___/_|_| |_\___/\_/\_/
    )";

/*static constexpr const char* molflowCliLogo = R"(
 __
 \ \      __  __     _  __ _
  \ \    |  \/  |___| |/ _| |_____ __ __
  / /    | |\/| / _ \ |  _| / _ \ V  V /
 /_/_____|_|  |_\___/_|_| |_\___/\_/\_/
   |_____|
    )";*/

int main(int argc, char** argv) {

#if defined(WIN32) || defined(__APPLE__)
    setlocale(LC_ALL, "C");
#else
    std::setlocale(LC_ALL, "C");
#endif

#if defined(USE_MPI)
    MFMPI::mpi_initialize();
#endif

    Log::console_msg_master(0, "%s\n", molflowCliLogo);

    SimulationManager simManager{};
    simManager.interactiveMode = true;
    std::shared_ptr<SimulationModel> model = std::make_shared<SimulationModel>();
    GlobalSimuState globState{};


    if(Initializer::initFromArgv(argc, argv, &simManager, model)){
#if defined(USE_MPI)
        MPI_Finalize();
#endif
        return 41;
    }

#if defined(USE_MPI)
    MFMPI::mpi_transfer_simu();
#endif

    if(Initializer::initFromFile(&simManager, model, &globState)){
#if defined(USE_MPI)
        MPI_Finalize();
#endif
        return 42;
    }
    size_t oldHitsNb = globState.globalHits.globalHits.nbMCHit;
    size_t oldDesNb = globState.globalHits.globalHits.nbDesorbed;

    // Get autosave file name
    std::string autoSave = Initializer::getAutosaveFile();


    //simManager.ReloadHitBuffer();
    //simManager.IncreasePriority();
    Log::console_msg_master(1,"[%s] Commencing simulation for %lu seconds from %lu desorptions.\n", Util::getTimepointString().c_str(), Settings::simDuration, globState.globalHits.globalHits.nbDesorbed);

    try {
        simManager.StartSimulation();
    }
    catch (std::runtime_error& e) {
        Log::console_error("Starting simulation: %s\n",e.what());
#if defined(USE_MPI)
        MPI_Finalize();
#endif
        return 43;
    }

    Chronometer simTimer;
    simTimer.Start();
    double elapsedTime = 0.0;

    bool endCondition = false;
    do {
        ProcessSleep(1000);

        elapsedTime = simTimer.Elapsed();
        if(model->otfParams.desorptionLimit != 0)
            endCondition = globState.globalHits.globalHits.nbDesorbed/* - oldDesNb*/ >= model->otfParams.desorptionLimit;

        if(endCondition){
            /*std::stringstream outFile;
            outFile << SettingsIO::outputPath << "/" <<"desorped_" << model->otfParams.desorptionLimit << "_" <<
                 std::filesystem::path(SettingsIO::outputFile).filename().string();*/
            std::string outFile = std::filesystem::path(SettingsIO::outputPath)
                    .append("desorped_")
                    .concat(std::to_string(model->otfParams.desorptionLimit))
                    .concat("_")
                    .concat(std::filesystem::path(SettingsIO::outputFile).filename().string()).string();
            try {
                std::filesystem::copy_file(SettingsIO::workFile, outFile, std::filesystem::copy_options::overwrite_existing);
                FlowIO::WriterXML writer;
                writer.SaveSimulationState(outFile, model, globState);
            } catch(std::filesystem::filesystem_error& e) {
                Log::console_error("Warning: Could not create file: %s\n", e.what());
            }

            // if there is a next des limit, handle that
            if(!Settings::desLimit.empty()) {
                model->otfParams.desorptionLimit = Settings::desLimit.front();
                Settings::desLimit.pop_front();
                simManager.ForwardOtfParams(&model->otfParams);
                endCondition = false;
                Log::console_msg_master(1, " Handling next des limit %z\n", model->otfParams.desorptionLimit);

                try {
                    ProcessSleep(1000);
                    simManager.StartSimulation();
                }
                catch (std::runtime_error& e) {
                    Log::console_error("ERROR: Starting simulation: %s\n", e.what());
                    endCondition = true;
                }
            }
        }
        else if(Settings::autoSaveDuration && (uint64_t)(elapsedTime)%Settings::autoSaveDuration==0){ // autosave every x seconds
            Log::console_msg_master(1,"[%.0lfs] Creating auto save file %s\n", elapsedTime, autoSave.c_str());
            FlowIO::WriterXML writer;
            writer.SaveSimulationState(autoSave, model, globState);
        }

        if(Settings::outputDuration && (uint64_t)(elapsedTime)%Settings::outputDuration==0){ // autosave every x seconds
            if((uint64_t)elapsedTime / Settings::outputDuration <= 1){
                Log::console_msg_master(1, "\n%-6s %-14s %-20s %-20s %-20s %-20s %-20s %-20s\n",
                                        "Node#", "Time",
                                        "#Hits (run)", "#Hits (total)","Hit/sec",
                                        "#Des (run)", "#Des (total)","Des/sec");
                Log::console_msg_master(1, "%s\n",std::string(6+14+20+20+20+20+20+20,'-').c_str());
            }
            Log::console_msg(1,"%-6d %-14.2lf %-20zu %-20zu %-20.2lf %-20zu %-20zu %-20.2lf\n",
                             MFMPI::world_rank, elapsedTime,
                             globState.globalHits.globalHits.nbMCHit - oldHitsNb, globState.globalHits.globalHits.nbMCHit,
                             (double) (globState.globalHits.globalHits.nbMCHit - oldHitsNb) /
                             (elapsedTime),
                             globState.globalHits.globalHits.nbDesorbed - oldDesNb, globState.globalHits.globalHits.nbDesorbed,
                             (double) (globState.globalHits.globalHits.nbDesorbed - oldDesNb) /
                             (elapsedTime));
        }

        // Check for potential time end
        if(Settings::simDuration > 0) {
            endCondition |= (elapsedTime >= (double) Settings::simDuration);
        }
    } while(!endCondition);
    simTimer.Stop();
    elapsedTime = simTimer.Elapsed();

    // Terminate simulation
    simManager.StopSimulation();
    simManager.KillAllSimUnits();
    Log::console_msg(1,"[%d][%s] Simulation finished!\n", MFMPI::world_rank, Util::getTimepointString().c_str());

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    Log::console_msg_master(1, "\n%-6s %-14s %-20s %-20s %-20s %-20s %-20s %-20s\n",
                            "Node#", "Time",
                            "#Hits (run)", "#Hits (total)","Hit/sec",
                            "#Des (run)", "#Des (total)","Des/sec");
    Log::console_msg_master(1, "%s\n",std::string(6+14+20+20+20+20+20+20,'-').c_str());
#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    //TODO: Send output to master node
    if(elapsedTime > 1e-4) {
        // Global result print --> TODO: ()
        Log::console_msg(1,"%-6d %-14.2lf %-20zu %-20zu %-20.2lf %-20zu %-20zu %-20.2lf\n",
               MFMPI::world_rank, elapsedTime,
               globState.globalHits.globalHits.nbMCHit - oldHitsNb, globState.globalHits.globalHits.nbMCHit,
               (double) (globState.globalHits.globalHits.nbMCHit - oldHitsNb) /
               (elapsedTime),
               globState.globalHits.globalHits.nbDesorbed - oldDesNb, globState.globalHits.globalHits.nbDesorbed,
                         (double) (globState.globalHits.globalHits.nbDesorbed - oldDesNb) /
                         (elapsedTime));
    }

#if defined(USE_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
    MFMPI::mpi_receive_states(model, globState);
    if(MFMPI::world_rank != 0){
        // Cleanup all files from nodes tmp path
        if (SettingsIO::outputPath.find("tmp") != std::string::npos) {
            std::filesystem::remove_all(SettingsIO::outputPath);
        }
        if(std::filesystem::exists(autoSave)){
            std::filesystem::remove(autoSave);
        }
    }
    // Finalize the MPI environment.
    MPI_Finalize();
    if(MFMPI::world_rank != 0){
        return 0;
    }
    Log::console_msg_master(1, "%s\n",std::string(6+14+20+20+20+20+20+20,'=').c_str());
#endif //USE_MPI

    if(elapsedTime > 1e-4) {
        Log::console_msg(1,"%-6s %-14.2lf %-20zu %-20zu %-20.2lf %-20zu %-20zu %-20.2lf\n",
                         "x", elapsedTime,
                         globState.globalHits.globalHits.nbMCHit - oldHitsNb, globState.globalHits.globalHits.nbMCHit,
                         (double) (globState.globalHits.globalHits.nbMCHit - oldHitsNb) /
                         (elapsedTime),
                         globState.globalHits.globalHits.nbDesorbed - oldDesNb, globState.globalHits.globalHits.nbDesorbed,
                         (double) (globState.globalHits.globalHits.nbDesorbed - oldDesNb) /
                         (elapsedTime));
    }

    if(MFMPI::world_rank == 0){
        // Export results
        //  a) Use existing autosave as base
        //  b) Create copy of input file
        // update geometry info (in case of param sweep)
        // and simply update simulation results
        bool createZip = std::filesystem::path(SettingsIO::outputFile).extension() == ".zip";
        SettingsIO::outputFile = std::filesystem::path(SettingsIO::outputFile).replace_extension(".xml").string();

        std::string fullOutFile = std::filesystem::path(SettingsIO::outputPath).append(SettingsIO::outputFile).string();
        if(std::filesystem::exists(autoSave)){
            std::filesystem::rename(autoSave, fullOutFile);
        }
        else if(!SettingsIO::overwrite){
            // Copy full file description first, in case outputFile is different
            if(!SettingsIO::workFile.empty() && std::filesystem::exists(SettingsIO::workFile)){
                try {
                    std::filesystem::copy_file(SettingsIO::workFile, fullOutFile,
                                               std::filesystem::copy_options::overwrite_existing);
                } catch (std::filesystem::filesystem_error &e) {
                    Log::console_error("Could not copy file to preserve initial file layout: %s\n", e.what());
                }
            }
        }
        FlowIO::WriterXML writer;
        pugi::xml_document newDoc;
        newDoc.load_file(fullOutFile.c_str());
        writer.SaveGeometry(newDoc, model, false, true);
        writer.SaveSimulationState(fullOutFile, model, globState);

        if(SettingsIO::isArchive){
            Log::console_msg_master(3, "Compressing xml to zip...\n");

            //Zipper library
            std::string fileNameWithZIP = std::filesystem::path(fullOutFile).replace_extension(".zip").string();
            if (std::filesystem::exists(fileNameWithZIP)) { // should be workFile == inputFile
                try {
                    std::filesystem::remove(fileNameWithZIP);
                }
                catch (std::exception &e) {
                    Log::console_error("Error compressing to \n{}\nMaybe file is in use:\n{}",fileNameWithZIP, e.what());
                }
            }
            ZipFile::AddFile(fileNameWithZIP, fullOutFile, FileUtils::GetFilename(fullOutFile));
            //At this point, if no error was thrown, the compression is successful
            try {
                std::filesystem::remove(fullOutFile);
            }
            catch (std::exception &e) {
                Log::console_error("Error removing\n{}\nMaybe file is in use:\n{}",fullOutFile,e.what());
            }
        }
    }

    // Cleanup
    SettingsIO::cleanup_files();

    return 0;
}