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


#include <thread>
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
#include <process.h>
#elif not(defined(__MACOSX__) || defined(__APPLE__))
#include <cstring> //memset on unix
#endif

#include "SimulationManager.h"
//#include "Buffer_shared.h" // TODO: Move SHCONTROL to seperate file or SMP.h
#include "SMP.h"
#include "ProcessControl.h"

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cereal/archives/binary.hpp>
#include <../src/Simulation/Simulation.h>
#include <Helper/ConsoleLogger.h>

SimulationManager::SimulationManager(int pid) {
    simulationChanged = true; // by default, always init simulation process the first time
    interactiveMode = true;
    isRunning = false;
    hasErrorStatus = false;
    allProcsDone = false;

    useCPU = false;
    nbThreads = 0;

    useGPU = false;

    useRemote = false;

    if(pid > -1)
        mainProcId = pid;
    else {
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
        mainProcId = _getpid();
#else
        mainProcId = ::getpid();
#endif
    }
}

SimulationManager::~SimulationManager() {
    KillAllSimUnits();
    for(auto& unit : simUnits){
        delete unit;
    }
}

//! Refresh proc status by looking for those that can be safely killed and remove them
int SimulationManager::refreshProcStatus() {
    int nbDead = 0;
    for(auto proc = simHandles.begin(); proc != simHandles.end() ; ){
        if(!(*proc).first.joinable()){
            auto myHandle = (*proc).first.native_handle();
#if defined(_WIN32) && defined(_MSC_VER)
            TerminateThread(myHandle, 1);
#else
            //Linux
            pthread_cancel(myHandle);
#endif
            proc = this->simHandles.erase(proc);
            ++nbDead;
        }
        else{
            ++proc;
        }
    }
    return nbDead;
}

/*!
 * @brief Starts the simulation on all available simulation units
 * @return 0=start successful, 1=PROCESS_DONE state entered
 */
int SimulationManager::StartSimulation() {

    LoadSimulation();

    if(interactiveMode) {
        refreshProcStatus();
        if (simHandles.empty())
            throw std::logic_error("No active simulation handles!");

        /*if(simulationChanged){
            if (ExecuteAndWait(COMMAND_LOAD, PROCESS_READY, 0, 0)) {
                throw std::runtime_error(MakeSubProcError("Subprocesses could not start the simulation"));
            }
            simulationChanged = false;
        }*/
        if (ExecuteAndWait(COMMAND_START, PROCESS_RUN, 0, 0)) {
            throw std::runtime_error(MakeSubProcError("Subprocesses could not start the simulation"));
        }
    }
    else {
        /*if(simulationChanged){
            this->procInformation.masterCmd  = COMMAND_LOAD; // TODO: currently needed to not break the loop
            for(auto& con : simController){
                con.Load();
            }
            simulationChanged = false;
        }*/
        this->procInformation.masterCmd = COMMAND_START; // TODO: currently needed to not break the loop
        for(auto& con : simController){
            con.Start();
        }
    }

    if(allProcsDone){
        isRunning = false;
        return 1;
    }
    isRunning = true;
    return 0;
}

//! Call simulation controllers to stop running simulations
int SimulationManager::StopSimulation() {
    isRunning = false;
    if(interactiveMode) {
        refreshProcStatus();
        if (simHandles.empty())
            return 1;

        if (ExecuteAndWait(COMMAND_PAUSE, PROCESS_READY, 0, 0))
            throw std::runtime_error(MakeSubProcError("Subprocesses could not stop the simulation"));
    }
    else {
        /*if (ExecuteAndWait(COMMAND_PAUSE, PROCESS_READY, 0, 0))
            throw std::runtime_error(MakeSubProcError("Subprocesses could not stop the simulation"));
        */
        return 1;
    }

    return 0;
}

int SimulationManager::LoadSimulation(){
    if(interactiveMode) {
        if (simHandles.empty())
            throw std::logic_error("No active simulation handles!");

        if (ExecuteAndWait(COMMAND_LOAD, PROCESS_READY, 0, 0)) {
            std::string errString = "Failed to send geometry to sub process:\n";
            errString.append(GetErrorDetails());
            throw std::runtime_error(errString);
        }
        simulationChanged = false;
    }
    else{
        bool errorOnLoad = false;
        for(auto& con : simController){
            errorOnLoad |= con.Load();
        }
        if(errorOnLoad){
            std::cerr << "Failed to load simulation!" << std::endl;
            return 1;
        }
        simulationChanged = false;
    }

    return 0;
}

//! Create simulation handles for CPU simulations with n Threads
int SimulationManager::CreateCPUHandle() {
    uint32_t processId = mainProcId;

    //Get number of cores
    if(nbThreads == 0) {
#if defined(DEBUG)
        nbThreads = 1;
#else
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
        SYSTEM_INFO sysinfo;
        GetSystemInfo(&sysinfo);
        nbThreads = (size_t) sysinfo.dwNumberOfProcessors;
#else
        nbThreads = (unsigned int) sysconf(_SC_NPROCESSORS_ONLN);
#endif // WIN
#endif // DEBUG
    }

    if(!simUnits.empty()){
        for(auto& sim : simUnits){
            delete sim;
        }
    }
    size_t nbSimUnits = 1; // nbThreads
    try{
        simUnits.resize(nbSimUnits);
        procInformation.Resize(nbThreads);
        simController.clear();
        simHandles.clear();
    }
    catch (const std::exception &e){
        std::cerr << "[SimManager] Invalid resize/clear " << nbThreads<< std::endl;
        throw std::runtime_error(e.what());
    }

    for(size_t t = 0; t < nbSimUnits; ++t){
        simUnits[t] = new Simulation();
    }
    simController.emplace_back(SimulationController{processId, 0, nbThreads,
                                                    simUnits.back(), &procInformation});
    if(interactiveMode) {
        simHandles.emplace_back(
                /*StartProc(arguments, STARTPROC_NOWIN),*/
                std::thread(&SimulationController::controlledLoop, &simController[0], NULL, nullptr),
                SimType::simCPU);
        /*simUnits.emplace_back(Simulation{nbThreads});
        procInformation.emplace_back(SubDProcInfo{});
        simController.emplace_back(SimulationController{"molflow", processId, iProc, nbThreads, &simUnits.back(), &procInformation.back()});
        simHandles.emplace_back(
                *//*StartProc(arguments, STARTPROC_NOWIN),*//*
            std::thread(&SimulationController::controlledLoop,&simController.back(),NULL,nullptr),
            SimType::simCPU);*/
        auto myHandle = simHandles.back().first.native_handle();
#if defined(_WIN32) && defined(_MSC_VER)
        SetThreadPriority(myHandle, THREAD_PRIORITY_IDLE);
#else
        int policy;
        struct sched_param param{};
        pthread_getschedparam(myHandle, &policy, &param);
        param.sched_priority = sched_get_priority_min(policy);
        pthread_setschedparam(myHandle, policy, &param);
        //Check! Some documentation says it's always 0
#endif
    }

    return 0;
}

// return 1=error
int SimulationManager::CreateGPUHandle() {
    return 1;
}

// return 1=error
int SimulationManager::CreateRemoteHandle() {
    return 1;
}

/*!
 * @brief Creates Simulation Units and waits for their ready status
 * @return 0=all SimUnits are ready, else = ret Units are active, but not all could be launched
 */
int SimulationManager::InitSimUnits() {

    if(useCPU){
        // Launch nbCores subprocesses
        CreateCPUHandle();
    }
    if(useGPU){
        CreateGPUHandle();
        //procInformation.push_back(simController.back().procInfo);
    }
    if(useRemote){
        CreateRemoteHandle();
        //procInformation.push_back(simController.back().procInfo);
    }

    return WaitForProcStatus(PROCESS_READY);
}

/*!
 * @brief Creates Simulation Units and waits for their ready status
 * @return 0=all SimUnits are ready, else = ret Units are active, but not all could be launched
 */
int SimulationManager::InitSimulation(const std::shared_ptr<SimulationModel>& model, GlobalSimuState *globState) {
    if (!model->m.try_lock()) {
        return 1;
    }
    // Prepare simulation unit
    ResetSimulations();
    ForwardSimModel(model);
    ForwardGlobalCounter(globState, nullptr);

    //bool invalidLoad = LoadSimulation();
    model->m.unlock();

    /*if(invalidLoad){
        std::string errString = "Failed to send geometry to sub process:\n";
        errString.append(GetErrorDetails());
        throw std::runtime_error(errString);
        return 1;
    }*/

    return 0;
}

/*!
 * @brief Wait until all SimulationUnits are in procStatus or reach another endstate (error, done)
 * @param procStatus Process Status that should be waited for
 * @return 0 if wait is successful
 */
int SimulationManager::WaitForProcStatus(const uint8_t procStatus) {
    // Wait for completion
    bool finished = false;
    const int waitAmount = 250;
    int prevIncTime = 0; // save last time a wait increment has been set; allows for a dynamic/reasonable increase
    int waitTime = 0;
    int timeOutAt = 40000; // 10 sec; max time for an idle operation to timeout
    allProcsDone = true;
    hasErrorStatus = false;

    // struct, because vector of char arrays is forbidden w/ clang
    struct StateString {
        char s[128];
    };
    std::vector<StateString> prevStateStrings(procInformation.subProcInfo.size());
    std::vector<StateString> stateStrings(procInformation.subProcInfo.size());

    {
        for (size_t i = 0; i < procInformation.subProcInfo.size(); i++) {
            snprintf(prevStateStrings[i].s, 128, "%s", procInformation.subProcInfo[i].statusString);
        }
    }

    do {

        finished = true;

        for (size_t i = 0; i < procInformation.subProcInfo.size(); i++) {
            auto procState = procInformation.subProcInfo[i].slaveState;
            if(procStatus == PROCESS_KILLED) // explicitly ask for killed state
                finished = finished & (procState==PROCESS_KILLED);
            else
                finished = finished & (procState==procStatus || procState==PROCESS_ERROR || procState==PROCESS_DONE);
            if( procState==PROCESS_ERROR ) {
                hasErrorStatus = true;
            }
            else if(procState == PROCESS_STARTING){
                snprintf(stateStrings[i].s, 128, "%s", procInformation.subProcInfo[i].statusString);
                if(strcmp(prevStateStrings[i].s, stateStrings[i].s) != 0) { // if strings are different
                    timeOutAt += (waitTime + 10000 < timeOutAt) ? (waitTime - prevIncTime) : (timeOutAt - waitTime +
                                                                                10000); // if task properly started, increase allowed wait time
                    prevIncTime = waitTime;
                }
            }
            allProcsDone = allProcsDone & (procState == PROCESS_DONE);
        }

        if (!finished) {
            ProcessSleep(waitAmount);
            waitTime += waitAmount;
        }
    } while (!finished && waitTime<timeOutAt);

    return waitTime>=timeOutAt || hasErrorStatus; // 0 = finished, 1 = timeout
}

//! Forward a command to simulation controllers
int SimulationManager::ForwardCommand(const int command, const size_t param, const size_t param2) {
    // Send command

    procInformation.masterCmd = command;
    procInformation.cmdParam = param;
    procInformation.cmdParam2 = param2;
    for(auto & i : procInformation.subProcInfo) {
        auto procState = i.slaveState;
        if(procState == PROCESS_READY || procState == PROCESS_RUN || procState==PROCESS_ERROR || procState==PROCESS_DONE) { // check if it'' ready before sending a new command
            i.slaveState = PROCESS_STARTING;
        }
    }

    return 0;
}

/*!
 * @brief Shortcut function combining ForwardCommand() and WaitForProcStatus() into a single call
 * @param command execution command for every subprocess
 * @param procStatus status that every subprocess has to reach
 * @param param additional command parameter
 * @return 0=success, 1=fail
 */
int SimulationManager::ExecuteAndWait(const int command, const uint8_t procStatus, const size_t param,
                                      const size_t param2) {
    if(!ForwardCommand(command, param, param2)) { // execute
        if (!WaitForProcStatus(procStatus)) { // and wait
            return 0;
        }
        return 1;
    }
    return 2;
}

int SimulationManager::KillAllSimUnits() {
    if (!simHandles.empty()) {
        if (ExecuteAndWait(COMMAND_EXIT, PROCESS_KILLED)) { // execute
            // Force kill
            for (auto& con : simController)
                con.EmergencyExit();
            if (ExecuteAndWait(COMMAND_EXIT, PROCESS_KILLED)) {
                int i = 0;
                for (auto tIter = simHandles.begin(); tIter != simHandles.end(); ++i) {
                    if (procInformation.subProcInfo[i].slaveState != PROCESS_KILLED) {
                        auto nativeHandle = simHandles[i].first.native_handle();
#if defined(_WIN32) && defined(_MSC_VER)
                        //Windows
                        //TerminateThread(nativeHandle, 1);
#else
                        //Linux
                        int s;
                        s = pthread_cancel(nativeHandle);
                        if (s != 0)
                            Log::console_msg(1, "pthread_cancel: {}\n", s);
                        tIter->first.detach();
#endif
                        //assume that the process doesn't exist, so remove it from our management structure
                        try {
                            tIter = simHandles.erase(tIter);
                        }
                        catch (const std::exception& e) {
                            char tmp[512];
                            snprintf(tmp, 512, "Could not terminate sub processes: %s\n", e.what());
                            throw std::runtime_error(tmp); // proc couldn't be killed!?
                        }
                        catch (...) {
                            char tmp[512];
                            snprintf(tmp, 512, "Could not terminate sub processes\n");
                            throw std::runtime_error(tmp); // proc couldn't be killed!?
                        }
                    }
                }
            }
        }

        // Join all threads before clearing the list
        int i = 0;
        for (auto& t : simHandles) {
            if (procInformation.subProcInfo[i++].slaveState == PROCESS_KILLED) {
                t.first.join();
            }
            else {
                if (ExecuteAndWait(COMMAND_EXIT, PROCESS_KILLED))
                    exit(1);
                auto nativeHandle = t.first.native_handle();
#if defined(_WIN32) && defined(_MSC_VER)
                //Windows
                //TerminateThread(nativeHandle, 1);
#else
                //Linux
                //pthread_cancel(nativeHandle); // TODO: Causes crash?
#endif
            }
        }
        simHandles.clear();
    }
    nbThreads = 0;
    return 0;
}


int SimulationManager::ResetSimulations() {
    if(interactiveMode) {
        if (ExecuteAndWait(COMMAND_CLOSE, PROCESS_READY, 0, 0))
            throw std::runtime_error(MakeSubProcError("Subprocesses could not restart"));
    }
    else {
        for(auto& con : simController){
            con.Reset();
        }
    }
    return 0;
}

int SimulationManager::ResetHits() {
    isRunning = false;
    if(interactiveMode) {
        if (ExecuteAndWait(COMMAND_RESET, PROCESS_READY, 0, 0))
            throw std::runtime_error(MakeSubProcError("Subprocesses could not reset hits"));
    }
    else {
        for(auto& con : simController){
            con.Reset();
        }
    }
    return 0;
}

int SimulationManager::GetProcStatus(ProcComm &procInfoList) {
    if(simHandles.empty())
        return 1;

    procInfoList.subProcInfo = procInformation.subProcInfo;

    return 0;
}

int SimulationManager::GetProcStatus(size_t *states, std::vector<std::string>& statusStrings) {

    if(statusStrings.size() < procInformation.subProcInfo.size())
        return 1;

    for (size_t i = 0; i < procInformation.subProcInfo.size(); i++) {
        //states[i] = shMaster->procInformation[i].masterCmd;
        states[i] = procInformation.subProcInfo[i].slaveState;
        char tmp[128];
        strncpy(tmp, procInformation.subProcInfo[i].statusString, 127);
        tmp[127] = 0;
        statusStrings[i] = tmp;
    }
    return 0;
}

bool SimulationManager::GetLockedHitBuffer() {
    return true;
}

int SimulationManager::UnlockHitBuffer() {
    return 0;
}

/*!
 * @brief Actively check running state (evaluate done & error)
 * @return 1 if simulation processes are running, 0 else
 */
bool SimulationManager::GetRunningStatus(){

    bool done = true;
    for (auto & i : procInformation.subProcInfo) {
        auto procState = i.slaveState;
        done = done & (procState==PROCESS_ERROR || procState==PROCESS_DONE);
        if( procState==PROCESS_ERROR ) {
            hasErrorStatus = true;
        }
    }

    allProcsDone = done;
    isRunning = isRunning && !allProcsDone;

    return isRunning;
}

/*!
 * @brief Return error information or current running state in case of a hangup
 * @return char array containing proc status (and error message/s)
 */
std::string SimulationManager::GetErrorDetails() {

    ProcComm procInfo;
    GetProcStatus(procInfo);

    std::string err;
    //static char err[1024];
    //std::memset(err, 0, 1024);

    for (size_t i = 0; i < procInfo.subProcInfo.size(); i++) {
        char tmp[512]{'\0'};
        size_t state = procInfo.subProcInfo[i].slaveState;
        if (state == PROCESS_ERROR) {
            sprintf(tmp, "[Thread #%zd] %s: %s\n", i, prStates[state],
                    procInfo.subProcInfo[i].statusString);
        } else {
            sprintf(tmp, "[Thread #%zd] %s\n", i, prStates[state]);
        }
        // Append at most up to character 1024, strncat would otherwise overwrite in memory
        err.append(tmp);
        if(err.size() >= 1024) {
            err.resize(1024);
            break;
        }
        /*strncat(err, tmp, std::min(1024 - strlen(err), (size_t)512));
        if(strlen(err) >= 1024) {
            err[1023] = '\0';
            break;
        }*/
    }
    return err;
}

std::string SimulationManager::MakeSubProcError(const char *message) {
    std::string errString;
    if (!message){
        errString.append("Bad response from sub process(es):\n");
        errString.append(GetErrorDetails());
    }
    else{
        errString.append(message);
        errString.append(":\n");
        errString.append(GetErrorDetails());
    }
    return errString;
}


int SimulationManager::ShareWithSimUnits(void *data, size_t size, LoadType loadType) {
    if(loadType < LoadType::NLOADERTYPES){
        /*if(CreateLoaderDP(size))
            return 1;
        if(UploadToLoader(data, size))
            return 1;*/
    }


    switch (loadType) {
        case LoadType::LOADGEOM:{
            if (ExecuteAndWait(COMMAND_LOAD, PROCESS_READY, size, 0)) {
                //CloseLoaderDP();
                std::string errString = "Failed to send geometry to sub process:\n";
                errString.append(GetErrorDetails());
                throw std::runtime_error(errString);
            }
            //CloseLoaderDP();
            break;
        }
        case LoadType::LOADPARAM:{

            //if (ExecuteAndWait(COMMAND_UPDATEPARAMS, isRunning ? PROCESS_RUN : PROCESS_READY, size, isRunning ? PROCESS_RUN : PROCESS_READY)) {
            if (ExecuteAndWait(COMMAND_UPDATEPARAMS, PROCESS_READY, size, isRunning ? PROCESS_RUN : PROCESS_READY)) {
                //CloseLoaderDP();
                std::string errString = "Failed to send params to sub process:\n";
                errString.append(GetErrorDetails());
                throw std::runtime_error(errString);
            }
            //CloseLoaderDP();
            if(isRunning) { // restart
                StartSimulation();
            }
            break;
        }
        case LoadType::LOADHITS:{
        }
        default:{
            // Unspecified load type
            return 1;
        }
    }

    return 0;
}

void SimulationManager::ForwardGlobalCounter(GlobalSimuState *simState, ParticleLog *particleLog) {
    for(auto& simUnit : simUnits) {
        if(!simUnit->tMutex.try_lock_for(std::chrono::seconds(10)))
            return;
        simUnit->globState = simState;
        simUnit->globParticleLog = particleLog;
        simUnit->tMutex.unlock();
    }
}

// Create hard copy for local usage
void SimulationManager::ForwardSimModel(const std::shared_ptr<SimulationModel>& model) {
    for(auto& sim : simUnits)
        sim->model = model;
}

// Create hard copy for local usage and resie particle logger
void SimulationManager::ForwardOtfParams(OntheflySimulationParams *otfParams) {
    for(auto& sim : simUnits) {
        if(sim->model)
            sim->model->otfParams = *otfParams;
        sim->ReinitializeParticleLog();
        sim->ReinitializeConvergenceLog();
    }

}

/**
* \brief Saves current facet hit counter from cache to results
*/
void SimulationManager::ForwardFacetHitCounts(std::vector<FacetHitBuffer*>& hitCaches) {
    for(auto& simUnit : simUnits){
        if(simUnit->globState->facetStates.size() != hitCaches.size()) return;
        if(!simUnit->tMutex.try_lock_for(std::chrono::seconds(10)))
            return;
        for (size_t i = 0; i < hitCaches.size(); i++) {
            simUnit->globState->facetStates[i].momentResults[0].hits = *hitCaches[i];
        }
        simUnit->tMutex.unlock();
    }
}

int SimulationManager::IncreasePriority() {
#if defined(_WIN32) && defined(_MSC_VER)
    // https://cpp.hotexamples.com/de/examples/-/-/SetPriorityClass/cpp-setpriorityclass-function-examples.html
    HANDLE hProcess = GetCurrentProcess();
    SetPriorityClass(hProcess, HIGH_PRIORITY_CLASS);
#else

#endif

    for(auto& handle : simHandles) {
        auto myHandle = handle.first.native_handle();
#if defined(_WIN32) && defined(_MSC_VER)
        SetThreadPriority(myHandle, THREAD_PRIORITY_HIGHEST);
#else
        int policy;
        struct sched_param param{};
        pthread_getschedparam(myHandle, &policy, &param);
        param.sched_priority = sched_get_priority_min(policy);
        pthread_setschedparam(myHandle, policy, &param);
        //Check! Some documentation says it's always 0
#endif
    }

    return 0;
}

int SimulationManager::DecreasePriority() {
#if defined(_WIN32) && defined(_MSC_VER)
    HANDLE hProcess = GetCurrentProcess();
    SetPriorityClass(hProcess, IDLE_PRIORITY_CLASS);
#else

#endif
    for(auto& handle : simHandles) {
        auto myHandle = handle.first.native_handle();
#if defined(_WIN32) && defined(_MSC_VER)
        SetThreadPriority(myHandle, THREAD_PRIORITY_IDLE);
#else
        int policy;
            struct sched_param param{};
            pthread_getschedparam(myHandle, &policy, &param);
            param.sched_priority = sched_get_priority_min(policy);
            pthread_setschedparam(myHandle, policy, &param);
            //Check! Some documentation says it's always 0
#endif
    }
    return 0;
}

int SimulationManager::RefreshRNGSeed(bool fixed) {
    if(simUnits.empty() || nbThreads == 0)
        return 1;
    for(auto& sim : simUnits){
        sim->SetNParticle(nbThreads, fixed);
    }

    return 0;
}

/*!
 * @brief Calls for a rebuild of the acceleration structure on the simulation side
 * @return 0=start successful, 1=PROCESS_DONE state entered
 */
int SimulationManager::BuildAccelStructure() {
    if(interactiveMode) {
        refreshProcStatus();
        if (simHandles.empty())
            throw std::logic_error("No active simulation handles!");

        if (ExecuteAndWait(COMMAND_ACCEL, PROCESS_RUN, 0, 0)) {
            throw std::runtime_error(MakeSubProcError("Subprocess could not build acceleration data structure!"));
        }
    }
    else {
        this->procInformation.masterCmd = COMMAND_ACCEL; // TODO: currently needed to not break the loop
        for(auto& con : simController){
            con.RebuildAccel();
        }
    }

    return 0;
}