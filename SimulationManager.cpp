//
// Created by pbahr on 15/04/2020.
//


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

SimulationManager::SimulationManager(const std::string &appName , const std::string &dpName) {
    isRunning = false;
    hasErrorStatus = false;
    allProcsDone = false;

    useCPU = false;
    nbCores = 0;
    nbThreads = 0;

    useGPU = false;

    useRemote = false;

    dpLog = nullptr;

#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    uint32_t pid = _getpid();
    const char* dpPrefix = "";
#else
    uint32_t pid = ::getpid();
    const char *dpPrefix = "/"; // creates semaphore as /dev/sem/%s_sema
#endif

    sprintf(this->appName,"%s", appName.c_str());
    sprintf(this->logDpName,"%s", std::string(dpPrefix+dpName+"LOG"+std::to_string(pid)).c_str());

}

SimulationManager::~SimulationManager() {
    CLOSEDP(dpLog);
    KillAllSimUnits();
    for(int t = 0; t < nbThreads; ++t){
        delete simUnits[t];
        simUnits[t] = nullptr;
    }

    /*for(auto& handle : simHandles){

        handle.first.join();
    }*/
}

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

int SimulationManager::LoadInput(const std::string& fileName) {
    std::ifstream inputFile(fileName);
    inputFile.seekg(0, std::ios::end);
    size_t size = inputFile.tellg();
    std::string buffer(size, ' ');
    inputFile.seekg(0);
    inputFile.read(&buffer[0], size);

    try {
        ShareWithSimUnits((BYTE *) buffer.c_str(), buffer.size(), LoadType::LOADGEOM);
    }
    catch (std::runtime_error& e) {
        throw e;
    }
    return 0;
}

int SimulationManager::ResetStatsAndHits() {

    return 0;
}

int SimulationManager::ReloadLogBuffer(size_t logSize, bool ignoreSubs) {//Send simulation mode changes to subprocesses without reloading the whole geometry
    if (simHandles.empty())
        throw std::logic_error("No active simulation handles!");

    if(!ignoreSubs) {
        //if (ExecuteAndWait(COMMAND_RELEASEDPLOG, isRunning ? PROCESS_RUN : PROCESS_READY,isRunning ? PROCESS_RUN : PROCESS_READY)) {
        if (ExecuteAndWait(COMMAND_RELEASEDPLOG, PROCESS_READY,isRunning ? PROCESS_RUN : PROCESS_READY)) {
            throw std::runtime_error(MakeSubProcError("Subprocesses didn't release dpLog handle"));
        }
    }
    //To do: only close if parameters changed
    if (dpLog && logSize == dpLog->size) {
        ClearLogBuffer(); //Fills values with 0
    } else {
        CloseLogDP();
        if (logSize)
            if(CreateLogDP(logSize)) {
                throw std::runtime_error(
                        "Failed to create 'dpLog' dataport.\nMost probably out of memory.\nReduce number of logged particles in Particle Logger.");
            }
    }

    return 0;
}

int SimulationManager::ReloadHitBuffer(size_t hitSize) {
    return ResetHits();
}

/*!
 * @brief Starts the simulation on all available simulation units
 * @return 0=start successful, 1=PROCESS_DONE state entered
 */
int SimulationManager::StartSimulation() {
    refreshProcStatus();
    if (simHandles.empty())
        throw std::logic_error("No active simulation handles!");

    if (ExecuteAndWait(COMMAND_START, PROCESS_RUN, 0, 0)){ // TODO: 0=MC_MODE, AC_MODE should be seperated completely
        throw std::runtime_error(MakeSubProcError("Subprocesses could not start the simulation"));
    }

    if(allProcsDone){
        isRunning = false;
        return 1;
    }
    isRunning = true;
    return 0;
}

int SimulationManager::StopSimulation() {
    isRunning = false;
    refreshProcStatus();
    if (simHandles.empty())
        return 1;

    if (ExecuteAndWait(COMMAND_PAUSE, PROCESS_READY, 0, 0))
        throw std::runtime_error(MakeSubProcError("Subprocesses could not stop the simulation"));
    return 0;
}

/*!
 * @brief Convenience function that stops a running simulation and starts a paused simulation
 * @return 0=success, 1=else
 * @todo add benchmark
 */
bool SimulationManager::StartStopSimulation(){
    if (isRunning) {

        // Stop
        //InnerStop(appTime);
        try {
            StopSimulation();
            //Update(appTime);
        }

        catch (std::runtime_error &e) {
            throw e;
        }
    } else {

        // Start
        try {
            //if (needsReload) RealReload(); //Synchronize subprocesses to main process

            StartSimulation();
        }
        catch (std::runtime_error &e) {
            throw e;
        }

        // Particular case when simulation ends before getting RUN state
        if (allProcsDone) {
            isRunning = false;
            //Update(appTime);
            //GLMessageBox::Display("Max desorption reached", "Information (Start)", GLDLG_OK, GLDLG_ICONINFO);
        }

    }
    return isRunning; // return previous state
}


int SimulationManager::TerminateSimHandles() {

    return 0;
}

int SimulationManager::CreateCPUHandle(uint16_t iProc) {
    uint32_t processId;

#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    processId = _getpid();
#else
    processId = ::getpid();
#endif //  WIN

    //Get number of cores
    if(nbThreads == 0) {
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
        SYSTEM_INFO sysinfo;
        GetSystemInfo(&sysinfo);
        nbThreads = (size_t) sysinfo.dwNumberOfProcessors;
#else
        nbThreads = (unsigned int) sysconf(_SC_NPROCESSORS_ONLN);
#endif
    }
    if(!simUnits.empty()){
        for(auto& sim : simUnits){
            delete sim;
        }
    }
    try{
        simUnits.resize(nbThreads);
        procInformation.Resize(nbThreads);
        simController.clear();
        simHandles.clear();
    }
    catch (std::exception& e){
        std::cerr << "[SimManager] Invalid resize/clear "<<iProc<<" / " << nbThreads<< std::endl;
        throw std::runtime_error(e.what());
    }

    for(int t = 0; t < nbThreads; ++t){
        simUnits[t] = new Simulation();
    }
    simController.emplace_back(SimulationController{"molflow", processId, iProc++, nbThreads,
                                                    reinterpret_cast<std::vector<SimulationUnit *> *>(&simUnits), &procInformation});
    simHandles.emplace_back(
            /*StartProc(arguments, STARTPROC_NOWIN),*/
            std::thread(&SimulationController::controlledLoop,&simController[0],NULL,nullptr),
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

//TODO: This is Molflow only
int SimulationManager::CreateLogDP(size_t logDpSize) {
    //size_t logDpSize = sizeof(size_t) + ontheflyParams.logLimit * sizeof(ParticleLoggerItem);
    dpLog = CreateDataport(logDpName, logDpSize);
    if (!dpLog) {
        //progressDlg->SetVisible(false);
        //SAFE_DELETE(progressDlg);
        return 1;
        //throw Error("Failed to create 'dpLog' dataport.\nMost probably out of memory.\nReduce number of logged particles in Particle Logger.");
    }

    return 0;
}

/*!
 * @brief Creates Simulation Units and waits for their ready status
 * @return 0=all SimUnits are ready, else = ret Units are active, but not all could be launched
 */
int SimulationManager::InitSimUnits() {

    if(useCPU){
        // Launch nbCores subprocesses
        auto nbActiveProcesses = simHandles.size();
        for(int iProc = 0; iProc < nbCores; ++iProc) {
            if(CreateCPUHandle(iProc + nbActiveProcesses)){ // abort initialization when creation fails
                nbCores = simHandles.size();
                return simHandles.size();
            }
            //procInformation.push_back(simController.back().procInfo);
        }
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
    int timeOutAt = 10000; // 10 sec; max time for an idle operation to timeout
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
    }
    return 1;
}

int SimulationManager::KillAllSimUnits() {
    if( !simHandles.empty() ) {
        if(ExecuteAndWait(COMMAND_EXIT, PROCESS_KILLED)){ // execute
            // Force kill

            for(size_t i=0;i<simHandles.size();i++) {
                if (procInformation.subProcInfo[i].slaveState != PROCESS_KILLED){
                    auto nativeHandle = simHandles[i].first.native_handle();
#if defined(_WIN32) && defined(_MSC_VER)
                    //Windows
				    TerminateThread(nativeHandle, 1);
#else
                    //Linux
                    pthread_cancel(nativeHandle);
#endif
                    //assume that the process doesn't exist, so remove it from our management structure
                    simHandles.erase((simHandles.begin()+i));
                    throw std::runtime_error(MakeSubProcError("Could not terminate sub processes")); // proc couldn't be killed!?

                }
            }
        }

        for(size_t i=0;i<simHandles.size();i++) {
            if (procInformation.subProcInfo[i].slaveState == PROCESS_KILLED) {
                simHandles[i].first.join();
            }
            else{
                auto nativeHandle = simHandles[i].first.native_handle();
#if defined(_WIN32) && defined(_MSC_VER)
                //Windows
                TerminateThread(nativeHandle, 1);
#else
                //Linux
                pthread_cancel(nativeHandle);
#endif
            }
        }
        simHandles.clear();
    }
    nbCores = 0;
    return 0;
}

int SimulationManager::ClearLogBuffer() {
    if (!dpLog) {
        return 1;
    }
    AccessDataport(dpLog);
    memset(dpLog->buff, 0, dpLog->size); //Also clears hits, leaks
    ReleaseDataport(dpLog);
    return 0;
}

int SimulationManager::ResetSimulations() {
    if (ExecuteAndWait(COMMAND_CLOSE, PROCESS_READY, 0, 0))
        throw std::runtime_error(MakeSubProcError("Subprocesses could not restart"));
    return 0;
}

int SimulationManager::ResetHits() {
    isRunning = false;
    if (ExecuteAndWait(COMMAND_RESET, PROCESS_READY, 0, 0))
        throw std::runtime_error(MakeSubProcError("Subprocesses could not reset hits"));
    return 0;
}

int SimulationManager::CloseLogDP() {
    CLOSEDP(dpLog);
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

BYTE *SimulationManager::GetLockedLogBuffer() {
    if (dpLog && AccessDataport(dpLog)) {
        return (BYTE*)dpLog->buff;
    }
    return nullptr;
}

int SimulationManager::UnlockLogBuffer() {
    if (dpLog && ReleaseDataport(dpLog)) {
        return 0;
    }
    return 1;
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
const char *SimulationManager::GetErrorDetails() {

    ProcComm procInfo;
    GetProcStatus(procInfo);

    static char err[1024];
    strcpy(err, "");

    for (size_t i = 0; i < procInfo.subProcInfo.size(); i++) {
        char tmp[512];
        size_t state = procInfo.subProcInfo[i].slaveState;
        if (state == PROCESS_ERROR) {
            sprintf(tmp, "[#%zd] Process [PID %zu] %s: %s\n", i, procInfo.subProcInfo[i].procId, prStates[state],
                    procInfo.subProcInfo[i].statusString);
        } else {
            sprintf(tmp, "[#%zd] Process [PID %zu] %s\n", i, procInfo.subProcInfo[i].procId, prStates[state]);
        }
        strncat(err, tmp, 512);
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
        case LoadType::LOADAC:{
            if (ExecuteAndWait(COMMAND_LOADAC, PROCESS_RUNAC, size, 0)) {
                //CloseLoaderDP();
                std::string errString = "Failed to send AC geometry to sub process:\n";
                errString.append(GetErrorDetails());
                throw std::runtime_error(errString);
            }
            //CloseLoaderDP();
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

void SimulationManager::ForwardGlobalCounter(GlobalSimuState *simState) {
    for(auto& simUnit : simUnits) {
        if(!simUnit->tMutex.try_lock_for(std::chrono::seconds(10)))
            return;
        simUnit->globState = simState;
        simUnit->tMutex.unlock();
    }
}

// Create hard copy for local usage
void SimulationManager::ForwardSimModel(SimulationModel *model) {
    for(auto& sim : simUnits)
        sim->model = *model;
}

// Create hard copy for local usage
void SimulationManager::ForwardOtfParams(OntheflySimulationParams *otfParams) {
    for(auto& sim : simUnits)
        sim->model.otfParams = *otfParams;
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