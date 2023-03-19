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

#include <cmath>
#include <sstream>
#include "SimulationController.h"
#include "../src/Simulation/Simulation.h"
#include "ProcessControl.h"
#include <Helper/ConsoleLogger.h>
#include <Helper/OutputHelper.h>

#include <omp.h>

#if not defined(_MSC_VER)
#include <sys/time.h>
#include <sys/resource.h>
#endif



#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
#include <process.h>
#endif

#define WAITTIME    500  // Answer in STOP mode

void InitOMP(){
    // "Warm up" threads, to remove overhead for performance benchmarks
    double randomCounter = 0;
#pragma omp parallel default(none) shared(randomCounter)
    {
        double local_result = 0;
#pragma omp for
        for (int i=0; i < 1000; i++) {
            local_result += 1;
        }
#pragma omp critical
        randomCounter += local_result;
    }
    DEBUG_PRINT("[OMP] Init: %f\n", randomCounter);
}

SimThread::SimThread(ProcComm *procInfo, SimulationUnit *sim, size_t threadNum) {
    this->threadNum = threadNum;
    this->procInfo = procInfo;
    this->status = nullptr;
    timeLimit = 0.0;
    simulation = sim;
    stepsPerSec = 1.0;
    particle = nullptr;
    simEos = false;
    localDesLimit = 0;
}

SimThread::~SimThread() = default;


// todo: fuse with runSimulation()
// Should allow simulation for N steps opposed to T seconds
int SimThread::advanceForSteps(size_t desorptions) {
    size_t nbStep = (stepsPerSec <= 0.0) ? 250 : std::ceil(stepsPerSec + 0.5);

    double timeStart = omp_get_wtime();
    double timeEnd = timeStart;
    do {
        simEos = !particle->SimulationMCStep(nbStep, threadNum, desorptions);
        timeEnd = omp_get_wtime();

        const double elapsedTimeMs = (timeEnd - timeStart); //std::chrono::duration<float, std::ratio<1, 1>>(end_time - start_time).count();
        if (elapsedTimeMs > 1e-6)
            stepsPerSec = (10.0 * nbStep) / (elapsedTimeMs); // every 1.0 second
        else
            stepsPerSec = (100.0 * nbStep); // in case of fast initial run

        desorptions -= particle->tmpState.globalHits.globalHits.nbDesorbed;
    } while (desorptions);


    return 0;
}

// run until end or until autosaveTime check
int SimThread::advanceForTime(double simDuration) {

    size_t nbStep = (stepsPerSec <= 0.0) ? 250 : std::ceil(stepsPerSec + 0.5);

    double timeStart = omp_get_wtime();
    double timeEnd = timeStart;
    do {
        simEos = !particle->SimulationMCStep(nbStep, threadNum, 0);
        timeEnd = omp_get_wtime();

        const double elapsedTimeMs = (timeEnd - timeStart); //std::chrono::duration<float, std::ratio<1, 1>>(end_time - start_time).count();
        if (elapsedTimeMs > 1e-6)
            stepsPerSec = (10.0 * nbStep) / (elapsedTimeMs); // every 1.0 second
        else
            stepsPerSec = (100.0 * nbStep); // in case of fast initial run

    } while (simDuration > timeEnd - timeStart);


    return 0;
}

/**
* \brief Simulation loop for an individual thread
 * \return 0 when simulation end has been reached via desorption limit, 1 otherwise
 */

bool SimThread::runLoop() {
    enum EoS : int {
        CONT,
        STOP,
        HALT
    };

    bool eos;
    bool lastUpdateOk = false;

    //printf("Lim[%zu] %lu --> %lu\n",threadNum, localDesLimit, simulation->globState->globalHits.globalHits.hit.nbDesorbed);

    double timeStart = omp_get_wtime();
    double timeLoopStart = timeStart;
    double timeEnd;

    const size_t batch_size = particle->model->otfParams.batch_size / particle->model->otfParams.nbProcess;
    size_t current_batch = batch_size;
    do {
        setSimState(getSimStatus());
        size_t desorptions = localDesLimit;//(localDesLimit > 0 && localDesLimit > particle->tmpState.globalHits.globalHits.hit.nbDesorbed) ? localDesLimit - particle->tmpState.globalHits.globalHits.hit.nbDesorbed : 0;
        //printf("Pre[%zu] %lu + %lu / %lu\n",threadNum, desorptions, particle->tmpState.globalHits.globalHits.hit.nbDesorbed, localDesLimit);

        // 1s step
        size_t nbStep = (stepsPerSec <= 0.0) ? 250.0 : std::ceil(stepsPerSec + 0.5);
        int res = runSimulation(nbStep, &current_batch, desorptions);      // Run during 1 sec
        simEos = (res != CONT) ? true : false;

        //printf("Pos[%zu][%d] %lu + %lu / %lu\n",threadNum, simEos, desorptions, particle->tmpState.globalHits.globalHits.hit.nbDesorbed, localDesLimit);

        timeEnd = omp_get_wtime();

        bool forceQueue = timeEnd-timeLoopStart > 60 || threadNum == 0; // update after 60s of no update or when thread 0 is called
        if (procInfo->activeProcs.front() == threadNum || forceQueue) {
            size_t readdOnFail = 0;
            if(simulation->model->otfParams.desorptionLimit > 0){
                if(localDesLimit > particle->tmpState.globalHits.globalHits.nbDesorbed) {
                    localDesLimit -= particle->tmpState.globalHits.globalHits.nbDesorbed;
                    readdOnFail = particle->tmpState.globalHits.globalHits.nbDesorbed;
                }
                else localDesLimit = 0;
            }

            size_t timeOut = lastUpdateOk ? 0 : 100; //ms
            lastUpdateOk = particle->UpdateHits(simulation->globState, simulation->globParticleLog,
                                                timeOut); // Update hit with 20ms timeout. If fails, probably an other subprocess is updating, so we'll keep calculating and try it later (latest when the simulation is stopped).

            if(!lastUpdateOk) // if update failed, the desorption limit is invalid and has to be reverted
                localDesLimit += readdOnFail;

            if(procInfo->activeProcs.front() == threadNum)
                procInfo->NextSubProc();
            timeLoopStart = omp_get_wtime();
        }
        else{
            lastUpdateOk = false;
        }

        if(current_batch == 0) {
            current_batch = batch_size;

            if (!lastUpdateOk) {
                //printf("[%zu] Updating on finish!\n",threadNum);
                setSimState("Final update...");
                particle->UpdateHits(simulation->globState, simulation->globParticleLog,
                                     20000); // Update hit with 20ms timeout. If fails, probably an other subprocess is updating, so we'll keep calculating and try it later (latest when the simulation is stopped).)
            }

            if(res == HALT && desorptions <= 0)
                simEos = false;

#pragma omp atomic
            procInfo->counting_barrier++;

            /*while(procInfo->counting_barrier != procInfo->subProcInfo.size()) {
//#pragma omp flush(procInfo->counting_barrier)
                Sleep(1000); // sleep for 1 millisecond
            }*/
                // Synchronise, save convergence results and restart

            if(procInfo->activeProcs.front() == threadNum) {
                    simulation->FetchConvValues();
                    if (simulation->model->otfParams.formula_ptr &&
                        simulation->model->otfParams.formula_ptr->CheckConvergence()) {
                        simulation->isConverged = true;
                    }
            }
#pragma omp atomic
            procInfo->counting_barrier--;

        }
        //printf("[%zu] PUP: %lu , %lu , %lu\n",threadNum, desorptions,localDesLimit, particle->tmpState.globalHits.globalHits.hit.nbDesorbed);
        eos = simEos || simulation->IsConverged() ||  (this->particle->model->otfParams.timeLimit != 0 ? timeEnd-timeStart >= this->particle->model->otfParams.timeLimit : false) || (procInfo->masterCmd != COMMAND_START) || (procInfo->subProcInfo[threadNum].slaveState == PROCESS_ERROR);
    } while (!eos);

    procInfo->RemoveAsActive(threadNum);
#pragma omp atomic
    procInfo->counting_barrier++;

    if (!lastUpdateOk) {
        //printf("[%zu] Updating on finish!\n",threadNum);
        setSimState("Final update...");
        particle->UpdateHits(simulation->globState, simulation->globParticleLog,
                             20000); // Update hit with 20ms timeout. If fails, probably an other subprocess is updating, so we'll keep calculating and try it later (latest when the simulation is stopped).)
    }
    return simEos;
}

void SimThread::setSimState(char *msg) const {
    procInfo->subProcInfo[threadNum].slaveState = PROCESS_RUN;
    snprintf(procInfo->subProcInfo[threadNum].statusString, 128, "%s", msg);
}

void SimThread::setSimState(const std::string& msg) const {
    procInfo->subProcInfo[threadNum].slaveState = PROCESS_RUN;
    snprintf(procInfo->subProcInfo[threadNum].statusString, 128, "%s", msg.c_str());
}

[[nodiscard]] char *SimThread::getSimStatus() const {
    static char ret[128];
    size_t count = particle->totalDesorbed + particle->tmpState.globalHits.globalHits.nbDesorbed;;

    size_t max = 0;
    if (simulation->model->otfParams.nbProcess)
        max = (simulation->model->otfParams.desorptionLimit / simulation->model->otfParams.nbProcess)
                + ((this->threadNum < simulation->model->otfParams.desorptionLimit % simulation->model->otfParams.nbProcess) ? 1 : 0);

    if (max != 0) {
        double percent = (double) (count) * 100.0 / (double) (max);
        sprintf(ret, "MC %zd/%zd (%.1f%%)", count, max, percent);
    } else {
        sprintf(ret, "MC %zd", count);
    }
    return ret;
}

/**
* \brief A "single (1sec)" MC step of a simulation run for a given thread
 * \return 0 when simulation continues, 1 when desorption limit is reached, 2 on sync
 */
int SimThread::runSimulation(size_t nb_steps_max, size_t *nb_des_max, size_t desorptions) {

    enum EoS : int {
        CONT,
        STOP,
        HALT
    };

    {
        char msg[128];
        snprintf(msg, 128, "%s [%zu event/s]", getSimStatus(), nb_steps_max);
        setSimState(msg);
    }

    // Check end of simulation
    int eos = CONT;
    size_t remainingDes = 0;

    if (particle->model->otfParams.desorptionLimit > 0) {
        if (desorptions <= particle->tmpState.globalHits.globalHits.nbDesorbed){
            //lastHitFacet = nullptr; // reset full particle status or go on from where we left
            eos = STOP;
        }
        else {
            //if(particle->tmpState.globalHits.globalHits.hit.nbDesorbed <= (particle->model->otfParams.desorptionLimit - simulation->globState->globalHits.globalHits.hit.nbDesorbed)/ particle->model->otfParams.nbProcess)
                remainingDes = desorptions - particle->tmpState.globalHits.globalHits.nbDesorbed;
        }
    }
    //auto start_time = std::chrono::high_resolution_clock::now();
    if(eos == CONT) {

        auto des_cache = particle->tmpState.globalHits.globalHits.nbDesorbed;
        double start_time = omp_get_wtime();
        //fmt::print("Sim start with Steps {} and Des: {} : {}\n", nb_steps_max, *nb_des_max, *nb_des_max > remainingDes ? *nb_des_max : remainingDes);

        bool res = particle->SimulationMCStep(nb_steps_max, threadNum, *nb_des_max > remainingDes ? *nb_des_max : remainingDes);
        eos = (res) ? CONT : (*nb_des_max > remainingDes ? HALT : STOP);
        double end_time = omp_get_wtime();

        if(*nb_des_max) {
            //fmt::print("Des: {} : {}\n", *nb_des_max, particle->tmpState.globalHits.globalHits.nbDesorbed-des_cache);
            *nb_des_max -= std::min(*nb_des_max, particle->tmpState.globalHits.globalHits.nbDesorbed - des_cache);
        }
    //auto end_time = std::chrono::high_resolution_clock::now();


        if (res) // don't update on end, this will give a false ratio (SimMCStep could return actual steps instead of plain "false"
        {
            const double elapsedTimeMs = (end_time - start_time); //std::chrono::duration<float, std::ratio<1, 1>>(end_time - start_time).count();
            if (elapsedTimeMs != 0.0)
                stepsPerSec = (1.0 * nb_steps_max) / (elapsedTimeMs); // every 1.0 second
            else
                stepsPerSec = (100.0 * nb_steps_max); // in case of fast initial run
        }
    }
#if defined(_DEBUG)
    //printf("Running: stepPerSec = %lf [%lu]\n", stepsPerSec, threadNum);
#endif

    return eos;
}

SimulationController::SimulationController(size_t parentPID, size_t procIdx, size_t nbThreads,
                                           SimulationUnit *simulationInstance, ProcComm *pInfo) {
    this->prIdx = procIdx;
    this->parentPID = parentPID;
    if (nbThreads == 0)
#if defined(DEBUG)
        this->nbThreads = 1;
#else
        this->nbThreads = omp_get_max_threads();
#endif
    else
        this->nbThreads = nbThreads;

    simulation = simulationInstance; // TODO: Find a nicer way to manager derived simulationunit for Molflow and Synrad
    procInfo = pInfo; // Set from outside to link with general structure

    loadOk = false;

    // resetControls()
    stepsPerSec = 0.0;
    endState = false;
    lastHitUpdateOK = true;
    benchmarkADS = false;

    SetRuntimeInfo();

    //DEBUG_PRINT("Controller created (%zd bytes)\n", sizeof(SHCONTROL));
    SetReady(false);
}

SimulationController::~SimulationController() = default;

SimulationController::SimulationController(SimulationController &&o) noexcept {
    simulation = o.simulation;
    procInfo = o.procInfo;

    o.simulation = nullptr;
    o.procInfo = nullptr;

    stepsPerSec = o.stepsPerSec;
    endState = o.endState;
    lastHitUpdateOK = o.lastHitUpdateOK;

    prIdx = o.prIdx;
    parentPID = o.parentPID;
    nbThreads = o.nbThreads;
    loadOk = o.loadOk;
    benchmarkADS = o.benchmarkADS;

    simThreads.reserve(nbThreads);
    for (size_t t = 0; t < nbThreads; t++) {
        simThreads.emplace_back(
                SimThread(procInfo, simulation, t));
        simThreads.back().particle = simulation->GetParticle(t);
    }
}

int SimulationController::resetControls() {
    lastHitUpdateOK = true;
    endState = false;

    stepsPerSec = 0.0;

    return 0;
}

int SimulationController::SetRuntimeInfo() {

    // Update runtime information
    for (auto &pInfo : procInfo->subProcInfo)
        GetProcInfo(pInfo.procId, &pInfo.runtimeInfo);

    return 0;
}

int SimulationController::ClearCommand() {

    procInfo->masterCmd = COMMAND_NONE;
    procInfo->cmdParam = 0;
    procInfo->cmdParam2 = 0;
    SetState(PROCESS_READY,GetSimuStatus());

    return 0;
}

int SimulationController::SetState(size_t state, const char *status, bool changeState, bool changeStatus) {

    if (changeState) {
        DEBUG_PRINT("setstate %zd \n", state);
        //master->procInformation[prIdx].masterCmd = state;
        for (auto &pInfo : procInfo->subProcInfo) {
            pInfo.slaveState = state;
        }
    }
    if (changeStatus) {
        for (auto &pInfo : procInfo->subProcInfo) {
            strncpy(pInfo.statusString, status, 127);
            pInfo.statusString[127] = '\0';
        }
    }

    if(state == PROCESS_ERROR){
        procInfo->masterCmd = PROCESS_WAIT;
    }
    return 0;
}

int SimulationController::SetState(size_t state, const std::vector<std::string> &status, bool changeState, bool changeStatus) {

    if (changeState) {
        DEBUG_PRINT("setstate %zd \n", state);
        //master->procInformation[prIdx].masterCmd = state;
        for (auto &pInfo : procInfo->subProcInfo) {
            pInfo.slaveState = state;
        }
    }
    if (changeStatus) {
        if(procInfo->subProcInfo.size() != status.size()){
            for (auto &pInfo : procInfo->subProcInfo) {
                strncpy(pInfo.statusString, "invalid state", 127);
                pInfo.statusString[127] = '\0';
            }
        }
        else {
            size_t pInd = 0;
            for (auto &pInfo : procInfo->subProcInfo) {
                strncpy(pInfo.statusString, status[pInd++].c_str(), 127);
                pInfo.statusString[127] = '\0';
            }
        }
    }

    if(state == PROCESS_ERROR){
        procInfo->masterCmd = PROCESS_WAIT;
    }
    return 0;
}

std::vector<std::string> SimulationController::GetSimuStatus() {

    std::vector<std::string> ret_vec(nbThreads);
    if (!simulation) {
        ret_vec.assign(nbThreads, "[NONE]");
    }
    else{
        size_t threadId = 0;
        auto* sim = simulation;
        for(size_t p = 0; p < nbThreads; ++p) {
            auto* particle = (simulation)->GetParticle(p);
            if(particle == nullptr) break;
            if(!particle->tmpState.initialized){
                ret_vec[threadId]= "[NONE]";
            }
            else {
                size_t count = 0;
                count = particle->totalDesorbed + particle->tmpState.globalHits.globalHits.nbDesorbed;
                size_t max = 0;
                if (sim->model->otfParams.nbProcess)
                    max = sim->model->otfParams.desorptionLimit / sim->model->otfParams.nbProcess + ((threadId < sim->model->otfParams.desorptionLimit % sim->model->otfParams.nbProcess) ? 1 : 0);

                char tmp[128];
                if (max != 0) {
                    double percent = (double) (count) * 100.0 / (double) (max);
                    sprintf(tmp, "MC %zd/%zd (%.1f%%)", count, max, percent);
                } else {
                    sprintf(tmp, "MC %zd", count);
                }
                ret_vec[threadId] = tmp;
            }
            ++threadId;
        }
    }
    return ret_vec;
}

void SimulationController::SetErrorSub(const char *message) {
    Log::console_error("Error: {}\n", message);
    SetState(PROCESS_ERROR, message);
}

void SimulationController::SetStatus(char *status) {
    for (auto &pInfo : procInfo->subProcInfo) {
        strncpy(pInfo.statusString, status, 127);
        pInfo.statusString[127] = '\0';
    }
}

void SimulationController::SetReady(const bool loadOk) {

    if (loadOk)
        SetState(PROCESS_READY, this->GetSimuStatus());
    else
        SetState(PROCESS_READY, "(No geometry)");

    ClearCommand();
}

// Return Error code if subproc info is out of sync
size_t SimulationController::GetLocalState() const {
    size_t locState = procInfo->subProcInfo[0].slaveState;
    for (auto &pInfo : procInfo->subProcInfo) {
        if (locState != pInfo.slaveState)
            return PROCESS_ERROR;
    }
    return locState;
}

// Main loop
int SimulationController::controlledLoop(int argc, char **argv) {
    endState = false;
    loadOk = false;
    while (!endState) {
        bool eos = false;
        switch (procInfo->masterCmd) {

            case COMMAND_LOAD: {
                Load();
                break;
            }
            case COMMAND_UPDATEPARAMS: {
                SetState(PROCESS_WAIT, GetSimuStatus());
                DEBUG_PRINT("[%d] COMMAND: UPDATEPARAMS (%zd,%zd)\n", prIdx, procInfo->cmdParam, procInfo->cmdParam2);
                if (UpdateParams()) {
                    SetReady(loadOk);
                    //SetState(procInfo->cmdParam2, GetSimuStatus());
                } else {
                    SetState(PROCESS_ERROR, "Could not update parameters");
                }
                break;
            }
            case COMMAND_START: {
                Start();
                break;
            }
            case COMMAND_PAUSE: {
                DEBUG_PRINT("[%d] COMMAND: PAUSE (%zd,%zu)\n", prIdx, procInfo->cmdParam, procInfo->cmdParam2);
                if (!lastHitUpdateOK) {
                    // Last update not successful, retry with a longer timeout
                    if ((GetLocalState() != PROCESS_ERROR)) {
                        SetState(PROCESS_STARTING, "Updating hits...", true, true);
                        //lastHitUpdateOK = simulation->UpdateHits(prIdx, 60000);
                        SetState(PROCESS_STARTING, GetSimuStatus(), true, true);
                    }
                }
                SetReady(loadOk);
                break;
            }
            case COMMAND_RESET: {
                Reset();
                break;
            }
            case COMMAND_EXIT: {
                DEBUG_PRINT("[%d] COMMAND: EXIT (%zd,%zu)\n", prIdx, procInfo->cmdParam, procInfo->cmdParam2);
                endState = true;
                break;
            }
            case COMMAND_CLOSE: {
                DEBUG_PRINT("[%d] COMMAND: CLOSE (%zd,%zu)\n", prIdx, procInfo->cmdParam, procInfo->cmdParam2);
                auto* sim = simulation;
                //for (auto &sim : *simulation)
                    sim->ClearSimulation();
                loadOk = false;
                SetReady(loadOk);
                break;
            }
            default: {
                ProcessSleep(WAITTIME);
                break;
            }
        }
    }
    SetState(PROCESS_KILLED, "Process terminated peacefully");
    return 0;
}

/**
* \brief Call a rebuild of the ADS
 * \return 0> error code, 0 when ok
 */
int SimulationController::RebuildAccel() {
    if (simulation->RebuildAccelStructure()) {
        return 1;
    }
    return 0;
}

/**
* \brief Load and init simulation geometry and initialize threads after previous sanity check
 * Setup inidividual particles per thread and local desorption limits
 * \return true on error, false when ok
 */
bool SimulationController::Load() {
    DEBUG_PRINT("[%d] COMMAND: LOAD (%zd,%zu)\n", prIdx, procInfo->cmdParam, procInfo->cmdParam2);
    SetState(PROCESS_STARTING, "Loading simulation");

    auto sane = simulation->SanityCheckModel(false);
    if(!sane.first) {
        SetState(PROCESS_STARTING, "Loading simulation");
        bool loadError = false;

        // Init particles / threads
        simulation->SetNParticle(nbThreads, false);
        if (simulation->LoadSimulation(procInfo->subProcInfo[0].statusString)) {
            loadError = true;
        }

        if (!loadError) { // loadOk = Load();
            loadOk = true;

            {//if(nbThreads != simThreads.size()) {
                simThreads.clear();
                simThreads.reserve(nbThreads);
                for (size_t t = 0; t < nbThreads; t++) {
                    simThreads.emplace_back(
                            SimThread(procInfo, simulation, t));
                    simThreads.back().particle = simulation->GetParticle(t);
                }
            }

            InitOMP();

            // Calculate remaining work
            size_t desPerThread = 0;
            size_t remainder = 0;
            size_t des_global = simulation->globState->globalHits.globalHits.nbDesorbed;
            if (des_global > 0) {
                desPerThread = des_global / nbThreads;
                remainder = des_global % nbThreads;
            }
            for (auto &thread : simThreads) {
                thread.particle->totalDesorbed = desPerThread;
                thread.particle->totalDesorbed += (thread.threadNum < remainder) ? 1 : 0;
            }

            SetRuntimeInfo();
        }
    }
    else {
        loadOk = false;
        SetState(PROCESS_ERROR, sane.second->c_str());
    }
    SetReady(loadOk);

    return !loadOk;
}

/**
* \brief Update on the fly parameters when called from the GUI
 * \return true on success
 */
bool SimulationController::UpdateParams() {
    // Load geometry
    auto* sim = simulation;
    if (sim->model->otfParams.enableLogging) {
        Log::console_msg(3, "Logging with size limit {}\n",
               sizeof(size_t) + sim->model->otfParams.logLimit * sizeof(ParticleLoggerItem));
    }
    sim->ReinitializeParticleLog();
    return true;
}

/**
* \brief Start the simulation after previous sanity checks
 * \return 0> error code, 0 when ok
 */
int SimulationController::Start() {

    // Check simulation model and geometry one last time
    auto sane = simulation->SanityCheckModel(true);
    if(sane.first){
        loadOk = false;
    }

    for(auto& thread : simThreads){
        if(!thread.particle)
            loadOk = false;
    }

    if(!loadOk) {
        if(sane.second)
            SetState(PROCESS_ERROR, sane.second->c_str());
        else
            SetState(PROCESS_ERROR, GetSimuStatus());
        return 1;
    }

    if(simulation->model->otfParams.raySampling) {
        simulation->globState->UpdateBatteryFrequencies();
        for (auto &freq: simThreads) {
            freq.particle->tmpState.hitBattery.rays = simulation->globState->hitBattery.rays;
            freq.particle->tmpState.hitBattery.grays = simulation->globState->hitBattery.grays;
        }
    }
    else{
        for (auto &thr: simThreads) {
            thr.particle->tmpState.StopBatteryChange();
        }
    }

#if defined(USE_OLD_BVH)
    if(simulation->model->structures.empty()/* && RebuildAccel()*/){
        loadOk = false;
        SetState(PROCESS_ERROR, "Failed building (old) acceleration structure!");
        return 1;
    }
#endif
    if(simulation->model->accel.empty()/* && RebuildAccel()*/){
        loadOk = false;
        SetState(PROCESS_ERROR, "Failed building acceleration structure!");
        return 1;
    }

    bool setPayload = false;
    for(auto& accel : simulation->model->accel){
        if(dynamic_cast<KdTreeAccel*>(accel.get())){
            if(dynamic_cast<KdTreeAccel*>(accel.get())->restartFromNode)
                setPayload = true;
        }
    }

    for (auto &thr: simThreads) {
        if(setPayload) {
            if (!thr.particle->particle.pay)
                thr.particle->particle.pay = new RopePayload;
        }
        else
            thr.particle->particle.pay = nullptr;
        /*else if (thr.particle->particle.pay)
            delete thr.particle->particle.pay;*/
    }

    if (simulation->model->otfParams.desorptionLimit > 0) {
        if (simulation->totalDesorbed >=
            simulation->model->otfParams.desorptionLimit /
            simulation->model->otfParams.nbProcess) {
            ClearCommand();
            SetState(PROCESS_DONE, GetSimuStatus());
        }
    }

    if (GetLocalState() != PROCESS_RUN) {
        DEBUG_PRINT("[%d] COMMAND: START (%zd,%zu)\n", prIdx, procInfo->cmdParam, procInfo->cmdParam2);
        SetState(PROCESS_RUN, GetSimuStatus());
    }


    bool eos = false;
    bool lastUpdateOk = true;
    if (loadOk) {
        procInfo->InitActiveProcList();

        // Calculate remaining work
        size_t desPerThread = 0;
        size_t remainder = 0;
        if(simulation->model->otfParams.desorptionLimit > 0){
            if(simulation->model->otfParams.desorptionLimit > (simulation->globState->globalHits.globalHits.nbDesorbed)) {
                size_t limitDes_global = simulation->model->otfParams.desorptionLimit;
                desPerThread = limitDes_global / nbThreads;
                remainder = limitDes_global % nbThreads;
            }
        }
        for(auto& thread : simThreads){
            size_t localDes = (desPerThread > thread.particle->totalDesorbed) ? desPerThread - thread.particle->totalDesorbed : 0;
            thread.localDesLimit = (thread.threadNum < remainder) ? localDes + 1 : localDes;
        }

        // Start with a benchmark
        if(simulation->model->otfParams.benchmarkADS)
            simulation->FindBestADS();

        int simuEnd = 0; // bool atomic is not supported by MSVC OMP implementation
#pragma omp parallel num_threads(nbThreads) default(none) firstprivate(/*stepsPerSec,*/ lastUpdateOk, eos) shared(/*procInfo,*/  simuEnd/*, simulation,simThreads*/)
        {
            // Set OpenMP thread priority on Windows whenever we start a simulation run
#if defined(_WIN32) && defined(_MSC_VER)
            SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_IDLE);
#endif
            eos = simThreads[omp_get_thread_num()].runLoop();

            if(eos) {
#pragma omp atomic
                simuEnd += 1;
            }
        }

        if (simuEnd) {
            if (GetLocalState() != PROCESS_ERROR && GetLocalState() != PROCESS_KILLED) {
                // Max desorption reached
                ClearCommand();
                SetState(PROCESS_DONE, GetSimuStatus());
                DEBUG_PRINT("[%d] COMMAND: PROCESS_DONE (Max reached)\n", prIdx);
            }
        }
        else {
            if (GetLocalState() != PROCESS_ERROR && GetLocalState() != PROCESS_KILLED) {
                // Time limit reached
                ClearCommand();
                SetState(PROCESS_DONE, GetSimuStatus());
                DEBUG_PRINT("[%d] COMMAND: PROCESS_DONE (Stopped)\n", prIdx);
            }
        }
    } else {
        SetErrorSub("No geometry loaded");
        ClearCommand();
    }
    return 0;
}

int SimulationController::Reset() {
    DEBUG_PRINT("[%d] COMMAND: RESET (%zd,%zu)\n", prIdx, procInfo->cmdParam, procInfo->cmdParam2);
    SetState(PROCESS_STARTING, "Resetting local cache...", true, true);
    resetControls();

    auto *sim = simulation;
    sim->ResetSimulation();
    SetReady(loadOk);

    return 0;
}

void SimulationController::EmergencyExit(){
    for(auto& t : simThreads)
        t.particle->allQuit = true;
}