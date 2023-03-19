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

#ifndef MOLFLOW_PROJ_PROCESSCONTROL_H
#define MOLFLOW_PROJ_PROCESSCONTROL_H

#include <cstddef> //size_t
#include <vector>
#include <mutex>
#include <list>

#define PROCESS_STARTING 0   // Loading state
#define PROCESS_RUN      1   // Running state
#define PROCESS_READY    2   // Waiting state
#define PROCESS_KILLED   3   // Process killed
#define PROCESS_ERROR    4   // Process in error
#define PROCESS_DONE     5   // Simulation ended
#define PROCESS_WAIT     6   // Command fully executed

#define COMMAND_NONE     10  // No change
#define COMMAND_LOAD     11  // Load geometry
#define COMMAND_START    12  // Start simu
#define COMMAND_PAUSE    13  // Pause simu
#define COMMAND_RESET    14  // Reset simu
#define COMMAND_EXIT     15  // Exit
#define COMMAND_CLOSE    16  // Release handles
#define COMMAND_UPDATEPARAMS 17 //Update simulation mode (low flux, fluxwise/powerwise, displayed regions)
#define COMMAND_ACCEL     18  // Build Accel

static const char *prStates[] = {

        "Not started",
        "Running",
        "Waiting",
        "Killed",
        "Error",
        "Done",
        "",
        "",
        "",
        "",
        "No command",
        "Loading",
        "Starting",
        "Stopping",
        "Resetting",
        "Exiting",
        "Closing",
        "Update params"
};

struct PROCESS_INFO{

    double cpu_time; // CPU time         (in second)
    size_t  mem_use;  // Memory usage     (in byte)
    size_t  mem_peak; // MAx Memory usage (in byte)

};

//! Process Communication class for handling inter process/thread communication
struct ProcComm {

    struct SubProcInfo {
        size_t procId;
        size_t slaveState;
        char statusString[128];
        PROCESS_INFO runtimeInfo;
    };
    int counting_barrier{0};
    size_t masterCmd;
    size_t cmdParam;
    size_t cmdParam2;

    std::list<size_t> activeProcs;
    std::mutex m;
    std::vector<SubProcInfo> subProcInfo;

    ProcComm();
    explicit ProcComm(size_t nbProcs) : ProcComm() {
        Resize(nbProcs);
    };
    void Resize(size_t nbProcs){
        subProcInfo.resize(nbProcs);
        InitActiveProcList();
    };

    void NextSubProc();

    ProcComm& operator=(const ProcComm & src);
    ProcComm& operator=(ProcComm && src) noexcept ;

    void RemoveAsActive(size_t id);

    void InitActiveProcList();
};

#endif //MOLFLOW_PROJ_PROCESSCONTROL_H
