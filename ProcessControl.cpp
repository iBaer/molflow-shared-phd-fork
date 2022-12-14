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

#include "ProcessControl.h"
#include <cstring>

ProcComm::ProcComm() : m() {
    masterCmd = 0;
    cmdParam = 0;
    cmdParam2 = 0;
}
/**
* \brief Assign operator
* \param src reference to source object
* \return address of this
*/
ProcComm& ProcComm::operator=(const ProcComm & src) {
    masterCmd = src.masterCmd;
    cmdParam = src.cmdParam;
    cmdParam2 = src.cmdParam2;
    subProcInfo = src.subProcInfo;
    return *this;
}

/**
* \brief Assign operator
* \param src reference to source object
* \return address of this
*/
ProcComm& ProcComm::operator=(ProcComm && src) noexcept {
    masterCmd = src.masterCmd;
    cmdParam = src.cmdParam;
    cmdParam2 = src.cmdParam2;
    subProcInfo = std::move(src.subProcInfo);
    return *this;
}

//! Moves first sub process to the back of the "active", for round robin fashion of communication for updates
void ProcComm::NextSubProc() {
    this->m.lock();
    activeProcs.emplace_back(activeProcs.front());
    activeProcs.pop_front();
    this->m.unlock();
}

//! Removes a process from the active list, in case it is finished
void ProcComm::RemoveAsActive(size_t id) {
    this->m.lock();
    for(auto proc = activeProcs.begin(); proc != activeProcs.end(); ++proc){
        if(id == (*proc)) {
            activeProcs.erase(proc);
            break;
        }
    }
    this->m.unlock();
}

//! Init list of active/simulating processes
void ProcComm::InitActiveProcList() {
    this->m.lock();
    activeProcs.clear();
    for(size_t id = 0; id < this->subProcInfo.size(); ++id)
        activeProcs.emplace_back(id);
    this->m.unlock();
}
