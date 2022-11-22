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

#ifndef MOLFLOW_PROJ_SIMULATIONUNIT_H
#define MOLFLOW_PROJ_SIMULATIONUNIT_H

#include "SMP.h"
#include "Buffer_shared.h"

class SimulationModel;
class GlobalSimuState;
struct ParticleLog;
struct ConvergenceData;

namespace MFSim {
    class Particle;
}

/**
* \brief Abstract Simulation unit that is implemented for the CPU based simulations for Synrad and Molflow
 */
class SimulationUnit {
public:
    SimulationUnit() : model(), totalDesorbed(0), m(){
        globState = nullptr;
        globParticleLog = nullptr;
    };
    SimulationUnit(const SimulationUnit& o) : model(o.model) , m() {
        globState = nullptr;
        globParticleLog = nullptr;
        totalDesorbed = o.totalDesorbed;
    };
    virtual ~SimulationUnit()= default;

    /*! Parse input and pre compute/prepare all necessary structures  */
    virtual size_t LoadSimulation(char *loadStatus) = 0;
    virtual int RebuildAccelStructure() = 0;

    virtual int ReinitializeParticleLog() = 0;
    virtual int ReinitializeConvergenceLog() = 0;
    virtual int FetchConvValues() = 0;
    virtual int IsConverged() = 0;
    virtual std::pair<int, std::optional<std::string>> SanityCheckModel(bool strictCheck) = 0;

    virtual void ResetSimulation() = 0;
    virtual void ClearSimulation() = 0;

    virtual size_t GetHitsSize() = 0;
    virtual MFSim::Particle * GetParticle(size_t i) = 0;
    virtual void SetNParticle(size_t n, bool fixedSeed) = 0;

    virtual void FindBestADS() = 0; //! benchmark to find best Acceleration Data Structure

public:
    std::shared_ptr<SimulationModel> model;
    //OntheflySimulationParams ontheflyParams;
    //GeomProperties sh;
    // Particle coordinates (MC)
    GlobalSimuState* globState;
    ParticleLog* globParticleLog; //Recorded particle log since last UpdateMCHits
    std::vector<ConvergenceData> convergenceValues; // One vector of nbDesorption,formulaValue pairs for each formula

    //GlobalSimuState tmpResults;

    bool isConverged{false};
    size_t totalDesorbed; // todo: should be a "sim counter"
    std::timed_mutex m;

};

#endif //MOLFLOW_PROJ_SIMULATIONUNIT_H
