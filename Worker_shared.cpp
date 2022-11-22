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
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
#define NOMINMAX
//#include <Windows.h>
//#include <Process.h>
//#include <direct.h>
#else

#endif

#include <cmath>
#include <cstdlib>
#include <algorithm>

#include "Worker.h"
#include "Facet_shared.h"
#include "GLApp/GLApp.h"
#include "GLApp/GLMessageBox.h"
#include "Helper/MathTools.h" //Min max
#include "GLApp/GLList.h"

#include "GLApp/GLUnitDialog.h"
#include "Interface/LoadStatus.h"
//#include "ProcessControl.h" // defines for process commands
//#include "SimulationManager.h"
//#include "Buffer_shared.h"

#include "SimulationFacet.h"

#if defined(MOLFLOW)
#include "../src/MolFlow.h"
#include "../src/MolflowGeometry.h"
#include "../src/Interface/FacetAdvParams.h"
#endif

#if defined(SYNRAD)
#include "../src/SynRad.h"
#include "../src/SynradGeometry.h"
#endif

//#include "File.h" //File utils (Get extension, etc)

/*
//Leak detection
#if defined(_DEBUG)
#define _CRTDBG_MAP_ALLOC
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif
*/
using namespace pugi;

#if defined(MOLFLOW)
extern MolFlow *mApp;
#endif

#if defined(SYNRAD)
extern SynRad*mApp;
#endif

Worker::~Worker() {
    simManager.KillAllSimUnits(); // kill in case of preemptive Molflow close, prevents updates on already deleted globalsimustate
    delete geom;
}

Geometry *Worker::GetGeometry() {
    return geom;
}

std::string Worker::GetCurrentFileName() const {
    return fullFileName;
}

// Given an absolute path (e.g. dir1\\dir2\\filename) it extracts and returns only the filename
std::string Worker::GetCurrentShortFileName() const {
    std::string::size_type filePos = fullFileName.rfind('/');
    if (filePos != std::string::npos)
        ++filePos;
    else {
        filePos = fullFileName.rfind('\\');
        if (filePos != std::string::npos)
            ++filePos;
        else
            filePos = 0;
    }

    return fullFileName.substr(filePos);
}

/*
char *Worker::GetShortFileName(char* longFileName) {

	static char ret[512];
	char *r = strrchr(longFileName, '/');
	if (!r) r = strrchr(longFileName, '\\');
	if (!r) strcpy(ret, longFileName);
	else   {
		r++;
		strcpy(ret, r);
	}

	return ret;

}
*/

void Worker::SetCurrentFileName(const char *fileName) {
    fullFileName = fileName;
}

void Worker::ExportTextures(const char *fileName, int grouping, int mode, bool askConfirm, bool saveSelected) {

    // Read a file
    FILE *f = fopen(fileName, "w");
    if (!f) {
        char tmp[256];
        sprintf(tmp, "Cannot open file for writing %s", fileName);
        throw Error(tmp);
    }

    //bool buffer_old = simManager.GetLockedHitBuffer();
#if defined(MOLFLOW)
    geom->ExportTextures(f, grouping, mode, globState, saveSelected);
#endif
#if defined(SYNRAD)
    geom->ExportTextures(f, grouping, mode, no_scans, globState, saveSelected);
#endif
    //simManager.UnlockHitBuffer();
    fclose(f);
}

void Worker::Stop_Public() {
    // Stop
    InnerStop(mApp->m_fTime);
    try {
        Stop();
        Update(mApp->m_fTime);
    }
    catch (const std::exception &e) {
        GLMessageBox::Display(e.what(), "Error (Stop)", GLDLG_OK, GLDLG_ICONERROR);
    }
}

void Worker::ReleaseHits() {
    simManager.UnlockHitBuffer();
}

bool Worker::GetHits() {
    try {
        if (needsReload) RealReload();
    }
    catch (const std::exception &e) {
        GLMessageBox::Display(e.what(), "Error (Stop)", GLDLG_OK, GLDLG_ICONERROR);
        return false;
    }
    return true;
}

ParticleLog & Worker::GetLog() {
    try {
        if(!particleLog.tMutex.try_lock_for(std::chrono::seconds(1)))
            throw std::runtime_error("Couldn't get log access");
        if (needsReload) RealReload();
    }
    catch (const std::exception &e) {
        GLMessageBox::Display(e.what(), "Error (GetLog)", GLDLG_OK, GLDLG_ICONERROR);
    }
    return particleLog;
}

void Worker::UnlockLog() {
    particleLog.tMutex.unlock();
}

void Worker::ThrowSubProcError(const char *message) {

    char errMsg[1024];
    if (!message)
        sprintf(errMsg, "Bad response from sub process(es):\n%s",GetErrorDetails().c_str());
    else
        sprintf(errMsg, "%s\n%s", message, GetErrorDetails().c_str());
    throw std::runtime_error(errMsg);

}

void Worker::Reload() {
    //Schedules a reload
    //Actual reloading is done by RealReload() method
    needsReload = true;
}

/*
void Worker::SetMaxDesorption(size_t max) {

	try {
		ResetStatsAndHits(0.0);

		desorptionLimit = max;
		Reload();

	}
	catch (const std::exception &e) {
		GLMessageBox::Display(e.what(), "Error", GLDLG_OK, GLDLG_ICONERROR);
	}

}
*/

std::string Worker::GetErrorDetails() {
    return simManager.GetErrorDetails();
}

void Worker::InnerStop(float appTime) {
    simuTimer.Stop();
}

void Worker::StartStop(float appTime) {
    if( IsRunning() )  {

        // Stop
        InnerStop(appTime);
        try {
            Stop();
            Update(appTime);

        }
        catch(const std::exception &e) {
            GLMessageBox::Display((char *)e.what(),"Error (Stop)",GLDLG_OK,GLDLG_ICONERROR);

        }
    }
    else {

        // Start
        try {
            if (needsReload) RealReload(); //Synchronize subprocesses to main process
            Start();
            simuTimer.Start();
        }
        catch (const std::exception &e) {
            //isRunning = false;
            GLMessageBox::Display(e.what(),"Error (Start)",GLDLG_OK,GLDLG_ICONERROR);
            return;
        }

        // Particular case when simulation ends before getting RUN state
        if(simManager.allProcsDone) {
            Update(appTime);
            GLMessageBox::Display("Max desorption reached","Information (Start)",GLDLG_OK,GLDLG_ICONINFO);
        }
    }

}

void Worker::ResetStatsAndHits(float appTime) {

    simuTimer.ReInit();
    //stopTime = 0.0f;
    //startTime = 0.0f;
    //simuTime = 0.0f;
    if (model->otfParams.nbProcess == 0)
        return;

    try {
        // First reset local states, then sim handles (to communicate cleared stats)
        ResetWorkerStats();
        simManager.ResetHits();

        if (needsReload) RealReload();
        Update(appTime);
    }
    catch (const std::exception &e) {
        GLMessageBox::Display(e.what(), "Error", GLDLG_OK, GLDLG_ICONERROR);
    }
}

void Worker::Stop() {
    try {
        if (simManager.StopSimulation()) {
            throw std::logic_error("No active simulation to stop!");
        }
    }
    catch (const std::exception &e) {
        throw Error(e.what());
    }
}

void Worker::InitSimProc() {

    simManager.useCPU = true;
    simManager.nbThreads = 0; // set to 0 to init max threads

    // Launch n subprocess
    if(simManager.InitSimUnits()) {
        throw Error("Initialising simulation unit failed!");
    }

    model->otfParams.nbProcess = simManager.nbThreads;

    //if (!mApp->loadStatus) mApp->loadStatus = new LoadStatus(this);
}

void Worker::SetProcNumber(size_t n) {

    // Kill all sub process
    try{
        simManager.KillAllSimUnits();
    }
    catch (const std::exception &e) {
        throw Error("Killing subprocesses failed!");
    }

    simManager.useCPU = true;
    simManager.nbThreads = std::clamp((size_t)n , (size_t)0 , MAX_PROCESS);

    // Launch n subprocess
    if ((model->otfParams.nbProcess = simManager.InitSimUnits())) {
        throw Error("Starting subprocesses failed!");
    }

    model->otfParams.nbProcess = simManager.nbThreads;

    //if (!mApp->loadStatus) mApp->loadStatus = new LoadStatus(this);
}

size_t Worker::GetPID(size_t prIdx) {
    return 0;
}

#ifdef MOLFLOW
void Worker::CalculateTextureLimits(){
    // first get tmp limit

    MolflowSimulationModel* mf_model = dynamic_cast<MolflowSimulationModel *>(model.get());
    TEXTURE_MIN_MAX limits[3];
    for(auto& lim : limits){
        lim.max.steady_state = lim.max.moments_only = 0;
        lim.min.steady_state = lim.min.moments_only = HITMAX;
    }

    for (const auto &sub : model->facets) {
        auto& subF = *sub;
        if (subF.sh.isTextured) {
            for (size_t m = 0; m < ( 1 + mf_model->tdParams.moments.size()); m++) {
                {
                    // go on if the facet was never hit before
                    auto &facetHitBuffer = globState.facetStates[subF.globalId].momentResults[m].hits;
                    if (facetHitBuffer.nbMCHit == 0 && facetHitBuffer.nbDesorbed == 0) continue;
                }

                //double dCoef = globState.globalHits.globalHits.hit.nbDesorbed * 1E4 * model->wp.gasMass / 1000 / 6E23 * MAGIC_CORRECTION_FACTOR;  //1E4 is conversion from m2 to cm2
                const double timeCorrection =
                        m == 0 ? model->wp.finalOutgassingRate : (model->wp.totalDesorbedMolecules) /
                                moments[m - 1].second;//model->tdParams.moments[m - 1].second;
                //model->wp.timeWindowSize;
                //Timecorrection is required to compare constant flow texture values with moment values (for autoscaling)
                const auto &texture = globState.facetStates[subF.globalId].momentResults[m].texture;
                const size_t textureSize = texture.size();
                for (size_t t = 0; t < textureSize; t++) {
                    //Add temporary hit counts

                    if (subF.largeEnough[t]) {
                        double val[3];  //pre-calculated autoscaling values (Pressure, imp.rate, density)

                        val[0] = texture[t].sum_v_ort_per_area * timeCorrection; //pressure without dCoef_pressure
                        val[1] = texture[t].countEquiv * subF.textureCellIncrements[t] * timeCorrection; //imp.rate without dCoef
                        val[2] = texture[t].sum_1_per_ort_velocity * subF.textureCellIncrements[t] * timeCorrection; //particle density without dCoef

                        //Global autoscale
                        for (int v = 0; v < 3; v++) {
                            if (m == 0) {
                                limits[v].max.steady_state = std::max(val[v], limits[v].max.steady_state);

                                if (val[v] > 0.0) {
                                    limits[v].min.steady_state = std::min(val[v], limits[v].min.steady_state);
                                }
                            }
                            //Autoscale ignoring constant flow (moments only)
                            else { //if (m != 0)
                                limits[v].max.moments_only = std::max(val[v],limits[v].max.moments_only);;

                                if (val[v] > 0.0)
                                    limits[v].min.moments_only = std::min(val[v],limits[v].min.moments_only);;
                            }
                        }
                    } // if largeenough
                }
            }
        }
    }

    double dCoef_custom[] = { 1.0, 1.0, 1.0 };  //Three coefficients for pressure, imp.rate, density
    //Autoscaling limits come from the subprocess corrected by "time factor", which makes constant flow and moment values comparable
    //Time correction factor in subprocess: MoleculesPerTP * nbDesorbed
    dCoef_custom[0] = 1E4 / (double)globState.globalHits.globalHits.nbDesorbed * mApp->worker.model->wp.gasMass / 1000 / 6E23*0.0100; //multiplied by timecorr*sum_v_ort_per_area: pressure
    dCoef_custom[1] = 1E4 / (double)globState.globalHits.globalHits.nbDesorbed;
    dCoef_custom[2] = 1E4 / (double)globState.globalHits.globalHits.nbDesorbed;

    // Add coefficient scaling
    for(int v = 0; v < 3; ++v) {
        limits[v].max.steady_state *= dCoef_custom[v];
        limits[v].max.moments_only *= dCoef_custom[v];
        limits[v].min.steady_state *= dCoef_custom[v];
        limits[v].min.moments_only *= dCoef_custom[v];
    }

    // Last put temp limits into global struct
    for(int v = 0; v < 3; ++v) {
        GetGeometry()->texture_limits[v].autoscale = limits[v];
    }
}
#endif

#ifdef SYNRAD
#define HITMAX 1E38
void Worker::CalculateTextureLimits(){
    // first get tmp limit
    TextureCell limitMin, limitMax;
    TEXTURE_MIN_MAX limits[3]; // count, flux, power
    for(auto& lim : limits){
        lim.max = 0;
        lim.min = HITMAX;
    }

    for (const auto &sub : model->facets) {
        auto& subF = *sub;
        if (subF.sh.isTextured) {
                {
                    // go on if the facet was never hit before
                    auto &facetHitBuffer = globState.facetStates[subF.globalId].momentResults[0].hits;
                    if (facetHitBuffer.nbMCHit == 0 && facetHitBuffer.nbDesorbed == 0) continue;
                }

                const auto &texture = globState.facetStates[subF.globalId].momentResults[0].texture;
                const size_t textureSize = texture.size();
                for (size_t t = 0; t < textureSize; t++) {
                    //Add temporary hit counts

                    if (subF.largeEnough[t]) { // TODO: For count it wasn't applied in Synrad so far
                        double val[3];  //pre-calculated autoscaling values (Pressure, imp.rate, density)

                        val[0] = texture[t].count;
                        val[1] = texture[t].flux /** subF.textureCellIncrements[t]*/;
                        val[2] = texture[t].power /** subF.textureCellIncrements[t]*/;

                        //Global autoscale
                        for (int v = 0; v < 3; v++) {
                            limits[v].max = std::max(val[v], limits[v].max);

                            if (val[v] > 0.0) {
                                limits[v].min = std::min(val[v], limits[v].min);
                            }
                        }
                    } // if largeenough
                }
        }
    }

    /*double dCoef_custom[] = { 1.0, 1.0, 1.0 };  //Three coefficients for pressure, imp.rate, density
    //Autoscaling limits come from the subprocess corrected by "time factor", which makes constant flow and moment values comparable
    //Time correction factor in subprocess: MoleculesPerTP * nbDesorbed
    dCoef_custom[0] = 1E4 / (double)globState.globalHits.globalHits.nbDesorbed * mApp->worker.model->wp.gasMass / 1000 / 6E23*0.0100; //multiplied by timecorr*sum_v_ort_per_area: pressure
    dCoef_custom[1] = 1E4 / (double)globState.globalHits.globalHits.nbDesorbed;
    dCoef_custom[2] = 1E4 / (double)globState.globalHits.globalHits.nbDesorbed;

    // Add coefficient scaling
    for(int v = 0; v < 3; ++v) {
        limits[v].max.steady_state *= dCoef_custom[v];
        limits[v].max.moments_only *= dCoef_custom[v];
        limits[v].min.steady_state *= dCoef_custom[v];
        limits[v].min.moments_only *= dCoef_custom[v];
    }*/

    // Last put temp limits into global struct
    for(int v = 0; v < 3; ++v) {
        GetGeometry()->texture_limits[v].autoscale = limits[v];
    }
}
#endif



void Worker::RebuildTextures() {

    if (mApp->needsTexture || mApp->needsDirection) {

        // Only reload if we are even rebuilding textures
        if (needsReload)
            RealReload();

        bool buffer_old = simManager.GetLockedHitBuffer();
        if (!buffer_old)
            return;


        try {
            CalculateTextureLimits();
            geom->BuildFacetTextures(globState,mApp->needsTexture,mApp->needsDirection);
        }
        catch (const std::exception &e) {
            simManager.UnlockHitBuffer();
            throw;
        }
    }
    simManager.UnlockHitBuffer();
}

size_t Worker::GetProcNumber() const {
    //return model->otfParams.nbProcess;
    return simManager.nbThreads;
}

bool Worker::IsRunning(){
    // In case a simulation ended prematurely escaping the check routines in FrameMove (only executed after at least one second)
    bool state = simManager.GetRunningStatus();;
    if(!state && simuTimer.isActive)
        simuTimer.Stop();
    return state;
}

void Worker::Update(float appTime) {
    //Refreshes interface cache:
    //Global statistics, leak/hits, global histograms
    //Facet hits, facet histograms, facet angle maps
    //No cache for profiles, textures, directions (plotted directly from shared memory hit buffer)


    if (needsReload) RealReload();

    // End of simulation reached (Stop GUI)
    if ((simManager.allProcsDone || simManager.hasErrorStatus) && (IsRunning() || (!IsRunning() && simuTimer.isActive)) && (appTime != 0.0f)) {
        InnerStop(appTime);
        if (simManager.hasErrorStatus) ThrowSubProcError();
    }

    // Retrieve hit count recording from the shared memory
    // Globals
    if(!globState.initialized || !simManager.nbThreads)
        return;
    mApp->changedSinceSave = true;


#if defined(MOLFLOW)
    bool needsAngleMapStatusRefresh = false;
#endif

    //for(auto& simUnit : simManager.simUnits) {
        /*if (!simUnit->tMutex.try_lock_for(std::chrono::milliseconds (100))) {
            continue;
        }*/
        // Global hits and leaks
        {
            size_t waitTime = (this->simManager.isRunning) ? 100 : 10000;
            if (!globState.tMutex.try_lock_for(std::chrono::milliseconds(waitTime))) {
                return;
            }
        }

        globState.stateChanged = false;
        globalHitCache = globState.globalHits;

#if defined(SYNRAD)

        if (globalHitCache.globalHits.nbDesorbed && model->wp.nbTrajPoints) {
            no_scans = (double)globalHitCache.globalHits.nbDesorbed / (double)model->wp.nbTrajPoints;
        }
        else {
            no_scans = 1.0;
        }
#endif


        //Copy global histogram
        //Prepare vectors to receive data
        RetrieveHistogramCache();

    //buffer = bufferStart; //Commented out as not used anymore

        // Refresh local facet hit cache for the displayed moment
        size_t nbFacet = geom->GetNbFacet();
        for (size_t i = 0; i < nbFacet; i++) {
            InterfaceFacet *f = geom->GetFacet(i);
#if defined(SYNRAD)
            //memcpy(&(f->facetHitCache), buffer + f->sh.hitOffset, sizeof(FacetHitBuffer));
#endif
#if defined(MOLFLOW)
            if (f->sh.anglemapParams.record) { //Recording, so needs to be updated
                if (f->selected && f->angleMapCache.empty())
                    needsAngleMapStatusRefresh = true; //Will update facetadvparams panel
                //Retrieve angle map from hits dp
                if(f->sh.desorbType != DES_ANGLEMAP) {
                    /*model->facets[i]->angleMap.pdf = globState.facetStates[i].recordedAngleMapPdf;
                    f->angleMapCache = model->facets[i]->angleMap.pdf;*/
                    f->angleMapCache = globState.facetStates[i].recordedAngleMapPdf;
                }
            }
#endif
        }

        //simUnit->tMutex.unlock();

    try {
        if (mApp->needsTexture || mApp->needsDirection) {
            CalculateTextureLimits();
            geom->BuildFacetTextures(globState, mApp->needsTexture, mApp->needsDirection);
        }
    }
    catch (const std::exception &e) {
        globState.tMutex.unlock();
        GLMessageBox::Display(e.what(), "Error building texture", GLDLG_OK, GLDLG_ICONERROR);
        return;
    }
    globState.tMutex.unlock();

#if defined(MOLFLOW)
    if (mApp->facetAdvParams && mApp->facetAdvParams->IsVisible() && needsAngleMapStatusRefresh)
        mApp->facetAdvParams->Refresh(geom->GetSelectedFacets());
#endif
}

void Worker::GetProcStatus(size_t *states, std::vector<std::string> &statusStrings) {

    if (model->otfParams.nbProcess == 0) return;
    simManager.GetProcStatus(states, statusStrings);
}

void Worker::GetProcStatus(ProcComm &procInfoList) {
    simManager.GetProcStatus(procInfoList);
}

void Worker::ChangePriority(int prioLevel) {
    if(prioLevel)
        simManager.IncreasePriority();
    else
        simManager.DecreasePriority();
}

[[maybe_unused]] std::vector<std::vector<double>> Worker::ImportCSV_double(FileReader *file) {
    std::vector<std::vector<double>> table;
    do {
        std::vector<double> currentRow;
        do {
            currentRow.push_back(file->ReadDouble());
            if (!file->IsEol()) file->ReadKeyword(",");
        } while (!file->IsEol());
        table.push_back(currentRow);
    } while (!file->IsEof());
    return table;
}



void Worker::ChangeSimuParams() { //Send simulation mode changes to subprocesses without reloading the whole geometry
    if (model->otfParams.nbProcess == 0 || !geom->IsLoaded()) return;
    if (needsReload) RealReload(); //Sync (number of) regions

    //auto *progressDlg = new GLProgress("Creating dataport...", "Passing simulation mode to workers");
    //progressDlg->SetVisible(true);
    //progressDlg->SetProgress(0.0);

    //To do: only close if parameters changed
    //progressDlg->SetMessage("Waiting for subprocesses to release log dataport...");

    particleLog.clear();

    //progressDlg->SetProgress(0.5);
    //progressDlg->SetMessage("Assembling parameters to pass...");

    std::string loaderString = SerializeParamsForLoader().str();
    try {
        simManager.ForwardOtfParams(&model->otfParams);

        if(simManager.ShareWithSimUnits((BYTE *) loaderString.c_str(), loaderString.size(),LoadType::LOADPARAM)){
            auto errString = fmt::format("Failed to send params to sub process:\n");
            //GLMessageBox::Display(errString.c_str(), "Warning
            //
            // (Updateparams)", GLDLG_OK, GLDLG_ICONWARNING);

            //progressDlg->SetVisible(false);
            //SAFE_DELETE(progressDlg);
            return;
            //throw std::runtime_error(errString.c_str());
        }
    }
    catch (const std::exception &e) {
        //GLMessageBox::Display(e.what(), "Error (LoadGeom)", GLDLG_OK, GLDLG_ICONERROR);
    }

    //progressDlg->SetVisible(false);
    //SAFE_DELETE(progressDlg);

/*#if defined(SYNRAD)
    //Reset leak and hit cache
    ResetWorkerStats();
    Update(0.0f);
#endif*/
}

/**
* \brief Extract a 7z file and return the file handle
* \param fileName name of the input file
* \param geomName name of the geometry file
* \return handle to opened decompressed file
*/
FileReader *Worker::ExtractFrom7zAndOpen(const std::string &fileName, const std::string &geomName) {
    std::ostringstream cmd;
    std::string sevenZipName;

#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    //Necessary push/pop trick to support UNC (network) paths in Windows command-line
    auto CWD = FileUtils::get_working_path();
    cmd << "cmd /C \"pushd \"" << CWD << "\"&&";
    sevenZipName += "7za.exe";
#else //Linux, MacOS
    sevenZipName = "7za"; //so that Exist() check fails and we get an error message on the next command
    std::string possibleLocations[] = {"./7za", //use 7za binary shipped with Molflow
                                       "/usr/bin/7za", //use p7zip installed system-wide
                                       "/usr/local/bin/7za", //use p7zip installed for user
                                       "/opt/homebrew/bin/7za"}; //homebrew on M1 mac
    for(auto& path : possibleLocations){
        if (FileUtils::Exist(path)) {
            sevenZipName = path;
        }
    }
#endif

    if (!FileUtils::Exist(sevenZipName)) {
        throw Error("7-zip compressor not found, can't extract file.");
    }
    cmd << sevenZipName << " x -t7z -aoa \"" << fileName << "\" -otmp";

#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    cmd << "&&popd\"";
#endif
    system(cmd.str().c_str());

    std::string toOpen, prefix;
    std::string shortFileName = FileUtils::GetFilename(fileName);
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    prefix = CWD + "\\tmp\\";
#else
    prefix = "tmp/";
#endif
    toOpen = prefix + geomName;
    if (!FileUtils::Exist(toOpen))
        //Inside the zip, try original filename with extension changed from geo7z to geo
        toOpen = prefix + (shortFileName).substr(0,shortFileName.length() - 2);


    return new FileReader(toOpen); //decompressed file opened
}


/**
* \brief Function that updates the global hit counter with the cached value + releases the mutex
* Send total hit counts to subprocesses
*/
void Worker::SendToHitBuffer() {
    simManager.ForwardGlobalCounter(&globState, &particleLog);
}

/**
* \brief Saves current facet hit counter from cache to results
*/
void Worker::SendFacetHitCounts() {
    size_t nbFacet = geom->GetNbFacet();
    std::vector<FacetHitBuffer*> facetHitCaches;
    for (size_t i = 0; i < nbFacet; i++) {
        InterfaceFacet *f = geom->GetFacet(i);
        facetHitCaches.push_back(&f->facetHitCache);
    }
    simManager.ForwardFacetHitCounts(facetHitCaches);
}

void Worker::RetrieveHistogramCache()
{
	//Copy histograms from hit dataport to local cache
	//The hit dataport contains histograms for all moments, the cache only for the displayed
    //dpHitStartAddress is the beginning of the dpHit buffer

    //GLOBAL HISTOGRAMS
	//Prepare vectors to receive data
#if defined(MOLFLOW)
    globalHistogramCache = globState.globalHistograms[displayedMoment];
#endif
#if defined(SYNRAD)
    if(!globState.globalHistograms.empty())
        globalHistogramCache = globState.globalHistograms[0];
#endif
    //FACET HISTOGRAMS
    for (size_t i = 0;i < geom->GetNbFacet();i++) {
        InterfaceFacet* f = geom->GetFacet(i);
#if defined(MOLFLOW)
        f->facetHitCache = globState.facetStates[i].momentResults[displayedMoment].hits;
        f->facetHistogramCache = globState.facetStates[i].momentResults[displayedMoment].histogram;
#endif
#if defined(SYNRAD)
        f->facetHitCache = globState.facetStates[i].momentResults[0].hits;
        f->facetHistogramCache = globState.facetStates[i].momentResults[0].histogram;
#endif
    }
}

// returns -1 on error, 0 on success
int Worker::ReloadSim(bool sendOnly, GLProgress *progressDlg) {
    // Send and Load geometry
    progressDlg->SetMessage("Waiting for subprocesses to load geometry...");
    try {
        if (!InterfaceGeomToSimModel()) {
            std::string errString = "Failed to send geometry to sub process!\n";
            /*GLMessageBox::Display(errString.c_str(), "Warning (LoadGeom)", GLDLG_OK, GLDLG_ICONWARNING);

            progressDlg->SetVisible(false);
            SAFE_DELETE(progressDlg);*/
            throw std::runtime_error(errString.c_str());
        }

        progressDlg->SetMessage("Initialising physical properties for model->..");
        if(model->PrepareToRun()){
            throw std::runtime_error("Error preparing simulation");
        }

        progressDlg->SetMessage("Constructing memory structure to store results...");
        if (!sendOnly) {
            globState.Resize(model);
        }

        progressDlg->SetMessage("Forwarding simulation model...");
        simManager.ForwardSimModel(model);
        progressDlg->SetMessage("Forwarding global simulation state...");
        simManager.ForwardGlobalCounter(&globState, &particleLog);
    }
    catch (const std::exception &e) {
        GLMessageBox::Display(e.what(), "Error (LoadGeom)", GLDLG_OK, GLDLG_ICONERROR);
        throw;
    }

    return 0;
}