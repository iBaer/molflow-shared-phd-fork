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

#ifndef MOLFLOW_PROJ_IMGUIWINDOW_H
#define MOLFLOW_PROJ_IMGUIWINDOW_H


#include <GLApp/GLApp.h>
#include <imgui/imgui.h>
#include "AppUpdater.h"

#if defined(MOLFLOW)
#include "../../src/MolFlow.h"
#else
#include "../../src/SynRad.h"
#endif

/*! Window manager for the Imgui GUI, right now it is rendered above the SDL2 GUI */
class ImguiWindow {
public:
    explicit ImguiWindow(GLApplication* app) {this->app = app;};
    void init();
    static void destruct();
    //void render();
    void renderSingle(); //!< Main function, calling a frame rendering cycle handling all other Imgui windows

    GLApplication* app;

    bool ToggleMainHub();
    bool ToggleMainMenu();
    bool ToggleSimSidebar();
    bool ToggleDemoWindow();
    bool ToggleGlobalSettings();
    bool ToggleAABBWindow();
    bool ToggleConvWindow();
protected:
    // Window states (visible or not)
    bool show_main_hub{false}; //!< Hub managing all other windows
    bool show_app_main_menu_bar{false}; //!< Window main menu bar at the top
    bool show_app_sidebar{false}; //!< Simulation sidebar for various 3d viewer, facet and simulation settings
    bool show_demo_window{false}; //!< Debug only: ImGui Demo Window to test all ImGui functionalities
    bool show_global_settings{false}; //!< Global Settings window
    bool show_perfo{false}; //!< Plot showing history of simulation performance
    bool show_conv{false}; //!< Plot showing convergence of simulation stats
    bool show_aabb{false};
    bool show_select{false};

    double start_time; // to keep track how long the ImGui GUI is running
#if defined(MOLFLOW)
    static void restartProc(int nbProc, MolFlow *mApp);
#else
    static void restartProc(int nbProc, SynRad *mApp);
#endif
};


#endif //MOLFLOW_PROJ_IMGUIWINDOW_H
