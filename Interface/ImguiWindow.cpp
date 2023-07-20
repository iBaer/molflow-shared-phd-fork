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

#include "ImguiExtensions.h"
#include "ImguiWindow.h"
#include "ImguiMenu.h"
#include "ImguiAABB.h"

#if defined(MOLFLOW)
#include "../../src/MolflowGeometry.h"
#else
#include "../../src/SynradGeometry.h"
#endif
#include "Facet_shared.h"
#include "../../src/Interface/Viewer3DSettings.h"

#include "imgui/imgui.h"
#include "imgui/imgui_impl_opengl2.h"
#include "imgui/imgui_impl_sdl.h"
#include "ImguiGlobalSettings.h"
#include "ImguiShowNeighbors.h"
#include "ImguiGeomAnalysis.h"
#include "ImguiPerformancePlot.h"
#include "ImguiSidebar.h"
#include "ImguiConvergencePlot.h"

#include <imgui/imgui_internal.h>
#include <imgui/IconsFontAwesome5.h>
#include <future>
#include <implot/implot.h>
#include <Helper/FormatHelper.h>
#include <GeometryTools.h>

// Varius toggle functions for individual window components
bool ImguiWindow::ToggleMainHub(){
    show_main_hub = !show_main_hub;
    return show_main_hub;
}
bool ImguiWindow::ToggleMainMenu(){
    show_app_main_menu_bar = !show_app_main_menu_bar;
    return show_app_main_menu_bar;
}
bool ImguiWindow::ToggleSimSidebar(){
    show_app_sidebar = !show_app_sidebar;
    return show_app_sidebar;
}
bool ImguiWindow::ToggleDemoWindow(){
    show_demo_window = !show_demo_window;
    return show_demo_window;
}
bool ImguiWindow::ToggleGlobalSettings(){
    show_global_settings = !show_global_settings;
    return show_global_settings;
}
bool ImguiWindow::ToggleAABBWindow(){
    show_aabb = !show_aabb;
    return show_aabb;
}
bool ImguiWindow::ToggleConvWindow(){
    show_conv = !show_conv;
    return show_conv;
}

bool ImguiWindow::ToggleGeomAnalysis(){
    show_analysis = !show_analysis;
    return show_analysis;
}


// --- Toggle functions ---

// Setup Dear ImGui context and various default values (font, colors etc.)
void ImguiWindow::init() {
    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO &io = ImGui::GetIO();
    (void) io;

    // Enable Multi Viewport
    //ImGui::GetIO().ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;

    // io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable
    // Keyboard Controls io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad; //
    // Enable Gamepad Controls

    // Setup Dear ImGui style
    // ImGui::StyleColorsDark();
    // ImGui::StyleColorsClassic();
    ImGui::StyleColorsLight();

    // Setup Platform/Renderer backends
    ImGui_ImplSDL2_InitForOpenGL(app->mainScreen, app->mainContext);
    ImGui_ImplOpenGL2_Init();

    // Load Fonts
    // - If no fonts are loaded, dear imgui will use the default font. You can
    // also load multiple fonts and use ImGui::PushFont()/PopFont() to select
    // them.
    // - AddFontFromFileTTF() will return the ImFont* so you can store it if you
    // need to select the font among multiple.
    // - If the file cannot be loaded, the function will return nullptr. Please
    // handle those errors in your application (e.g. use an assertion, or display
    // an error and quit).
    // - The fonts will be rasterized at a given size (w/ oversampling) and stored
    // into a texture when calling ImFontAtlas::Build()/GetTexDataAsXXXX(), which
    // ImGui_ImplXXXX_NewFrame below will call.
    // - Read 'docs/FONTS.md' for more instructions and details.
    // - Remember that in C/C++ if you want to include a backslash \ in a string
    // literal you need to write a double backslash \\ !
    // io.Fonts->AddFontDefault();
    //static const ImWchar sym_ranges[] = {0x2000, 0x3000, 0}; // some symbols
    static const ImWchar sym_ranges[] = {0x0000, 0xFFFF, 0}; // epsilon etc
    ImFontConfig sym_config;
    sym_config.MergeMode = true;
    sym_config.PixelSnapH = true;
    io.Fonts->AddFontFromFileTTF("DroidSans.ttf", 16.0f);
    io.Fonts->AddFontFromFileTTF("FreeMono.ttf", 16.0f, &sym_config, sym_ranges); // vector arrow

    // merge in icons from Font Awesome
    static const ImWchar icons_ranges[] = {ICON_MIN_FA, ICON_MAX_FA, 0};
    ImFontConfig icons_config;
    icons_config.MergeMode = true;
    icons_config.PixelSnapH = true;
    io.Fonts->AddFontFromFileTTF(FONT_ICON_FILE_NAME_FAS, 16.0f, &icons_config, icons_ranges);

    io.Fonts->AddFontFromFileTTF("DroidSans.ttf", 14.0f);
    io.Fonts->AddFontFromFileTTF("FreeMono.ttf", 14.0f, &sym_config, sym_ranges); // vector arrow
    io.Fonts->AddFontFromFileTTF(FONT_ICON_FILE_NAME_FAS, 14.0f, &icons_config, icons_ranges);

    /*io.Fonts->AddFontFromFileTTF("FreeMono.ttf", 16.0f);
    io.Fonts->AddFontFromFileTTF("FreeMono.ttf", 16.0f, &sym_config, sym_ranges);*/

// use FONT_ICON_FILE_NAME_FAR if you want regular instead of solid

    // io.Fonts->AddFontFromFileTTF("../../misc/fonts/Cousine-Regular.ttf", 15.0f);
    // io.Fonts->AddFontFromFileTTF("../../misc/fonts/DroidSans.ttf", 16.0f);
    // io.Fonts->AddFontFromFileTTF("../../misc/fonts/ProggyTiny.ttf", 10.0f);
    // ImFont* font =
    // io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\ArialUni.ttf", 18.0f,
    // nullptr, io.Fonts->GetGlyphRangesJapanese()); IM_ASSERT(font != nullptr);

    show_main_hub = false;
    show_demo_window = false;
    show_global_settings = false;
    show_app_main_menu_bar = false;
    show_app_sidebar = false;
    show_aabb = false;
    show_perfo = false;
    show_conv = false;
    show_select = false;
    show_analysis = false;

    start_time = ImGui::GetTime();
}

// Gracefully clears and shutsdown Dear ImGui context
void ImguiWindow::destruct() {
    // Cleanup
    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplSDL2_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
}

// TODO: When Imgui will be the main window/main GUI, use a full render cycle instead of a single frame rendering call
/*void ImguiWindow::render() {
  ImGuiIO &io = ImGui::GetIO();
  (void)io;

  // Our state
  bool show_demo_window = false;
  bool show_global_settings = true;
  ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

  // Main loop
  bool done = false;
  while (!done) {
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
      ImGui_ImplSDL2_ProcessEvent(&event);
      if (event.type == SDL_QUIT)
        done = true;
      if (event.type == SDL_WINDOWEVENT &&
          event.window.event == SDL_WINDOWEVENT_CLOSE &&
          event.window.windowID == SDL_GetWindowID(app->mainScreen))
        done = true;
    }

    // Start the Dear ImGui frame
    ImGui_ImplOpenGL2_NewFrame();
    ImGui_ImplSDL2_NewFrame(app->mainScreen);
    ImGui::NewFrame();

    // Show contents

    // Rendering
    ImGui::Render();
    glViewport((int)io.DisplaySize.x * 0.25, (int)io.DisplaySize.y * 0.25,
               (int)io.DisplaySize.x / 2.0, (int)io.DisplaySize.y / 2.0);
    // glClearColor(clear_color.x * clear_color.w, clear_color.y *
    // clear_color.w, clear_color.z * clear_color.w, clear_color.w);
    // glClear(GL_COLOR_BUFFER_BIT);
    // glUseProgram(0); // You may want this if using this code in an OpenGL 3+
    // context where shaders may be bound
    ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
    SDL_GL_SwapWindow(app->mainScreen);
  }
}*/

// Function for an individual frame rendering step
// If active, renders the individual components and handles the corresponding user actions
void ImguiWindow::renderSingle() {
    ImGuiIO &io = ImGui::GetIO();
    (void) io;

#if defined(MOLFLOW)
    MolFlow *mApp = (MolFlow *) app;
#else
    SynRad *mApp = (SynRad *) app;
#endif
    if (mApp) {
        bool nbProcChanged = false;
        bool recalcOutg = false;
        bool changeDesLimit = false;
        bool redrawAabb = false;
        bool rebuildAabb = false;
        static int nbProc = mApp->worker.GetProcNumber();


        // Start the Dear ImGui frame
        ImGui_ImplOpenGL2_NewFrame();
        ImGui_ImplSDL2_NewFrame(/*app->mainScreen*/);
        ImGui::NewFrame();

        if (show_app_main_menu_bar)
            ShowAppMainMenuBar();

        static bool open_viewer_window = false;
        if (show_app_sidebar)
            ShowAppSidebar(&show_app_sidebar, mApp, mApp->worker.GetGeometry(), &show_global_settings,
                           &open_viewer_window);

        // 1. Show the big demo window (Most of the sample code is in
        // ImGui::ShowDemoWindow()! You can browse its code to learn more about Dear
        // ImGui!).
        if ((io.KeyCtrl && io.KeyShift && io.KeyAlt && ImGui::IsKeyDown(ImGuiKey_D)))
            show_demo_window = !show_demo_window;
        if (show_demo_window) {
            ImGui::ShowDemoWindow(&show_demo_window);
            ImPlot::ShowDemoWindow(&show_demo_window);
        }

        // 2. Show Molflow x ImGui Hub window
        if (show_main_hub) {
            ImGui::Begin("[BETA] _Molflow ImGui Suite_"); // Create a window called "Hello, world!"
            // and append into it.

#if defined(DEBUG)
// only show in debug mode
            ImGui::Checkbox(
                    "Demo Window",
                    &show_demo_window); // Edit bools storing our window open/close state
#endif
            ImGui::Checkbox("Global settings", &show_global_settings);
            ImGui::Checkbox("Acceleration Structure", &show_aabb);
            ImGui::Checkbox("Menu bar", &show_app_main_menu_bar);
            ImGui::Checkbox("Sidebar", &show_app_sidebar);
            ImGui::Checkbox("Performance Plot", &show_perfo);
            ImGui::Checkbox("Convergence Plot", &show_conv);
            ImGui::Checkbox("Select neighbors", &show_select);
            ImGui::Checkbox("Geometrie analysis", &show_analysis);

            ImGui::Text("Avg %.3f ms/frame (%.1f FPS)",
                        1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
            auto now_time = ImGui::GetTime();
            ImGui::Text("Application time %.3f s [%.3f s]",
                        ImGui::GetTime(), difftime(now_time, start_time));
            ImGui::End();
        }

        // 3. Show window plotting the simulation performance
        if (show_perfo) {
            ShowPerfoPlot(&show_perfo, mApp);
        }

        // 3. Show window plotting the simulation convergence
        if (show_conv) {
            ShowConvPlot(&show_conv, mApp);
        }

        // 3. Show global settings
        if (show_global_settings) {
            ShowGlobalSettings(mApp, &show_global_settings, nbProcChanged, recalcOutg, changeDesLimit, nbProc);
            ImGui::End();
        }

        // 4. Create placeholder for async launches
        std::promise<int> p;
        static std::future<int> future_int;
        static bool active_prev_state;
        if (!future_int.valid())
            future_int = p.get_future();

        if (show_select) {
            ShowNeighborSelect(mApp, &show_select);
        }

        if (show_analysis) {
            ShowGeomAnalysis(mApp, &show_analysis);
        }

        // 5. Show Window for ADS configuration/visualisation
        if (show_aabb) {
            static ImguiAABBVisu visu{};
            visu.ShowAABB(mApp, &show_aabb, redrawAabb, rebuildAabb);

            if (rebuildAabb) {
                if (mApp->worker.IsRunning())
                    mApp->StartStopSimulation();

                future_int = std::async(std::launch::async, &SimulationModel::BuildAccelStructure, mApp->worker.model,
                                        &mApp->worker.globState, mApp->worker.model->wp.accel_type,
                                        mApp->worker.model->wp.splitMethod, mApp->worker.model->wp.bvhMaxPrimsInNode, mApp->worker.model->wp.hybridWeight);

                active_prev_state = true;
                mApp->wereEvents = true;
            }
            // Evaluate running progress
            if (future_int.wait_for(std::chrono::seconds(0)) != std::future_status::ready) {
                static float progress = 0.0f, load_time = 0.0f;
                ImGui::Loader(progress, load_time);

                active_prev_state = true;
                mApp->wereEvents_imgui = true;
            }
            else if (active_prev_state) {
                redrawAabb = true;
                active_prev_state = false;
                mApp->wereEvents_imgui = true;
            }
        }

        // Rendering
        ImGui::Render();
        glViewport(0, 0, (int) io.DisplaySize.x, (int) io.DisplaySize.y);
        ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
        // SDL_GL_SwapWindow(app->mainScreen);

        // Update and Render additional Platform Windows
        /*if (io.ConfigFlags & ImGuiConfigFlags_ViewportsEnable)
        {
            ImGui::UpdatePlatformWindows();
            ImGui::RenderPlatformWindowsDefault();
        }*/

        // Handle button events at the end, because some functions call GL Repaint,
        // which doesnt work well with ImGui if frame is not finalized
        if (nbProcChanged) {
            restartProc(nbProc, mApp);
        } else if (recalcOutg) {
            if (mApp->AskToReset()) {
                try {
                    mApp->worker.RealReload();
                } catch (std::exception &e) {
                    if (ImGui::BeginPopupModal("Error", nullptr,
                                               ImGuiWindowFlags_AlwaysAutoResize)) {
                        ImGui::Text("Recalculation failed: Couldn't reload Worker:\n%s",
                                    e.what());
                        ImGui::Separator();

                        if (ImGui::Button("OK", ImVec2(120, 0))) {
                            ImGui::CloseCurrentPopup();
                        }
                        ImGui::SetItemDefaultFocus();
                        ImGui::EndPopup();
                    }
                }
            }
        } else if (changeDesLimit) {
            mApp->worker.ChangeSimuParams(); // Sync with subprocesses
        } else if (redrawAabb) {
            mApp->worker.GetGeometry()->BuildGLList();
            redrawAabb = false;
        }
        if(open_viewer_window){
            open_viewer_window = false;
            if (!mApp->viewer3DSettings)
                mApp->viewer3DSettings = new Viewer3DSettings();
            mApp->viewer3DSettings->SetVisible(!mApp->viewer3DSettings->IsVisible());
            mApp->viewer3DSettings->Reposition();
            auto curViewer = mApp->curViewer;
            auto viewer = mApp->viewer[curViewer];
            mApp->viewer3DSettings->Refresh(mApp->worker.GetGeometry(), viewer);
        }
    }
}

/**
 * \brief Function to apply changes to the number of processes.
 */
#if defined(MOLFLOW)
void ImguiWindow::restartProc(int nbProc, MolFlow *mApp) {
#else
void ImguiWindow::restartProc(int nbProc, SynRad *mApp) {
#endif
    try {
        mApp->worker.Stop_Public();
        mApp->worker.SetProcNumber(nbProc);
        mApp->worker.RealReload(true);
        mApp->SaveConfig();
    } catch (Error &e) {
        if (ImGui::BeginPopupModal("Error", nullptr,
                                   ImGuiWindowFlags_AlwaysAutoResize)) {
            ImGui::Text("%s", e.what());
            ImGui::Separator();

            if (ImGui::Button("OK", ImVec2(120, 0))) {
                ImGui::CloseCurrentPopup();
            }
            ImGui::SetItemDefaultFocus();
            ImGui::EndPopup();
        }
    }
}