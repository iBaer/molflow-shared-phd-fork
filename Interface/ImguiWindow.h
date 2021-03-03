//
// Created by pascal on 2/26/21.
//

#ifndef MOLFLOW_PROJ_IMGUIWINDOW_H
#define MOLFLOW_PROJ_IMGUIWINDOW_H


#include <GLApp/GLApp.h>
#include <imgui/imgui.h>
#include "AppUpdater.h"
#include "../../src/MolFlow.h"

class ImguiWindow {
public:
    ImguiWindow(GLApplication* app) {this->app = app;};
    void init();
    void destruct();
    void render();
    void renderSingle();

    GLApplication* app;

protected:
    // Our state
    bool show_demo_window;
    bool show_another_window;
    ImVec4 clear_color;

    static void restartProc(int nbProc, MolFlow *mApp);

    static void ProcessControllTable(MolFlow *mApp);

    //static bool InputRightSide(const char *desc, double *const &val, const char *format);
};


#endif //MOLFLOW_PROJ_IMGUIWINDOW_H
