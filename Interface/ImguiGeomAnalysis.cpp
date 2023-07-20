//
// Created by pbahr on 11/9/21.
//

#include <vector>
#include <imgui.h>
#include <GeometryTools.h>
#include "ImguiGeomAnalysis.h"
#include "fmt/ranges.h"

#if defined(MOLFLOW)
#include "../../src/MolFlow.h"
#endif

#if defined(SYNRAD)
#include "../src/SynRad.h"
#endif

void ShowGeomAnalysis(MolFlow *mApp, bool *show_analysis){
    ImGui::Begin("Geometry analysis", show_analysis); // Create a window called "Hello, world!"
    // and append into it.

    if(ImGui::Button("Use all functions")){
        GeometryTools::AnalyseGeometry(mApp->worker.GetGeometry());
    }

    static int func = 2;
    static bool run_cleanup = false;
    ImGui::Checkbox("Apply cleanup for common edges", &run_cleanup);
    ImGui::SliderInt("Select function", &func, 0, 5);

    if(ImGui::Button("Use selected function")){
        GeometryTools::AnalyseGeometry(mApp->worker.GetGeometry(), func, run_cleanup);
    }

    static int comp_func = 2;
    ImGui::SliderInt("Compare against", &comp_func, 0, 5);

    if(ImGui::Button("Compare all functions")){
        GeometryTools::CompareAlgos(mApp->worker.GetGeometry(), comp_func);
    }

    if(ImGui::Button("Compare selected function")){
        GeometryTools::CompareAlgorithm(mApp->worker.GetGeometry(), func, comp_func);
    }

    ImGui::End();
}