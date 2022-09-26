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

#include "ImguiConvergencePlot.h"
#include "imgui/imgui.h"
#include <implot/implot.h>

#include "Interface.h"

// Utility classes from implot_demo.cpp
// utility structure for realtime plot
struct ScrollingBuffer {
    int MaxSize;
    int Offset;
    ImVector<ImVec2> Data;
    ScrollingBuffer(int max_size = 20) {
        MaxSize = max_size;
        Offset  = 0;
        Data.reserve(MaxSize);
    }
    void AddPoint(float x, float y) {
        if (Data.size() < MaxSize)
            Data.push_back(ImVec2(x,y));
        else {
            Data[Offset] = ImVec2(x,y);
            Offset =  (Offset + 1) % MaxSize;
        }
    }
    void Erase() {
        if (Data.size() > 0) {
            Data.shrink(0);
            Offset  = 0;
        }
    }
};

// Demonstrate creating a simple static window with no decoration
// + a context-menu to choose which corner of the screen to use.
void ShowConvPlot(bool *p_open, Interface *mApp) {
    ImGuiIO &io = ImGui::GetIO();

    // Always center this window when appearing
    ImVec2 center = ImGui::GetMainViewport()->GetCenter();
    ImGui::SetNextWindowPos(center, ImGuiCond_Appearing,
                            ImVec2(0.5f, 0.5f));
    ImGui::SetNextWindowSize(ImVec2(400, 400), ImGuiCond_FirstUseEver);

    static ImGuiWindowFlags flags =/*
            ImGuiWindowFlags_AlwaysAutoResize |*/
            ImGuiWindowFlags_NoSavedSettings;

    if (ImGui::Begin("Convergence Plot", p_open, flags)) {

        // Fill an array of contiguous float values to plot
        // Tip: If your float aren't contiguous but part of a structure, you can pass a pointer to your first float
        // and the sizeof() of your structure in the "stride" parameter.
        static ScrollingBuffer values(100);
        //static ScrollingBuffer values_des;
        //static ScrollingBuffer tvalues;

        static auto refresh_time = ImGui::GetTime();
        if (!true || refresh_time == 0.0) // force
            refresh_time = ImGui::GetTime();
        auto now_time = ImGui::GetTime();

        auto formulas = mApp->formula_ptr;
        auto formula_id = 0;
        if(!formulas->formulas_n.empty()) {
            if (mApp->worker.IsRunning() && difftime(now_time, refresh_time) > 1.0) {
                //static float phase = 0.0f;

                if (mApp->worker.globalHitCache.globalHits.nbDesorbed > 0) {
                    auto& conv_vec = formulas->convergenceValues[formula_id].conv_vec;

                    if(!conv_vec.empty()) {
                        auto &conv_last = conv_vec.back();
                        if (values.Data.Size == 0 || conv_last.first >= values.Data[std::max(0,(values.Offset - 1) % values.Data.Size)].x) {
                            //for (int j = std::max(0,(int)conv_vec.size()-1000); j < conv_vec.size(); j++) {// limit data points to last 1000
                            values.AddPoint(conv_last.first, conv_last.second);
                        }
                    }
                }

                //phase += 0.10f * values_offset;
                refresh_time = now_time;
            }

            // Plots can display overlay texts
            // (in this example, we will display an average value)
            {

                float max_val = -9999999.9f;
                float min_val = 9999999.9f;
                for (int i = 0; i < values.Data.size(); ++i) {
                    if (values.Data[i].y > max_val) {
                        max_val = values.Data[i].y;
                    }
                    if (values.Data[i].y < min_val) {
                        min_val = values.Data[i].y;
                    }
                    /*if (values_des[i] > max_val) {
                        max_val = values_des[i];
                    }
                    if (values_des[i] < min_val) {
                        min_val = values_des[i];
                    }*/
                }
                /*char overlay[32];
                sprintf(overlay, "avg %f hit/s", average);*/
                //ImGui::PlotLines(""*//*"Hit/s"*//*, values, IM_ARRAYSIZE(values), values_offset, overlay, min_val * 0.95f, max_val * 1.05f,ImVec2(0, 80.0f));

                float rel = (max_val - min_val) / (max_val + min_val);
                ImPlot::SetNextAxisLimits(ImAxis_Y1, std::max(0.0f, min_val * (1.0f-0.2f*rel)), max_val * (1.0f+0.2f*rel), ImGuiCond_Always);
                if (ImPlot::BeginPlot("##Conv", "Numer of desorptions", "Value (formula)", ImVec2(-1, -1),
                                      ImPlotAxisFlags_None,
                                      ImPlotAxisFlags_AutoFit /*| ImPlotAxisFlags_Time*//*, ImPlotAxisFlags_AutoFit*/)) {


                    ImPlot::PushStyleVar(ImPlotStyleVar_FillAlpha, 0.25f);
                    if(values.Data.Size > 1) {
                        ImPlot::PlotLine("Hit/s", &values.Data[0].x, &values.Data[0].y, values.Data.size(), 0,
                                         values.Offset, 2 * sizeof(float));
                        ImPlot::PlotShaded("Hit/s", &values.Data[0].x, &values.Data[0].y, values.Data.size(), -INFINITY, 0, values.Offset, 2 * sizeof(float));
                    }
                    //ImPlot::PlotLine("Des/s", tvalues, values_des, IM_ARRAYSIZE(values_des), values_offset);
                    //ImPlot::PlotShaded("Des/s", tvalues, values_des, IM_ARRAYSIZE(values_des), -INFINITY, values_offset);
                    ImPlot::PopStyleVar();

                    ImPlot::EndPlot();
                }

                if(ImGui::Button("Reset Data")){
                    values = {};
                }


                static int val_ind = 0;
                if(values.Data.Size >= 1) {
                    ImGui::SliderInt("Index", &val_ind, 0, values.Data.size());
                    ImGui::InputFloat2("##", &values.Data[val_ind].x);
                }
            }
        }
    }
    ImGui::End();
}
