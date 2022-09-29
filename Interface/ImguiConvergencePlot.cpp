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
        if (!Data.empty()) {
            Data.shrink(0);
            Offset  = 0;
        }
    }
};

// Demonstrate creating a simple static window with no decoration
// + a context-menu to choose which corner of the screen to use.
void ShowConvPlot(bool *p_open, Interface *mApp) {
    ImGuiIO &io = ImGui::GetIO();

    static ImGuiWindowFlags flags =/*
            ImGuiWindowFlags_AlwaysAutoResize |*/
            ImGuiWindowFlags_NoSavedSettings;

    auto formulas = mApp->formula_ptr;

    static int c_prec = formulas->epsilon;
    static int cb_len = formulas->cb_length;
    static float c_rate = 0.0;
    static float sp1 = 0.0;
    static float sp2 = 0.0;
    static bool use_abs_err = formulas->useAbsEps;

    static bool plot_ascbr = true;
    static int plot_limit_y = 0;
    /*static bool plot_limit_y = true; // 1
    static bool plot_limit_y2 = false; // 2
    static bool plot_limit_y_ascbr = false; // 3*/

    if (ImGui::Begin("Convergence Settings", p_open, flags | ImGuiWindowFlags_AlwaysAutoResize)) {
        if(ImGui::BeginTabBar("Conv Crit")){
            if(ImGui::BeginTabItem("Plot options")) {
                ImGui::RadioButton("Unlimit Y", &plot_limit_y, 0);
                ImGui::RadioButton("Limit Y", &plot_limit_y, 1);
                ImGui::RadioButton("Limit Y 200%", &plot_limit_y,2);
                ImGui::RadioButton("Limit Y to ASCBR", &plot_limit_y,3);

                ImGui::EndTabItem();
            }
            if(ImGui::BeginTabItem("ASCBR")) {

                ImGui::Checkbox("Plot ASCBR Intervals", &plot_ascbr);
                ImGui::Spacing();

                bool apply_change = false;
                if(ImGui::SliderInt("Precision {10^-d}: d=", &c_prec, 0, 31)
                ||ImGui::SliderInt("CB length: C=", &cb_len, 0, 31)){
                    apply_change = true;
                }

               if(ImGui::Checkbox("Absolute precision", &use_abs_err)){
                    apply_change = true;
                }

               if(!formulas->formulas_n.empty()) {
                   ImGui::Spacing();
                   static int fid = 0;
                   ImGui::SliderInt("Formula #", &fid, 0, formulas->formulas_n.size()-1);
                   if (mApp->worker.IsRunning() && formulas->FetchNewConvValue()) {
                       sp1 = formulas->ApproxShapeParameter(fid, 1);
                       sp2 = formulas->ApproxShapeParameter(fid, 0);
                       c_rate = formulas->GetConvRate(fid);
                   }

                   ImGui::Text("%s", fmt::format("Conv. rate: {:e}", c_rate).c_str());
                   ImGui::SameLine();
                   ImGui::Text("%s", fmt::format("Shape Param 1: {:5.4f}", sp1).c_str());
                   ImGui::SameLine();
                   ImGui::Text("%s", fmt::format("Shape Param 2: {:5.4f}", sp2).c_str());

                   if (!formulas->convergenceValues.empty()) {
                       auto &formula = formulas->convergenceValues[fid];
                       ImGui::Text("%s", fmt::format("Current chain: {}", formula.chain_length).c_str());
                       ImGui::SameLine();
                       ImGui::Text("%s", fmt::format("Boundaries: [{:e} , {:e}]", formula.lower_bound,
                                                     formula.upper_bound).c_str());
                   }
               }// active formula

                if(apply_change){

                    formulas->useAbsEps = use_abs_err;
                    formulas->cb_length = cb_len;
                    formulas->epsilon = c_prec;

                    for(int id = 0; id < formulas->convergenceValues.size(); ++id)
                        formulas->RestartASCBR(id);
                }
                ImGui::EndTabItem();
            }
            ImGui::EndTabBar();
        } //tabbar
        ImGui::End();
    }

    // Always center this window when appearing
    ImVec2 center = ImGui::GetMainViewport()->GetCenter();
    ImGui::SetNextWindowPos(center, ImGuiCond_Appearing,
                            ImVec2(0.5f, 0.5f));
    ImGui::SetNextWindowSize(ImVec2(400, 400), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Convergence Plot", p_open, flags)) {

        // Fill an array of contiguous float values to plot
        // Tip: If your float aren't contiguous but part of a structure, you can pass a pointer to your first float
        // and the sizeof() of your structure in the "stride" parameter.
        static ScrollingBuffer values(1000);
        static ScrollingBuffer sqrtval1(101);
        static ScrollingBuffer sqrtval2(101);

        static ScrollingBuffer ascbr_lower(2);
        static ScrollingBuffer ascbr_upper(2);

        //static ScrollingBuffer values_des;
        //static ScrollingBuffer tvalues;

        static auto refresh_time = ImGui::GetTime();
        if (!true || refresh_time == 0.0) // force
            refresh_time = ImGui::GetTime();
        auto now_time = ImGui::GetTime();

        auto formula_id = 0;
        if(!formulas->formulas_n.empty()) {
            if (mApp->worker.IsRunning() && difftime(now_time, refresh_time) > 1.0) {
                //static float phase = 0.0f;

                if (mApp->worker.globalHitCache.globalHits.nbDesorbed > 0) {
                    auto& conv_vec = formulas->convergenceValues[formula_id].conv_vec;

                    if(!conv_vec.empty()) {

                        auto &conv_last = conv_vec.back();
                        if (values.Data.Size == 0 || (static_cast<float>(conv_last.first) > 0 && static_cast<float>(conv_last.first) >= values.Data[std::max(0, (values.Offset - 1) % values.MaxSize)].x)) {
                            //for (int j = std::max(0,(int)conv_vec.size()-1000); j < conv_vec.size(); j++) {// limit data points to last 1000
                            values.AddPoint(static_cast<float>(conv_last.first), static_cast<float>(conv_last.second));
                        }

                        if (values.Data.Size > 12) {
                            int first = std::max(0, (values.Offset) % values.Data.Size);
                            int mid = std::max(0, (first + 10) % values.Data.Size);
                            int last = values.Data.Size == values.MaxSize ? std::max(0, (values.Offset - 2) % values.Data.Size) : values.Data.Size - 1;

                            float c_1 = sqrtf(values.Data[first].x) * values.Data[mid].y;
                            float c_2 = sqrtf(values.Data[last].x) * values.Data[last].y;
                            float d_1 = values.Data[last].y - c_1 / sqrtf(values.Data[last].x);
                            float d_2 = values.Data[last].y + c_1 / sqrtf(values.Data[last].x);

                            sqrtval1.Erase();
                            sqrtval2.Erase();
                            for(int i = 0; i <= 100; i++) {
                                int n_des = values.Data[first].x + i * 0.01 *(values.Data[last].x - values.Data[first].x);
                                sqrtval1.AddPoint(
                                        static_cast<float>(n_des),
                                        static_cast<float>(c_1 / sqrtf(n_des) + d_1));
                                sqrtval2.AddPoint(
                                        static_cast<float>(n_des),
                                        static_cast<float>(-1.0 * c_1 / sqrtf(n_des) + d_2));

                            }
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
                for (int i = values.Data.size() > 20 ? 5 : 0; i < values.Data.size(); ++i) {
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
                //ImPlot::SetNextAxisLimits(ImAxis_Y1, std::max(0.0f, min_val * (1.0f-0.2f*rel)), max_val * (1.0f+0.2f*rel), ImGuiCond_Always);
                /*if (ImPlot::BeginPlot("##Conv", "Numer of desorptions", "Value (formula)", ImVec2(-1, -1),
                                      ImPlotAxisFlags_None,
                                      ImPlotAxisFlags_AutoFit *//*| ImPlotAxisFlags_Time*//**//*, ImPlotAxisFlags_AutoFit*//*)) {*/
                if (ImPlot::BeginPlot("##Conv", ImVec2(-1, -1), ImPlotFlags_None)) {
                    ImPlot::SetupAxes("Number of desorptions", "Value (formula)", ImPlotAxisFlags_AutoFit, plot_limit_y ? ImPlotAxisFlags_AutoFit : ImPlotAxisFlags_None);
                    //ImPlot::SetupAxesLimits(0, 0, std::max(0.0f, min_val * (1.0f-0.2f*rel)), max_val * (1.0f+0.2f*rel), ImPlotCond_Once);

                    ImPlot::PushStyleVar(ImPlotStyleVar_FillAlpha, 0.25f);
                    if(values.Data.Size > 1) {
                        float xmean = values.Data.back().y;
                        float dist_to_mean = std::max(fabs(min_val - xmean), fabs(max_val - xmean));
                        if(plot_limit_y == 1)
                            ImPlot::SetupAxisLimits(ImAxis_Y1, std::max(0.0f, min_val * (1.0f-0.2f*rel)), max_val * (1.0f+0.2f*rel), ImPlotCond_Always);
                        else if(plot_limit_y == 2)
                            ImPlot::SetupAxisLimits(ImAxis_Y1, xmean - 2.0f * dist_to_mean, xmean +  2.0f * dist_to_mean, ImPlotCond_Always);
                        else if(plot_ascbr && plot_limit_y == 3)
                            ImPlot::SetupAxisLimits(ImAxis_Y1, ascbr_lower.Data[0].y - dist_to_mean, ascbr_upper.Data[0].y + dist_to_mean, ImPlotCond_Always);

                        ImPlot::PlotLine(formulas->formulas_n.front()->GetName(), &values.Data[0].x, &values.Data[0].y, values.Data.size(), 0,
                                         values.Offset, 2 * sizeof(float));
                        ImPlot::PlotShaded(formulas->formulas_n.front()->GetName(), &values.Data[0].x, &values.Data[0].y, values.Data.size(), -INFINITY, 0, values.Offset, 2 * sizeof(float));
                    }
                    if(sqrtval1.Data.Size > 1) {
                        ImPlot::PlotLine("c1 1/sqrt(N)", &sqrtval1.Data[0].x, &sqrtval1.Data[0].y, sqrtval1.Data.size(), 0,
                                         0, 2 * sizeof(float));
                        ImPlot::PlotLine("c2 1/sqrt(N)", &sqrtval2.Data[0].x, &sqrtval2.Data[0].y, sqrtval2.Data.size(), 0,
                                         0, 2 * sizeof(float));
                    }
                    if(plot_ascbr){
                        ascbr_lower.AddPoint(static_cast<float>(values.Data.front().x),
                                             static_cast<float>(formulas->convergenceValues.front().lower_bound));
                        ascbr_lower.AddPoint(static_cast<float>(values.Data.back().x),
                                             static_cast<float>(formulas->convergenceValues.front().lower_bound));
                        ascbr_upper.AddPoint(static_cast<float>(values.Data.front().x),
                                             static_cast<float>(formulas->convergenceValues.front().upper_bound));
                        ascbr_upper.AddPoint(static_cast<float>(values.Data.back().x),
                                             static_cast<float>(formulas->convergenceValues.front().upper_bound));
                        ImPlot::PlotLine("ASCBR", &ascbr_lower.Data[0].x, &ascbr_lower.Data[0].y, ascbr_lower.Data.size(), 0,
                                         0, 2 * sizeof(float));
                        ImPlot::PlotLine("ASCBR", &ascbr_upper.Data[0].x, &ascbr_upper.Data[0].y, ascbr_upper.Data.size(), 0,
                                         0, 2 * sizeof(float));
                    }

                    //ImPlot::PlotLine("Des/s", tvalues, values_des, IM_ARRAYSIZE(values_des), values_offset);
                    //ImPlot::PlotShaded("Des/s", tvalues, values_des, IM_ARRAYSIZE(values_des), -INFINITY, values_offset);
                    ImPlot::PopStyleVar();

                    ImPlot::EndPlot();
                }

                if(ImGui::Button("Reset Data")){
                    values.Erase();
                }
                static int prune_n = 0;
                ImGui::SliderInt("N", &prune_n, 0, values.Data.size());
                if(ImGui::Button("Prune every N")){
                    int skip_last_n = 0;
                    for (int formulaId = 0; formulaId < formulas->formulas_n.size(); ++formulaId) {
                        formulas->pruneEveryN(4, formulaId, skip_last_n);
                    }
                    for (int i = values.Data.size() - prune_n - skip_last_n; i > 0; i = i - prune_n)
                        values.Data.erase(values.Data.begin() + i);
                }
                ImGui::SameLine();
                if(ImGui::Button("Prune first N")){
                    for (int formulaId = 0; formulaId < formulas->formulas_n.size(); ++formulaId) {
                        formulas->pruneFirstN(100, formulaId);
                    }
                    values.Data.erase(values.Data.begin(), values.Data.begin() + std::min(prune_n, values.Data.size()));
                }

                static int val_ind = 0;
                if(values.Data.Size >= 1) {
                    ImGui::SliderInt("Index", &val_ind, 0, values.Data.size());
                    ImGui::InputFloat2("##", &values.Data[val_ind].x, "%e");
                }
            }
        }
        ImGui::End();
    }

    //
}
