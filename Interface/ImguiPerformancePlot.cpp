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

#include "ImguiPerformancePlot.h"
#include "imgui/imgui.h"
#include <implot/implot.h>

#include "Interface.h"

// Demonstrate creating a simple static window with no decoration
// + a context-menu to choose which corner of the screen to use.
void ShowPerfoPlot(bool *p_open, Interface *mApp) {
    ImGuiIO &io = ImGui::GetIO();

    // Always center this window when appearing
    ImVec2 center = ImGui::GetMainViewport()->GetCenter();
    ImGui::SetNextWindowPos(center, ImGuiCond_Appearing,
                            ImVec2(0.5f, 0.5f));
    ImGui::SetNextWindowSize(ImVec2(400, 400), ImGuiCond_FirstUseEver);

    static ImGuiWindowFlags flags =/*
            ImGuiWindowFlags_AlwaysAutoResize |*/
            ImGuiWindowFlags_NoSavedSettings;

    if (ImGui::Begin("Performance Plot", p_open, flags)) {

        // Fill an array of contiguous float values to plot
        // Tip: If your float aren't contiguous but part of a structure, you can pass a pointer to your first float
        // and the sizeof() of your structure in the "stride" parameter.
        static float values[20] = {0.0f};
        static float values_des[20] = {0.0f};
        static float tvalues[20] = {0.0f};
        static int values_offset = 0;
        static auto refresh_time = ImGui::GetTime();
        if (!true || refresh_time == 0.0) // force
            refresh_time = ImGui::GetTime();
        auto now_time = ImGui::GetTime();
        if (mApp->worker.IsRunning() && difftime(now_time, refresh_time) > 1.0 &&
            mApp->hps.eventsAtTime.size() >= 2)
        {
            //static float phase = 0.0f;
            values[values_offset] = static_cast<float>(mApp->hps.avg());
            values_des[values_offset] = static_cast<float>(mApp->dps.avg());
            tvalues[values_offset] = static_cast<float>(now_time);
            if (values[values_offset] != values[(values_offset - 1) % IM_ARRAYSIZE(values)])
                values_offset = (values_offset + 1) % IM_ARRAYSIZE(values);
            //phase += 0.10f * values_offset;
            refresh_time = now_time;
        }

        // Plots can display overlay texts
        // (in this example, we will display an average value)
        {
            float average = 0.0f;
            for (float value: values)
                average += value;
            average /= (float) IM_ARRAYSIZE(values);

            float max_val = values[0];
            float min_val = values[0];
            for (int i = 1; i < IM_ARRAYSIZE(values); ++i) {
                if (values[i] > max_val) {
                    max_val = values[i];
                }
                if (values[i] < min_val) {
                    min_val = values[i];
                }
                if (values_des[i] > max_val) {
                    max_val = values_des[i];
                }
                if (values_des[i] < min_val) {
                    min_val = values_des[i];
                }
            }
            char overlay[32];
            sprintf(overlay, "avg %f hit/s", average);
            //ImGui::PlotLines(""*//*"Hit/s"*//*, values, IM_ARRAYSIZE(values), values_offset, overlay, min_val * 0.95f, max_val * 1.05f,ImVec2(0, 80.0f));

            ImPlot::SetNextAxisLimits(std::max(0.0f, min_val * 0.8f),max_val * 1.2f, ImGuiCond_Always);
            if (ImPlot::BeginPlot("##Perfo", "time (s)", "performance (hit/s)", ImVec2(-1, -1), ImPlotFlags_None,
                                  ImPlotAxisFlags_AutoFit/* | ImPlotAxisFlags_Time*//*, ImPlotAxisFlags_AutoFit*/)) {
                ImPlot::PushStyleVar(ImPlotStyleVar_FillAlpha, 0.25f);
                ImPlot::PlotLine("Hit/s", tvalues, values, IM_ARRAYSIZE(values), values_offset);
                ImPlot::PlotShaded("Hit/s", tvalues, values, IM_ARRAYSIZE(values), -INFINITY, values_offset);
                ImPlot::PlotLine("Des/s", tvalues, values_des, IM_ARRAYSIZE(values_des), values_offset);
                ImPlot::PlotShaded("Des/s", tvalues, values_des, IM_ARRAYSIZE(values_des), -INFINITY, values_offset);
                ImPlot::PopStyleVar();
                ImPlot::EndPlot();
            }
        }
    }
    ImGui::End();
}
