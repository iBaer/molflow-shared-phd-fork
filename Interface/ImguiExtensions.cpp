//
// Created by pbahr on 8/2/21.
//

#ifndef MOLFLOW_PROJ_IMGUIEXTENSIONS_H
#define MOLFLOW_PROJ_IMGUIEXTENSIONS_H

#include "imgui/imgui.h"
#include <string>
#include <imgui/imgui_internal.h>

namespace ImGui {
// Make the UI compact because there are so many fields
    void PushStyleCompact() {
        ImGuiStyle &style = ImGui::GetStyle();
        ImGui::PushStyleVar(
                ImGuiStyleVar_FramePadding,
                ImVec2(style.FramePadding.x, (float) (int) (style.FramePadding.y * 0.60f)));
        ImGui::PushStyleVar(
                ImGuiStyleVar_ItemSpacing,
                ImVec2(style.ItemSpacing.x, (float) (int) (style.ItemSpacing.y * 0.60f)));
    }

    void PopStyleCompact() { ImGui::PopStyleVar(2); }

// Helper to display a little (?) mark which shows a tooltip when hovered.
// In your own code you may want to display an actual icon if you are using a
// merged icon fonts (see docs/FONTS.md)
    void HelpMarker(const char *desc) {
        ImGui::TextDisabled("(?)");
        if (ImGui::IsItemHovered()) {
            ImGui::BeginTooltip();
            ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
            ImGui::TextUnformatted(desc);
            ImGui::PopTextWrapPos();
            ImGui::EndTooltip();
        }
    }

    void PlaceAtWindowCenter(const char *str) {
        const std::string btnText = str;
        float font_size = ImGui::GetFontSize() * btnText.size() / 2;
        ImGui::NewLine();
        ImGui::SameLine((ImGui::GetWindowSize().x / 2) - font_size + (font_size / 2));
    }

    void PlaceAtRegionCenter(const char *str) {
        // float font_size = ImGui::GetFontSize() * btnText.size() / 2;
        float font_size = ImGui::CalcTextSize(str).x;
        ImGui::NewLine();
        ImGui::SameLine((ImGui::GetContentRegionAvail().x * 0.5f) - (font_size / 2));
    }

    void PlaceAtRegionRight(const char *str, bool sameLine) {
        // float font_size = ImGui::GetFontSize() * btnText.size() / 2;
        float font_size = ImGui::CalcTextSize(str).x;
        if (!sameLine)
            ImGui::NewLine();
        ImGui::SameLine(ImGui::GetContentRegionAvail().x -
                        (font_size + ImGui::GetStyle().FramePadding.x * 2));
    }

    bool InputRightSide(const char *desc, double *val, const char *format) {
        double tmp = *val;
        ImGui::AlignTextToFramePadding();
        ImGui::Text("%s:", desc);
        {
            // Move to right side
            ImGui::SameLine((ImGui::GetContentRegionAvail().x) - 100.0f);
            ImGui::PushItemWidth(100.0f);
            ImGui::PushID(desc);
            ImGui::InputDouble("", val, 0.00f, 0.0f, format);
            ImGui::PopID();
            ImGui::PopItemWidth();
        }

        return *val != tmp; // true if changed
    }

// Add spacing of checkbox width
    void AddCheckboxWidthSpacing() {
        ImGui::NewLine();
        ImGui::SameLine(ImGui::GetFrameHeightWithSpacing() +
                        ImGui::GetStyle().FramePadding.y * 4);
    }

    bool
    BufferingBar(const char *label, float value, const ImVec2 &size_arg, const ImU32 &bg_col, const ImU32 &fg_col) {
        ImGuiWindow *window = GetCurrentWindow();
        if (window->SkipItems)
            return false;

        ImGuiContext &g = *GImGui;
        const ImGuiStyle &style = g.Style;
        const ImGuiID id = window->GetID(label);

        ImVec2 pos = window->DC.CursorPos;
        ImVec2 size = size_arg;
        size.x -= style.FramePadding.x * 2;

        const ImRect bb(pos, ImVec2(pos.x + size.x, pos.y + size.y));
        ItemSize(bb, style.FramePadding.y);
        if (!ItemAdd(bb, id))
            return false;

        // Render
        const float circleStart = size.x * 0.7f;
        const float circleEnd = size.x;
        const float circleWidth = circleEnd - circleStart;

        window->DrawList->AddRectFilled(bb.Min, ImVec2(pos.x + circleStart, bb.Max.y), bg_col);
        window->DrawList->AddRectFilled(bb.Min, ImVec2(pos.x + circleStart * value, bb.Max.y), fg_col);

        const float t = g.Time;
        const float r = size.y / 2;
        const float speed = 1.5f;

        const float a = speed * 0;
        const float b = speed * 0.333f;
        const float c = speed * 0.666f;

        const float o1 = (circleWidth + r) * (t + a - speed * (int) ((t + a) / speed)) / speed;
        const float o2 = (circleWidth + r) * (t + b - speed * (int) ((t + b) / speed)) / speed;
        const float o3 = (circleWidth + r) * (t + c - speed * (int) ((t + c) / speed)) / speed;

        window->DrawList->AddCircleFilled(ImVec2(pos.x + circleEnd - o1, bb.Min.y + r), r, bg_col);
        window->DrawList->AddCircleFilled(ImVec2(pos.x + circleEnd - o2, bb.Min.y + r), r, bg_col);
        window->DrawList->AddCircleFilled(ImVec2(pos.x + circleEnd - o3, bb.Min.y + r), r, bg_col);

        return true;
    }

    bool Spinner(const char *label, float radius, int thickness, const ImU32 &color) {
        ImGuiWindow *window = GetCurrentWindow();
        if (window->SkipItems)
            return false;

        ImGuiContext &g = *GImGui;
        const ImGuiStyle &style = g.Style;
        const ImGuiID id = window->GetID(label);

        ImVec2 pos = window->DC.CursorPos;
        ImVec2 size((radius) * 2, (radius + style.FramePadding.y) * 2);

        const ImRect bb(pos, ImVec2(pos.x + size.x, pos.y + size.y));
        ItemSize(bb, style.FramePadding.y);
        if (!ItemAdd(bb, id))
            return false;

        // Render
        window->DrawList->PathClear();

        int num_segments = 30;
        int start = abs(ImSin(g.Time * 1.8f) * (num_segments - 5));

        const float a_min = IM_PI * 2.0f * ((float) start) / (float) num_segments;
        const float a_max = IM_PI * 2.0f * ((float) num_segments - 3) / (float) num_segments;

        const ImVec2 centre = ImVec2(pos.x + radius, pos.y + radius + style.FramePadding.y);

        for (int i = 0; i < num_segments; i++) {
            const float a = a_min + ((float) i / (float) num_segments) * (a_max - a_min);
            window->DrawList->PathLineTo(ImVec2(centre.x + ImCos(a + g.Time * 8) * radius,
                                                centre.y + ImSin(a + g.Time * 8) * radius));
        }

        window->DrawList->PathStroke(color, false, thickness);
        return true;
    }

    bool Loader(float& progress, float& time){
        static float /*progress = 0.0f,*/ progress_dir = 1.0f;
        time += ImGui::GetIO().DeltaTime;
        progress += progress_dir * (1.0f / 60.0f) * ImGui::GetIO().DeltaTime;
        if (progress >= +1.0f) {
            progress = +1.0f;
            progress_dir *= -1.0f;
        }
        if (progress <= -0.0f) {
            progress = -0.0f;
            progress_dir *= -1.0f;
        }
        //}
        const ImU32 col = ImGui::GetColorU32(ImGuiCol_ButtonHovered);
        const ImU32 bg = ImGui::GetColorU32(ImGuiCol_Button);
        const float size = 36.0f;
        const int thickness = 6.0f;
        const float radius = 0.5f * (size - thickness);

        const ImGuiViewport *viewport = ImGui::GetMainViewport();
        ImGui::SetNextWindowPos(true ?
                                ImVec2(viewport->Size.x - viewport->WorkSize.x * 0.5f, 0.34f * viewport->Size.y)
                                     : viewport->Pos);
        ImGui::SetNextWindowSize(ImVec2(
                164.0f+size+3.0f*ImGui::GetStyle().ItemInnerSpacing.x,
                size+2.0f*ImGui::GetStyle().WindowPadding.y)
                );
        /*use_work_area ? ImVec2(viewport->WorkSize.x * 0.25f, viewport->WorkSize.y) : viewport->Size);*/
        static ImGuiWindowFlags flags =
                ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize |
                ImGuiWindowFlags_NoSavedSettings;
        bool open = true;

        if (ImGui::Begin("Loader", &open, flags)) {
            ImGui::Spinner("##spinner", radius, thickness, col);
            ImGui::SameLine(0.0f, ImGui::GetStyle().ItemInnerSpacing.x);
            char buf[32];
            if(time >= 60)
                sprintf(buf, "%d:%2dm", (int)(time / 60.0), ((int)time % 60));
            else
                sprintf(buf, "%ds", (int)(time));
            ImGui::ProgressBar(progress, ImVec2(164.0f, size), buf);
            //ImGui::SameLine();
            ImGui::End();
        }

    }

}
#endif //MOLFLOW_PROJ_IMGUIEXTENSIONS_H
