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

#include "ImguiSidebar.h"
#include "ImguiExtensions.h"

#if defined(MOLFLOW)
#include "../../src/MolFlow.h"
#include "../../src/MolflowGeometry.h"
#else
#include "../../src/SynRad.h"
#include "../../src/SynradGeometry.h"
#endif

#include "Facet_shared.h"
#include "../../src/Interface/Viewer3DSettings.h"

#include "imgui/imgui.h"
#include <imgui/imgui_internal.h>
#include <Helper/FormatHelper.h>


// Dummy data structure that we use for the Table demo.
// (pre-C++11 doesn't allow us to instantiate ImVector<MyItem> template if this structure if defined inside the demo function)
namespace {
// We are passing our own identifier to TableSetupColumn() to facilitate identifying columns in the sorting code.
// This identifier will be passed down into ImGuiTableSortSpec::ColumnUserID.
// But it is possible to omit the user id parameter of TableSetupColumn() and just use the column index instead! (ImGuiTableSortSpec::ColumnIndex)
// If you don't use sorting, you will generally never care about giving column an ID!
    enum FacetDataColumnID {
        FacetDataColumnID_ID,
        FacetDataColumnID_Hits,
        FacetDataColumnID_Des,
        FacetDataColumnID_Abs
    };

    struct FacetData {
        int ID;
        size_t hits;
        size_t des;
        double abs;

        // We have a problem which is affecting _only this demo_ and should not affect your code:
        // As we don't rely on std:: or other third-party library to compile dear imgui, we only have reliable access to qsort(),
        // however qsort doesn't allow passing user data to comparing function.
        // As a workaround, we are storing the sort specs in a static/global for the comparing function to access.
        // In your own use case you would probably pass the sort specs to your sorting/comparing functions directly and not use a global.
        // We could technically call ImGui::TableGetSortSpecs() in CompareWithSortSpecs(), but considering that this function is called
        // very often by the sorting algorithm it would be a little wasteful.
        static const ImGuiTableSortSpecs *s_current_sort_specs;

        // Compare function to be used by qsort()
        static int IMGUI_CDECL CompareWithSortSpecs(const void *lhs, const void *rhs) {
            const FacetData *a = (const FacetData *) lhs;
            const FacetData *b = (const FacetData *) rhs;
            for (int n = 0; n < s_current_sort_specs->SpecsCount; n++) {
                // Here we identify columns using the ColumnUserID value that we ourselves passed to TableSetupColumn()
                // We could also choose to identify columns based on their index (sort_spec->ColumnIndex), which is simpler!
                const ImGuiTableColumnSortSpecs *sort_spec = &s_current_sort_specs->Specs[n];
                int delta = 0;
                switch (sort_spec->ColumnUserID) {
                    case FacetDataColumnID_ID:
                        delta = (a->ID - b->ID);
                        break;
                    case FacetDataColumnID_Hits:
                        delta = (a->hits > b->hits) ? 1 : (a->hits == b->hits) ? 0 : -1;
                        break;
                    case FacetDataColumnID_Des:
                        delta = (a->des > b->des) ? 1 : (a->des == b->des) ? 0 : -1;
                        break;
                    case FacetDataColumnID_Abs:
                        delta = (a->abs > b->abs) ? 1 : (a->abs == b->abs) ? 0 : -1;
                        break;
                    default:
                        IM_ASSERT(0);
                        break;
                }
                if (delta > 0)
                    return (sort_spec->SortDirection == ImGuiSortDirection_Ascending) ? +1 : -1;
                if (delta < 0)
                    return (sort_spec->SortDirection == ImGuiSortDirection_Ascending) ? -1 : +1;
            }

            // qsort() is instable so always return a way to differenciate items.
            // Your own compare function may want to avoid fallback on implicit sort specs e.g. a Name compare if it wasn't already part of the sort specs.
            return (a->ID - b->ID);
        }
    };

    const ImGuiTableSortSpecs *FacetData::s_current_sort_specs = nullptr;
}

// Sidebar containing 3d viewer settings, facet settings and simulation data
#if defined(MOLFLOW)
void ShowAppSidebar(bool *p_open, MolFlow *mApp, Geometry *geom, bool *show_global, bool *newViewer) {
#else
void ShowAppSidebar(bool *p_open, SynRad *mApp, Geometry *geom, bool *show_global, bool *newViewer) {
#endif
    const float PAD = 10.0f;
    static int corner = 0;
    static float TEXT_BASE_HEIGHT = ImGui::GetTextLineHeightWithSpacing();
    static float TEXT_BASE_WIDTH = ImGui::CalcTextSize("A").x;

    ImGuiIO &io = ImGui::GetIO();

    const ImGuiViewport *viewport = ImGui::GetMainViewport();
    static bool use_work_area = false; // experiment with work area vs normal area
    ImGui::SetNextWindowPos(use_work_area ?
                            ImVec2(viewport->Size.x - viewport->WorkSize.x * 0.25f, viewport->WorkPos.y)
                                          : viewport->Pos);
    // Set size initially (to be able for resize) to take the full height on the right side
    ImGui::SetNextWindowSize(
            use_work_area ? ImVec2(viewport->WorkSize.x * 0.25f, viewport->WorkSize.y) : ImVec2(0, viewport->Size.y), ImGuiCond_FirstUseEver);

    static ImGuiWindowFlags flags =
            ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;
    /*ImGuiWindowFlags_NoDecoration | *//*ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove *//*| ImGuiWindowFlags_NoResize *//*|
            ImGuiWindowFlags_NoSavedSettings;*/

    if (!use_work_area) // example from ShowExampleAppSimpleOverlay() [demo] for corner positioning
    {
        int corner = 1; // top right
        const float PAD = 0.0f;//const float PAD = 10.0f;
        const ImGuiViewport* viewport = ImGui::GetMainViewport();
        ImVec2 work_pos = viewport->WorkPos; // Use work area to avoid menu-bar/task-bar, if any!
        ImVec2 work_size = viewport->WorkSize;
        ImVec2 window_pos, window_pos_pivot;
        window_pos.x = (corner & 1) ? (work_pos.x + work_size.x - PAD) : (work_pos.x + PAD);
        window_pos.y = (corner & 2) ? (work_pos.y + work_size.y - PAD) : (work_pos.y + PAD);
        window_pos_pivot.x = (corner & 1) ? 1.0f : 0.0f;
        window_pos_pivot.y = (corner & 2) ? 1.0f : 0.0f;
        ImGui::SetNextWindowPos(window_pos, ImGuiCond_Always, window_pos_pivot);
        //ImGui::SetNextWindowViewport(viewport->ID);
        //ImGui::SetNextWindowViewport(viewport->ID);
        flags |= ImGuiWindowFlags_NoMove;
    }

    if (ImGui::Begin("[BETA] Molflow Sidebar", p_open, flags)) {
#if defined(DEBUG)
        if (ImGui::CollapsingHeader("[DEMO] Window flags")) {

            ImGui::Checkbox("Use work area instead of main area", &use_work_area);
            ImGui::SameLine();
            ImGui::HelpMarker(
                    "Main Area = entire viewport,\nWork Area = entire viewport minus sections used by the main menu bars, task bars etc.\n\nEnable the main-menu bar in Examples menu to see the difference.");

            ImGui::CheckboxFlags("ImGuiWindowFlags_NoBackground", &flags, ImGuiWindowFlags_NoBackground);
            ImGui::CheckboxFlags("ImGuiWindowFlags_NoDecoration", &flags, ImGuiWindowFlags_NoDecoration);
            ImGui::Indent();
            ImGui::CheckboxFlags("ImGuiWindowFlags_NoTitleBar", &flags, ImGuiWindowFlags_NoTitleBar);
            ImGui::CheckboxFlags("ImGuiWindowFlags_NoCollapse", &flags, ImGuiWindowFlags_NoCollapse);
            ImGui::CheckboxFlags("ImGuiWindowFlags_NoScrollbar", &flags, ImGuiWindowFlags_NoScrollbar);
            ImGui::Unindent();

            if (p_open && ImGui::Button("Close this window"))
                *p_open = false;
        }
#endif
        if (ImGui::CollapsingHeader("3D Viewer settings", ImGuiTreeNodeFlags_DefaultOpen)) {
            auto curViewer = mApp->curViewer;
            auto viewer = mApp->viewer[curViewer];
            if (ImGui::BeginTable("table_3dviewer", 3, ImGuiTableFlags_None)) {
                /*ImGui::TableSetupColumn("col1");
                ImGui::TableSetupColumn("col2");
                ImGui::TableSetupColumn("col3");
                ImGui::TableHeadersRow();*/
                ImGui::TableNextRow();
                {
                    ImGui::TableNextColumn();
                    ImGui::Checkbox("Rules", &viewer->showRule);
                    ImGui::TableNextColumn();
                    ImGui::Checkbox("Normals", &viewer->showNormal);
                    ImGui::TableNextColumn();
                    ImGui::Checkbox(fmt::format("{}", u8"u\u20d7,v\u20d7").c_str(), &viewer->showUV);
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    ImGui::Checkbox("Lines", &viewer->showLine);
                    ImGui::TableNextColumn();
                    ImGui::Checkbox("Leaks", &viewer->showLeak);
                    ImGui::TableNextColumn();
                    ImGui::Checkbox("Hits", &viewer->showHit);
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    ImGui::Checkbox("Volume", &viewer->showVolume);
                    ImGui::TableNextColumn();
                    ImGui::Checkbox("Texture", &viewer->showTexture);
                    ImGui::TableNextColumn();
                    ImGui::Checkbox("FacetIDs", &viewer->showFacetId);
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    if(ImGui::Button("<< View")){
                        *newViewer = true;
                    }
                    ImGui::TableNextColumn();
                    ImGui::Checkbox("Indices", &viewer->showIndex);
                    ImGui::TableNextColumn();
                    ImGui::Checkbox("VertexIDs", &viewer->showVertexId);
                }

                ImGui::EndTable();

            }
        }

        std::string title;
        if(geom->GetNbSelectedFacets() > 1) {
            title = fmt::format("Selected Facet ({} selected)", geom->GetNbSelectedFacets());
        }
        else if(geom->GetNbSelectedFacets() == 1){
            title = fmt::format("Selected Facet (#{})", geom->GetSelectedFacets().front());
        }
        else {
            title = fmt::format("Selected Facet (none)");
        }
        if (ImGui::CollapsingHeader(title.c_str(), ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x*0.35f);
            ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x*0.25f);
            size_t selected_facet_id = geom->GetNbSelectedFacets() ? geom->GetSelectedFacets().front() : 0;
            auto sel = geom->GetNbSelectedFacets() ? geom->GetFacet(selected_facet_id) : nullptr;
#if defined(MOLFLOW)
            if (ImGui::TreeNodeEx("Particles in", ImGuiTreeNodeFlags_DefaultOpen)) {
                static int des_idx = 0;
                if(sel) des_idx = sel->sh.desorbType;
                if (ImGui::Combo("Desorption", &des_idx, "None\0Uniform\0Cosine\0Cosine\u207f\0Recorded\0")) {
                    switch (des_idx) {
                        case 0:
                            fmt::print("none");
                            break;
                        case 1:
                            fmt::print("uni");
                            break;
                        case 2:
                            fmt::print("cos");
                            break;
                        case 3:
                            fmt::print("cos \u207f");
                            break;
                        case 4:
                            fmt::print("rec");
                            break;
                    }
                }

                static bool use_og_area = false;
                bool use_og = !use_og_area;
                static double og = 1.0;
                static double og_area = 1.0;

                if(sel) og = sel->sh.outgassing;
                if (ImGui::Checkbox(u8"Outgassing [mbar\u00b7l/s]", &use_og)) {
                    use_og_area = !use_og;
                }
                ImGui::SameLine();
                ImGui::InputDouble("##in", &og);
                ImGui::Checkbox(u8"Outg/area [mbar\u00b7l/s/cm\u00b2]", &use_og_area);
                ImGui::SameLine();
                ImGui::InputDouble("##ina", &og);
                ImGui::TreePop();
            }
#endif

            if (ImGui::TreeNodeEx("Particles out", ImGuiTreeNodeFlags_DefaultOpen)) {
                static double sf = 1.0;
                static double ps = 1.0;
                if(sel) sf = sel->sh.sticking;
                ImGui::InputRightSide("Sticking factor", &sf);
                ImGui::InputRightSide("Pumping speed [l/s]", &ps);
                ImGui::TreePop();
            }
            ImGui::PopItemWidth();
            {
                static int sides_idx = 0;
                if(sel) sides_idx = sel->sh.is2sided;
                if (ImGui::Combo("Sides", &sides_idx, "1 Sided\0 2 Sided\0")) {
                    switch (sides_idx) {
                        case 0:
                            fmt::print("1s");
                            break;
                        case 1:
                            fmt::print("2s");
                            break;
                    }
                }

                static double opacity = 1.0;
                if(sel) opacity = sel->sh.opacity;
                ImGui::InputRightSide("Opacity", &opacity);

#if defined(MOLFLOW)
                static double temp = 1.0;
                if(sel) temp = sel->sh.temperature;
                ImGui::InputRightSide(u8"Temperature [\u00b0\u212a]", &temp);
#endif

                static double area = 1.0;
                if(sel) area = sel->sh.area;
                ImGui::InputRightSide(u8"Area [cm\u00b2]", &area);

                static int prof_idx = 0;
                if(sel) prof_idx = sel->sh.profileType;
                if (ImGui::Combo("Profile", &prof_idx,
                                 "None\0Pressure u\0Pressure v\0Incident angle\0Speed distribution\0Orthogonal velocity\0 Tangential velocity\0")) {
                    switch (prof_idx) {
                        case 0:
                            fmt::print("None");
                            break;
                        case 1:
                            fmt::print("Pu");
                            break;
                        case 2:
                            fmt::print("Pv");
                            break;
                        case 3:
                            fmt::print("Angle");
                            break;
                        case 4:
                            fmt::print("Speed");
                            break;
                        case 5:
                            fmt::print("Ortho vel");
                            break;
                        case 6:
                            fmt::print("Tangen vel");
                            break;
                    }
                }
            }
            ImGui::PopItemWidth();
        }


        if (ImGui::CollapsingHeader("Simulation", ImGuiTreeNodeFlags_DefaultOpen)) {
            if(ImGui::Button("<< Sim")){
                *show_global = !*show_global;
            }
            ImGui::SameLine();
            {
                if (mApp->worker.IsRunning()) {
                    title = fmt::format("Pause");
                } else if (mApp->worker.globalHitCache.globalHits.nbMCHit > 0) {
                    title = fmt::format("Resume");
                } else {
                    title = fmt::format("Begin");
                }
            }
            if(ImGui::Button(title.c_str())){
                mApp->changedSinceSave = true;
                mApp->StartStopSimulation();
            }
            ImGui::SameLine();
            if(!mApp->worker.IsRunning() && mApp->worker.globalHitCache.globalHits.nbDesorbed > 0)
                ImGui::BeginDisabled();
            if(ImGui::Button("Reset")){
                mApp->changedSinceSave = true;
                mApp->ResetSimulation();
            }
            if(!mApp->worker.IsRunning() && mApp->worker.globalHitCache.globalHits.nbDesorbed > 0)
                ImGui::EndDisabled();

            ImGui::Checkbox("Auto update scene", &mApp->autoFrameMove);
            if(mApp->autoFrameMove){
                ImGui::BeginDisabled();
            }
            if(ImGui::Button("Update")){
                mApp->updateRequested = true;
                mApp->FrameMove();
            }
            if(mApp->autoFrameMove){
                ImGui::EndDisabled();
            }

            ImVec2 outer_size = ImVec2(std::max(0.0f, ImGui::GetContentRegionAvail().x), 0.0f);
            // for layouting Simulation stats are inside a table/grid
            if (ImGui::BeginTable("simugrid", 2, ImGuiTableFlags_None/*ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingFixedFit |
                                                 ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersOuter |
                                                 ImGuiTableFlags_BordersV | ImGuiTableFlags_Resizable |
                                                 ImGuiTableFlags_Reorderable | ImGuiTableFlags_Hideable |
                                                 ImGuiTableFlags_Sortable*/, outer_size)) {

                static std::string hit_stat;
                static std::string des_stat;
                bool runningState = mApp->worker.IsRunning();
                if ((mApp->worker.simuTimer.Elapsed() <= 2.0f) && runningState) {
                    hit_stat = "Starting...";
                    des_stat = "Starting...";
                } else {
                    double current_avg = 0.0;
                    if (!runningState) current_avg = mApp->hps_runtotal.avg();
                    else current_avg = (current_avg != 0.0) ? current_avg : mApp->hps.last();

                    hit_stat = fmt::format("{} ({})",
                                           Util::formatInt(mApp->worker.globalHitCache.globalHits.nbMCHit, "hit"),
                                           Util::formatPs(current_avg, "hit"));

                    current_avg = 0.0;
                    if (!runningState) current_avg = mApp->dps_runtotal.avg();
                    else current_avg = (current_avg != 0.0) ? current_avg : mApp->dps.last();

                    des_stat = fmt::format("{} ({})",
                                           Util::formatInt(mApp->worker.globalHitCache.globalHits.nbDesorbed, "des"),
                                           Util::formatPs(current_avg, "des"));
                }


                ImGui::TableSetupColumn("name", ImGuiTableColumnFlags_WidthFixed, 0.0f);
                ImGui::TableSetupColumn("field", ImGuiTableColumnFlags_WidthStretch, 0.0f);
                //ImGui::TableHeadersRow();

                static char inputName[128] = "";
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::Text("Hits");
                ImGui::TableNextColumn();
                strcpy(inputName, fmt::format("{}", hit_stat).c_str());
                ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::InputText("##hit", inputName, IM_ARRAYSIZE(inputName));
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::Text("Des.");
                ImGui::TableNextColumn();
                strcpy(inputName, fmt::format("{}", des_stat).c_str());
                ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::InputText("##des", inputName, IM_ARRAYSIZE(inputName));
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::Text("Leaks");
                ImGui::TableNextColumn();
                strcpy(inputName, fmt::format("{}", "None").c_str());
                ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::InputText("##leak", inputName, IM_ARRAYSIZE(inputName));
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::Text("Time");
                if (runningState)
                    strcpy(inputName,
                           fmt::format("Running: {}", Util::formatTime(mApp->worker.simuTimer.Elapsed())).c_str());
                else
                    strcpy(inputName,
                           fmt::format("Stopped: {}", Util::formatTime(mApp->worker.simuTimer.Elapsed())).c_str());
                ImGui::TableNextColumn();
                ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::InputText("##time", inputName, IM_ARRAYSIZE(inputName));

                ImGui::EndTable();
            }
        }

        try {
            // Facet list
            if (geom->IsLoaded()) {
                // TODO: Colors for moment highlighting

                // Create item list that is sortable etc.
                static ImVector<FacetData> items;
                if (items.Size != geom->GetNbFacet()) {
                    items.resize(geom->GetNbFacet(), FacetData());
                    for (int n = 0; n < items.Size; n++) {
                        InterfaceFacet *f = geom->GetFacet(n);
                        FacetData &item = items[n];
                        item.ID = n;
                        item.hits = f->facetHitCache.nbMCHit;
                        item.des = f->facetHitCache.nbDesorbed;
                        item.abs = f->facetHitCache.nbAbsEquiv;
                    }
                }
                else if (mApp->worker.IsRunning()){
                    for (int n = 0; n < items.Size; n++) {
                        InterfaceFacet *f = geom->GetFacet(n);
                        FacetData &item = items[n];
                        item.hits = f->facetHitCache.nbMCHit;
                        item.des = f->facetHitCache.nbDesorbed;
                        item.abs = f->facetHitCache.nbAbsEquiv;
                    }
                }

                static ImGuiTableFlags tFlags =
                        ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingFixedFit |
                        /*ImGuiTableFlags_RowBg | */ImGuiTableFlags_BordersOuter |
                        ImGuiTableFlags_BordersV | ImGuiTableFlags_Resizable |
                        ImGuiTableFlags_Reorderable | ImGuiTableFlags_Hideable |
                        ImGuiTableFlags_Sortable;

                ImVec2 outer_size = ImVec2(0.0f, std::max(ImGui::GetContentRegionAvail().y, TEXT_BASE_HEIGHT * 8.f));
                if (ImGui::BeginTable("facetlist", 4, tFlags, outer_size)) {
                    ImGui::TableSetupScrollFreeze(0, 1); // Make top row always visible
                    ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, 0.0f, FacetDataColumnID_ID);
                    ImGui::TableSetupColumn("Hits", ImGuiTableColumnFlags_WidthFixed, 0.0f, FacetDataColumnID_Hits);
                    ImGui::TableSetupColumn("Des", ImGuiTableColumnFlags_WidthStretch, 0.0f, FacetDataColumnID_Des);
                    ImGui::TableSetupColumn("Abs", ImGuiTableColumnFlags_WidthStretch, 0.0f, FacetDataColumnID_Abs);
                    ImGui::TableHeadersRow();

                    // Sort our data if sort specs have been changed!
                    if (ImGuiTableSortSpecs *sorts_specs = ImGui::TableGetSortSpecs())
                        if (sorts_specs->SpecsDirty) {
                            FacetData::s_current_sort_specs = sorts_specs; // Store in variable accessible by the sort function.
                            if (items.Size > 1)
                                qsort(&items[0], (size_t) items.Size, sizeof(items[0]),
                                      FacetData::CompareWithSortSpecs);
                            FacetData::s_current_sort_specs = nullptr;
                            sorts_specs->SpecsDirty = false;
                        }

                    // Demonstrate using clipper for large vertical lists
                    ImGuiListClipper clipper;
                    clipper.Begin(items.size());
                    while (clipper.Step()) {
                        for (int i = clipper.DisplayStart; i < clipper.DisplayEnd; i++) {
                            FacetData *item = &items[i];
                            //ImGui::PushID(item->ID);
                            ImGui::TableNextRow();
                            ImGui::TableSetColumnIndex(0);
                            ImGui::Text("%d", item->ID);
                            ImGui::TableNextColumn();
                            ImGui::Text("%zd", item->hits);
                            ImGui::TableNextColumn();
                            ImGui::Text("%zd", item->des);
                            ImGui::TableNextColumn();
                            ImGui::Text("%g", item->abs);
                            //ImGui::PopID();
                        }
                    }
                    ImGui::EndTable();
                }

            }
        }
        catch (const std::exception &e) {
            char errMsg[512];
            sprintf(errMsg, "%s\nError while updating facet hits", e.what());
            /*GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);*/
        }
    }
    ImGui::End();
}
