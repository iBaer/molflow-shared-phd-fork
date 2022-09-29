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
#include "ConvergencePlotter.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLCombo.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLToggle.h"
#include "Helper/MathTools.h"
#include "GLApp/GLList.h"
#include "GLApp/GLChart/GLChart.h"
#include "GLApp/GLParser.h"
#include "Formulas.h"
#include "Geometry_shared.h"
#include "Facet_shared.h"

#include <cmath>
#include <Helper/StringHelper.h>

#if defined(MOLFLOW)

#include "../../src/MolFlow.h"

#endif

#if defined(SYNRAD)
#include "../src/SynRad.h"
#endif

constexpr size_t max_vector_size() { return 16384; };

/**
* \brief Constructor with initialisation for Convergence plotter window (Tools/Convergence Plotter)
*/
ConvergencePlotter::ConvergencePlotter(Worker *appWorker, std::shared_ptr<Formulas> formulas)
        : GLWindow(), views{} {

    int wD = 625;
    int hD = 400;

    SetTitle("Convergence plotter");
    SetIconfiable(true);
    nbView = 0;
    worker = appWorker;
    formula_ptr = formulas;

    lastUpdate = 0.0f;

    chart = new GLChart(0);
    chart->SetBorder(BORDER_BEVEL_IN);
    chart->GetY1Axis()->SetGridVisible(true);
    chart->GetXAxis()->SetGridVisible(true);
    chart->GetY1Axis()->SetAutoScale(true);
    chart->GetY2Axis()->SetAutoScale(true);
    chart->GetY1Axis()->SetAnnotation(VALUE_ANNO);
    chart->GetXAxis()->SetAnnotation(VALUE_ANNO);
    chart->GetXAxis()->SetName("Number of desorptions");

    Add(chart);

    dismissButton = new GLButton(0, "Dismiss");
    Add(dismissButton);

    pruneEveryNButton = new GLButton(0, "Remove every 4th");
    Add(pruneEveryNButton);
    pruneFirstNButton = new GLButton(0, "Remove first 100");
    Add(pruneFirstNButton);

    addButton = new GLButton(0, "Add curve");
    Add(addButton);

    removeButton = new GLButton(0, "Remove curve");
    Add(removeButton);

    removeAllButton = new GLButton(0, "Remove all");
    Add(removeAllButton);

    profCombo = new GLCombo(0);
    profCombo->SetEditable(true);
    Add(profCombo);

    logYToggle = new GLToggle(0, "Log Y");
    Add(logYToggle);

    colorToggle = new GLToggle(0, "Colorblind mode");
    Add(colorToggle);

    fixedLineWidthText = new GLLabel("Change linewidth:");
    Add(fixedLineWidthText);
    fixedLineWidthButton = new GLButton(0, "-> Apply linewidth");
    Add(fixedLineWidthButton);
    fixedLineWidthField = new GLTextField(0, "2");
    fixedLineWidthField->SetEditable(true);
    Add(fixedLineWidthField);

    formulaText = new GLTextField(0, "");
    formulaText->SetEditable(true);
    Add(formulaText);

    formulaBtn = new GLButton(0, "-> Plot expression");
    Add(formulaBtn);

    panel2 = new GLTitledPanel("Convergence Band Rule");
    panel2->SetClosable(false);
    Add(panel2);
    convEpsText = new GLLabel("Convergence precision {10^-d}: d=");
    panel2->Add(convEpsText);
    convEpsField = new GLTextField(0, "5");
    convEpsField->SetEditable(true);
    panel2->Add(convEpsField);
    convBandLenText = new GLLabel("Convergence band length: C=");
    panel2->Add(convBandLenText);
    convBandLenField = new GLTextField(0, "51");
    convBandLenField->SetEditable(true);
    panel2->Add(convBandLenField);

    convApplyButton = new GLButton(0, "-> Apply changes");
    panel2->Add(convApplyButton);

    shapeParamField = new GLTextField(0, "0");
    panel2->Add(shapeParamField);
    shapeParamField2 = new GLTextField(0, "0");
    panel2->Add(shapeParamField2);
    absOrRelToggle = new GLToggle(0, "Absolute precision");
    absOrRelToggle->SetState(formula_ptr->useAbsEps);
    panel2->Add(absOrRelToggle);

    // Center dialog
    int wS, hS;
    GLToolkit::GetScreenSize(&wS, &hS);
    int xD = (wS - wD) / 2;
    int yD = (hS - hD) / 2;
    SetBounds(xD, yD, wD, hD);
    SetResizable(true);
    SetMinimumSize(wD, 220);

    RestoreDeviceObjects();

}

// copy constructor
ConvergencePlotter::ConvergencePlotter(const ConvergencePlotter& copy) : ConvergencePlotter(copy.worker, copy.formula_ptr) {
    //nbView = copy.nbView;
    worker = copy.worker;
    formula_ptr = copy.formula_ptr;

    lastUpdate = 0.0f;

    profCombo->SetSelectedIndex(copy.profCombo->GetSelectedIndex());
    logYToggle->SetState(copy.logYToggle->GetState());
    colorToggle->SetState(copy.colorToggle->GetState());

    if (!colorToggle->GetState())
        chart->SetColorSchemeDefault();
    else
        chart->SetColorSchemeColorblind();

    SetViews(copy.GetViews());
}


/**
* \brief Sets positions and sizes of all UI elements
* \param x x-coordinate of the element
* \param y y-coordinate of the element
* \param w width of the element
* \param h height of the element
*/
void ConvergencePlotter::SetBounds(int x, int y, int w, int h) {

    size_t lineHeightDiff = 45 + 25;
    panel2->SetBounds(5, h - lineHeightDiff - 25, w - 10, 25*2+19);
    absOrRelToggle->SetBounds(12, h - lineHeightDiff + 15, 105, 19);

    lineHeightDiff += 25;
    convEpsText->SetBounds(12, h - lineHeightDiff + 15, 170, 19);
    convEpsField->SetBounds(182, h - lineHeightDiff + 15, 30, 19);
    convBandLenText->SetBounds(220, h - lineHeightDiff + 15, 145, 19);
    convBandLenField->SetBounds(368, h - lineHeightDiff + 15, 30, 19);
    convApplyButton->SetBounds(406, h - lineHeightDiff + 15, 100, 19);
    shapeParamField->SetBounds(w-60, h - lineHeightDiff + 15, 50, 19);
    shapeParamField2->SetBounds(w-115, h - lineHeightDiff + 15, 50, 19);

    lineHeightDiff += 25;
    formulaText->SetBounds(7, h - lineHeightDiff, 350, 19);
    formulaBtn->SetBounds(360, h - lineHeightDiff, 120, 19);;
    dismissButton->SetBounds(w - 100, h - lineHeightDiff, 90, 19);

    lineHeightDiff += 25;
    colorToggle->SetBounds(7, h - lineHeightDiff, 105, 19);
    fixedLineWidthText->SetBounds(112, h - lineHeightDiff, 93, 19);
    fixedLineWidthField->SetBounds(206, h - lineHeightDiff, 30, 19);
    fixedLineWidthButton->SetBounds(240, h - lineHeightDiff, 100, 19);
    pruneEveryNButton->SetBounds(w - 215, h - lineHeightDiff, 100, 19);
    pruneFirstNButton->SetBounds(w - 110, h - lineHeightDiff, 100, 19);

    lineHeightDiff += 25;
    profCombo->SetBounds(7, h - lineHeightDiff, 160, 19);
    logYToggle->SetBounds(180, h - lineHeightDiff, 40, 19);
    addButton->SetBounds(w - 270, h - lineHeightDiff, 80, 19);
    removeButton->SetBounds(w - 180, h - lineHeightDiff, 80, 19);
    removeAllButton->SetBounds(w - 90, h - lineHeightDiff, 80, 19);

    chart->SetBounds(7, 5, w - 15, h - lineHeightDiff - 10);

    GLWindow::SetBounds(x, y, w, h);

}

/**
* \brief Resets convergence data e.g. when simulation resets
*/
void ConvergencePlotter::ResetData() {
    //Rebuild vector size
    for (auto &convVec : formula_ptr->convergenceValues) {
        convVec = ConvergenceData();
    }
}

/**
* \brief Refreshes all window components (combo, chart, values)
*/
void ConvergencePlotter::Refresh() {

    //Rebuild selection combo box
    size_t nbFormulas = formula_ptr->formulas_n.size(); // minimum 1 for custom input
    profCombo->Clear();
    if (nbFormulas) {
        profCombo->SetSize(nbFormulas);
        for (size_t i = 0; i < nbFormulas; i++) {
            char tmp[128];
            sprintf(tmp, "[%zd] %s", i + 1, formula_ptr->formulas_n[i]->GetExpression());
            profCombo->SetValueAt(i, tmp, (int) i);
        }
        profCombo->SetEditable(true);
    } else {
        profCombo->SetSize(1);
        char tmp[128];
        sprintf(tmp, "%s", "No formula found");
        profCombo->SetValueAt(0, tmp, (int) -1);
        profCombo->SetEditable(false);
    }
    profCombo->SetSelectedIndex(0);

    //Remove profiles that aren't present anymore
    for (int v = nbView - 1; v >= 0; v--) { //int because it can be -1, nbView is also int
        if (formula_ptr->formulas_n.empty()) {
            chart->GetY1Axis()->RemoveDataView(views[v]);
            SAFE_DELETE(views[v]);
            for (size_t j = v; j < nbView - 1; j++) views[j] = views[j + 1];
            nbView--;
            continue;
        }
        int formId = 0;
        while (views[v]->userData1 !=
               (int) std::hash<std::string>{}(formula_ptr->formulas_n[formId]->GetExpression())) {
            ++formId;
            if (formId >= formula_ptr->formulas_n.size()) {
                chart->GetY1Axis()->RemoveDataView(views[v]);
                SAFE_DELETE(views[v]);
                for (size_t j = v; j < nbView - 1; j++) views[j] = views[j + 1];
                nbView--;
                break;
            }
        }
    }

    //Update values
    refreshViews();
}


/**
* \brief Displays window with refreshed values
* \param w worker handle
*/
void ConvergencePlotter::Display(Worker *w) {


    SetWorker(w);
    formula_ptr->UpdateVectorSize();
    Refresh();
    SetVisible(true);

}

/**
* \brief Refreshes the view if needed
* \param appTime current time of the application
*/
void ConvergencePlotter::Update(float appTime) {
    bool hasChanged = formula_ptr->FetchNewConvValue();
    if (hasChanged || ((appTime - lastUpdate > 1.0f) && nbView)) {
        if (worker->IsRunning())
            refreshViews();
        lastUpdate = appTime;
        shapeParamField->SetText(formula_ptr->ApproxShapeParameter(0, 1));
        shapeParamField2->SetText(formula_ptr->ApproxShapeParameter(0,0));
    }
}

/**
* \brief Creates a plot from the expression given in the textbox of the form f(x)=EXPRESSION (e.g. 2*x+50)
*/
void ConvergencePlotter::plot() {

    GLParser *parser = new GLParser();
    parser->SetExpression(formulaText->GetText().c_str());
    if (!parser->Parse()) {
        GLMessageBox::Display(parser->GetErrorMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
        SAFE_DELETE(parser);
        return;
    }

    int nbVar = parser->GetNbVariable();
    if (nbVar == 0) {
        GLMessageBox::Display("Variable 'x' not found", "Error", GLDLG_OK, GLDLG_ICONERROR);
        SAFE_DELETE(parser);
        return;
    }
    if (nbVar > 1) {
        GLMessageBox::Display("Too much variables or unknown constant", "Error", GLDLG_OK, GLDLG_ICONERROR);
        SAFE_DELETE(parser);
        return;
    }
    VLIST *var = parser->GetVariableAt(0);
    if (!iequals(var->name, "x")) {
        GLMessageBox::Display("Variable 'x' not found", "Error", GLDLG_OK, GLDLG_ICONERROR);
        SAFE_DELETE(parser);
        return;
    }

    GLDataView *v;

    // Check that view is not already added
    bool found = false;
    int i = 0;
    while (i < nbView && !found) {
        found = (views[i]->userData1 == -1);
        if (!found) i++;
    }

    if (found) {
        v = views[i];
        v->SetName(formulaText->GetText().c_str());
        v->Reset();
    } else {

        if (nbView < 50) {
            v = new GLDataView();
            v->SetName(formulaText->GetText().c_str());
            v->userData1 = -1;
            chart->GetY1Axis()->AddDataView(v);
            views[nbView] = v;
            nbView++;
        } else {
            return;
        }
    }

    // Plot
    for (i = 0; i < 1000; i++) {
        double x = (double) i;
        double y;
        var->value = x;
        parser->Evaluate(&y);
        v->Add(x, y, false);
    }
    v->CommitChange();

    delete parser;

}

/**
* \brief Refreshes view by updating the data for the plot
*/
void ConvergencePlotter::refreshViews() {

    // Lock during update

    for (int i = 0; i < nbView; i++) {

        GLDataView *v = views[i];
        if (formula_ptr->formulas_n.empty()) return;
        int formId = 0;
        while (v->userData1 != (int) std::hash<std::string>{}(formula_ptr->formulas_n[formId]->GetExpression())) {
            ++formId;
            if (formId >= formula_ptr->formulas_n.size())
                return;
        }

        v->Reset();
        if (worker->globalHitCache.globalHits.nbDesorbed > 0) {
            auto& conv_vec = formula_ptr->convergenceValues[formId].conv_vec;
            for (int j = std::max(0,(int)conv_vec.size()-1000); j < conv_vec.size(); j++) // limit data points to last 1000
                v->Add(conv_vec[j].first, conv_vec[j].second, false);
        }
        v->CommitChange();

    }
}

/**
* \brief Adds view/plot for a specific facet
* \param formulaHash specific facet ID
* \return 0 if okay, 1 if already plotted
*/
int ConvergencePlotter::addView(int formulaHash) {

    char tmp[128];

    // Check that view is not already added
    bool found = false;
    int i = 0;
    while (i < nbView && !found) {
        found = (views[i]->userData1 == formulaHash);
        if (!found) i++;
    }
    if (found) {
        return 1;
    }
    if (nbView < MAX_VIEWS) {
        GLParser *formula = nullptr;
        for (GLParser *f : formula_ptr->formulas_n) {
            int str_hash = std::hash<std::string>{}(f->GetExpression());
            if (str_hash == formulaHash) {
                formula = f;
                break;
            }
        }
        if (!formula || formula_ptr->formulas_n.empty()) return 0;
        GLDataView *v = new GLDataView();
        sprintf(tmp, "%s", formula->GetExpression());
        v->SetName(tmp);
        //Look for first available color
        GLColor col = chart->GetFirstAvailableColor();
        int lineStyle = chart->GetFirstAvailableLinestyle(col);
        v->SetColor(col);
        v->SetMarkerColor(col);
        v->SetStyle(lineStyle);
        v->SetLineWidth(2);
        v->userData1 = formulaHash;

        chart->GetY1Axis()->AddDataView(v);
        views[nbView] = v;
        nbView++;
    }

    return 0;
}

/**
* \brief Removes view/plot for a specific facet
* \param formulaHash specific facet ID
* \return 0 if okay, 1 if not plotted
*/
int ConvergencePlotter::remView(int formulaHash) {

    bool found = false;
    int i = 0;
    while (i < nbView && !found) {
        found = (views[i]->userData1 == formulaHash);
        if (!found) i++;
    }
    if (!found) {
        return 1;
    }
    chart->GetY1Axis()->RemoveDataView(views[i]);
    SAFE_DELETE(views[i]);
    for (int j = i; j < nbView - 1; j++) views[j] = views[j + 1];
    nbView--;

    return 0;
}

/**
* \brief Resets the whole chart
*/
void ConvergencePlotter::Reset() {

    chart->GetY1Axis()->ClearDataView();
    for (int i = 0; i < nbView; i++) SAFE_DELETE(views[i]);
    nbView = 0;
}

/**
* \brief Function for processing various inputs (button, check boxes etc.)
* \param src Exact source of the call
* \param message Type of the source (button)
*/
void ConvergencePlotter::ProcessMessage(GLComponent *src, int message) {

    switch (message) {
        case MSG_BUTTON:
            if (src == dismissButton) {
                SetVisible(false);
            } else if (src == pruneEveryNButton) {
                for (int formulaId = 0; formulaId < formula_ptr->formulas_n.size(); ++formulaId) {
                    formula_ptr->pruneEveryN(4, formulaId, 0);
                }
                refreshViews();
            } else if (src == pruneFirstNButton) {
                for (int formulaId = 0; formulaId < formula_ptr->formulas_n.size(); ++formulaId) {
                    formula_ptr->pruneFirstN(100, formulaId);
                }
                refreshViews();
            } else if (src == addButton) {
                int idx = profCombo->GetSelectedIndex();
                if (idx >= 0 && !formula_ptr->formulas_n.empty()) { //Something selected (not -1)
                    if(formula_ptr->formulas_n.at(profCombo->GetUserValueAt(idx))->hasVariableEvalError){
                        GLMessageBox::Display("Formula can't be evaluated.", "Error", GLDLG_OK, GLDLG_ICONERROR);
                        break;
                    }
                    int str_hash = std::hash<std::string>{}(
                            formula_ptr->formulas_n.at(profCombo->GetUserValueAt(idx))->GetExpression());
                    if (addView(str_hash))
                        GLMessageBox::Display("Profile already plotted", "Info", GLDLG_OK, GLDLG_ICONINFO);
                    refreshViews();
                }
            } else if (src == removeButton) {

                int idx = profCombo->GetSelectedIndex();
                if (idx >= 0 && !formula_ptr->formulas_n.empty()) { //Something selected (not -1)
                    int str_hash = std::hash<std::string>{}(
                            formula_ptr->formulas_n.at(profCombo->GetUserValueAt(idx))->GetExpression());
                    if (remView(str_hash)) {
                        GLMessageBox::Display("Profile not plotted", "Error", GLDLG_OK, GLDLG_ICONERROR);
                    }
                    refreshViews();
                }
            } else if (src == removeAllButton) {
                Reset();
            } else if (src == formulaBtn) {
                plot();
            } else if (src == fixedLineWidthButton) {
                int linW;
                fixedLineWidthField->GetNumberInt(&linW);
                for (int viewId = 0; viewId < nbView; viewId++) {
                    GLDataView *v = views[viewId];
                    v->SetLineWidth(linW);
                }
            } else if (src == convApplyButton) {
                int newEps, newCBLen;
                convEpsField->GetNumberInt(&newEps);
                convBandLenField->GetNumberInt(&newCBLen);

                if(newEps != formula_ptr->epsilon || newCBLen != formula_ptr->cb_length){
                    formula_ptr->epsilon = newEps;
                    formula_ptr->cb_length = newCBLen;
                    for(int id = 0; id < formula_ptr->convergenceValues.size(); ++id)
                        formula_ptr->RestartASCBR(id);
                }
            }
            break;
        case MSG_COMBO:
            if (src == profCombo) {
                int profMode = profCombo->GetSelectedIndex();
            }
            break;
        case MSG_TOGGLE:
            if (src == logYToggle) {
                chart->GetY1Axis()->SetScale(logYToggle->GetState());
            } else if (src == colorToggle) {
                if (!colorToggle->GetState())
                    chart->SetColorSchemeDefault();
                else
                    chart->SetColorSchemeColorblind();

                const auto &colors = chart->GetColorScheme();
                for (int viewId = 0; viewId < nbView; viewId++) {
                    GLDataView *v = views[viewId];
                    auto col = colors[viewId % colors.size()];
                    int lineStyle = chart->GetFirstAvailableLinestyle(col);
                    v->SetColor(col);
                    v->SetMarkerColor(col);
                    v->SetStyle(lineStyle);

                    std::string facId(v->GetName());
                    facId = facId.substr(facId.find('#') + 1);
                }
            } else if(src == absOrRelToggle){
                formula_ptr->useAbsEps = absOrRelToggle->GetState();
                for(int id = 0; id < formula_ptr->convergenceValues.size(); ++id)
                    formula_ptr->RestartASCBR(id);
            }
            break;
    }
    //this->worker->GetGeometry()->SetPlottedFacets(plottedFacets);
    GLWindow::ProcessMessage(src, message);

}

/**
* \brief Adds views to the plotter if loaded form a file (XML)
* \param updatedViews vector containing the ids of the views
*/
void ConvergencePlotter::SetViews(const std::vector<int> &updatedViews) {
    Reset();
    for (int view : updatedViews)
        addView(view);
    //Refresh(); //Commented out: at this point, simulation results are not yet loaded
}

/**
* \brief Create and return a vector of view IDs
* \return vector containing the IDs of the views
*/
std::vector<int> ConvergencePlotter::GetViews() const {
    std::vector<int> v;
    v.reserve(nbView);
    for (size_t i = 0; i < nbView; i++)
        v.push_back(views[i]->userData1);
    return v;
}

/**
* \brief Returns bool for active/inactive logarithmic scaling
* \return bool that expresses if Y axis is logarithmically scaled
*/
bool ConvergencePlotter::IsLogScaled() {
    return chart->GetY1Axis()->GetScale();
}

/**
* \brief Sets type of scale (logarithmic or not)
* \param bool that expresses if Y axis should be logarithmically scaled or not
*/
void ConvergencePlotter::SetLogScaled(bool logScale) {
    chart->GetY1Axis()->SetScale(logScale);
    logYToggle->SetState(logScale);
}

/**
* \brief Sets worker handle for loading views before the full geometry
* \param w worker handle
*/
void ConvergencePlotter::SetWorker(Worker *w) { //for loading views before the full geometry

    worker = w;

}