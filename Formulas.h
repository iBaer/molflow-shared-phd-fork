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

#ifndef MOLFLOW_PROJ_FORMULAS_H
#define MOLFLOW_PROJ_FORMULAS_H

#include <memory>
#include <vector>
#include "GLApp/GLParser.h"
#include "FormulaEvaluator.h"

struct ConvergenceData {
    ConvergenceData() : conv_total(0.0), n_samples(0), upper_bound(0.0), lower_bound(0.0), chain_length(0) {};
    std::vector<std::pair<size_t,double>> conv_vec;
    double conv_total; /* for now unused accumulator, where convegence values are summed up */
    size_t n_samples;

    // ASCBR values
    double upper_bound;
    double lower_bound;
    size_t chain_length;
};

//! Defines a formula object that can be used to retrieve and store a parsed result
struct Formulas {

    Formulas(FormulaEvaluator* eval) : formulasChanged(true), sampleConvValues(true), epsilon(5), cb_length(51), useAbsEps(true){
        evaluator=eval;
    };
    ~Formulas(){
        for(auto f : formulas_n){
            SAFE_DELETE(f);
        }
        delete evaluator;
    };

    void AddFormula(const char *fName, const char *formula);
    void ClearFormulas();

    void UpdateVectorSize();
    bool InitializeFormulas();
    bool UpdateFormulaValues(size_t nbDesorbed);
    bool FetchNewConvValue();
    double GetConvRate(int formulaId);
    void RestartASCBR(int formulaId);
    bool CheckASCBR(int formulaId);
    double ApproxShapeParameter(int formulaId, int index_from);

    void pruneEveryN(size_t everyN, int formulaId, size_t skipLastN);
    void pruneFirstN(size_t n, int formulaId);

    std::vector<GLParser*> formulas_n;
    std::vector<std::pair<size_t,double>> lastFormulaValues;
    std::vector<ConvergenceData> convergenceValues; // One vector of nbDesorption,formulaValue pairs for each formula
    std::vector<std::vector<size_t>> freq_accum;
    bool formulasChanged;
    bool sampleConvValues;

    bool useAbsEps;
    int epsilon;
    size_t cb_length;
    FormulaEvaluator* evaluator;
};

#endif //MOLFLOW_PROJ_FORMULAS_H
