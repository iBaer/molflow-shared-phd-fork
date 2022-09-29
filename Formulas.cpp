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

#include "Formulas.h"
#include "Helper/ConsoleLogger.h"
#include <cmath>

// convergence constants
constexpr size_t max_vector_size() { return 65536; };
constexpr size_t d_precision() { return 5; };

//! Add a formula to the formula storage
void Formulas::AddFormula(const char *fName, const char *formula) {
    GLParser *f2 = new GLParser();
    f2->SetExpression(formula);
    f2->SetName(fName);
    f2->Parse();

    formulas_n.push_back(f2);
    UpdateVectorSize();
}

//! Clear formula storage
void Formulas::ClearFormulas() {
    for (auto &f : formulas_n)
        SAFE_DELETE(f);
    formulas_n.clear();
    UpdateVectorSize();
}

//! Initialize formulas by parsing string to values
bool Formulas::InitializeFormulas(){
    bool allOk = true;
    for (auto & i : formulas_n) {

        // Evaluate variables
        int nbVar = i->GetNbVariable();
        bool ok = true;
        for (int j = 0; j < nbVar && ok; j++) {
            VLIST *v = i->GetVariableAt(j);
            ok = evaluator->EvaluateVariable(v);
        }

        if (ok) {
            i->hasVariableEvalError = false;
        }
        else{
            i->SetVariableEvalError("Unknown formula specifier");
            allOk = false;
        }
    }
    return allOk;
};

/**
* \brief Resizes Vector storing convergence values if needed (change of expressions, etc.)
*/
void Formulas::UpdateVectorSize() {
    //Rebuild vector size
    size_t nbFormulas = formulas_n.size();
    if (convergenceValues.size() != nbFormulas) {
        convergenceValues.resize(nbFormulas);
    }
    if (lastFormulaValues.size() != nbFormulas) {
        lastFormulaValues.resize(nbFormulas);
    }
    if (freq_accum.size() != nbFormulas) {
        freq_accum.resize(nbFormulas, std::vector<size_t>(cb_length));
    }

}

//! Get updated formula values for later usage in corresponding windows
bool Formulas::UpdateFormulaValues(size_t nbDesorbed) {
    // First fetch new variable values
    InitializeFormulas();

    // Next evaluate each formula and cache values for later usage
    // in FormulaEditor, ConvergencePlotter
    if (!formulas_n.empty()) {
        for (int formulaId = 0; formulaId < formulas_n.size(); ++formulaId) {
            double r;
            if (formulas_n[formulaId]->Evaluate(&r)) {
                lastFormulaValues[formulaId] = std::make_pair(nbDesorbed, r);
            }
            else{
                formulas_n[formulaId]->SetVariableEvalError(formulas_n.at(formulaId)->GetErrorMsg());
            }
        }
        FetchNewConvValue();
    }

    return true;
}

//! Add new value to convergence vector
bool Formulas::FetchNewConvValue() {
    bool hasChanged = false;
    if (!formulas_n.empty()) {
        UpdateVectorSize();
        for (int formulaId = 0; formulaId < formulas_n.size(); ++formulaId) {
            if(formulas_n.at(formulaId)->hasVariableEvalError) {
                continue;
            }

            auto& currentValues = lastFormulaValues[formulaId];
            // Check against potential nan values to skip
            if(std::isnan(currentValues.second))
                continue;
            auto& conv_vec = convergenceValues[formulaId].conv_vec;
            if (conv_vec.size() >= max_vector_size()) {
                pruneEveryN(4, formulaId, 1000); // delete every 4th element for now
                hasChanged = true;
            }
            // Insert new value when completely new value pair inserted
            if(conv_vec.empty() || (currentValues.first != conv_vec.back().first && currentValues.second != conv_vec.back().second)){
                convergenceValues[formulaId].n_samples += 1;
                convergenceValues[formulaId].conv_total += currentValues.second;
                conv_vec.emplace_back(currentValues);
                hasChanged = true;
            }
            else if (currentValues.first == conv_vec.back().first && currentValues.second != conv_vec.back().second){
                // if (for some reason) the nbDesorptions dont change but the formula value, only update the latter
                convergenceValues[formulaId].conv_total -= conv_vec.back().second + currentValues.second;
                conv_vec.back().second = currentValues.second;
                hasChanged = true;
            }
            else if (currentValues.second == conv_vec.back().second){
                if(conv_vec.size() > 1 && currentValues.second == (conv_vec.rbegin()[1]).second){
                    // if the value remains constant, just update the nbDesorbptions
                    conv_vec.back().first = currentValues.first;
                    hasChanged = true;
                }
                else {
                    // if there was no previous point with the same value, add a new one (only a second one to save space)
                    convergenceValues[formulaId].n_samples += 1;
                    convergenceValues[formulaId].conv_total += currentValues.second;
                    conv_vec.emplace_back(currentValues);
                    hasChanged = true;
                }
            }
        }
    }

    return hasChanged;
};

/**
* \brief Removes every everyN-th element from the convergence vector in case the max size has been reached
* \param everyN specifies stepsize for backwards delete
* \param formulaId formula whose convergence values shall be pruned
* \param skipLastN skips last N data points e.g. to retain a sharp view, while freeing other data points
*/
void Formulas::pruneEveryN(size_t everyN, int formulaId, size_t skipLastN) {
    auto &conv_vec = convergenceValues[formulaId].conv_vec;
    for (int i = conv_vec.size() - everyN - skipLastN; i > 0; i = i - everyN)
        conv_vec.erase(conv_vec.begin() + i);
}

/**
* \brief Removes first n elements from the convergence vector
 * \param n amount of elements that should be removed from the front
 * \param formulaId formula whose convergence values shall be pruned
*/
void Formulas::pruneFirstN(size_t n, int formulaId) {
    auto &conv_vec = convergenceValues[formulaId].conv_vec;
    conv_vec.erase(conv_vec.begin(), conv_vec.begin() + std::min(n, conv_vec.size()));
}

/**
* \brief Removes first n elements from the convergence vector
 * \param formulaId formula whose convergence values shall be pruned
*/
double Formulas::GetConvRate(int formulaId) {

    auto &conv_vec = convergenceValues[formulaId].conv_vec;
    if(conv_vec.empty()) return 0.0;

    double sumValSquared = 0.0;
    double sumVal = 0.0;
    for(auto& convVal : conv_vec){
        sumValSquared += convVal.second * convVal.second;
        sumVal += convVal.second;
    }
    double invN = 1.0 / conv_vec.size();
    double var = invN * sumValSquared - std::pow(invN * sumVal,2);
    double convDelta = 3.0 * std::sqrt(var) / (conv_vec.size() * conv_vec.size());


    //double xmean = sumVal * invN;
    double xmean = conv_vec.back().second; // ^ it doesnt make sense to calculate a mean, as it' s already commulative
    double varn = 0.0;
    for(auto& convVal : conv_vec){
        varn += std::pow(convVal.second - xmean,2);
    }
    varn *= 1.0 / (conv_vec.size() - 1.0);
    double quant95Val = 1.96 * std::sqrt(varn) / std::sqrt(conv_vec.size());
    double quant99Val = 2.576 * std::sqrt(varn) / std::sqrt(conv_vec.size());
    double quant995Val = 3.0 * std::sqrt(varn) / std::sqrt(conv_vec.size());

    const double val_eps = useAbsEps ? 0.5 * std::pow(10.0, -1.0*epsilon) : 0.5 * xmean * std::pow(10.0, -1.0*epsilon);
    double upper_bound = xmean + val_eps;
    double lower_bound = xmean - val_eps;

   /* if(convDelta / xmean < 1.0e-6)
        Log::console_msg(2, "[1] Sufficient convergent reached: {:e}\n", convDelta / conv_vec.back().second);

*/
    if(xmean - quant995Val >= lower_bound && xmean + quant995Val <= upper_bound)
        Log::console_msg(2, "[CI 99.5%] Convergence reached: {:e} +- {:e} / {:e} = {:e} in [{:e} , {:e}]\n",
             conv_vec.back().second, varn, std::sqrt(varn), quant995Val, xmean - quant99Val, xmean + quant995Val);
    else if(xmean - quant99Val >= lower_bound && xmean + quant99Val <= upper_bound)
        Log::console_msg(2, "[CI 99.0%] Convergence reached: {:e} in [{:e} , {:e}]\n", conv_vec.back().second, xmean - quant99Val, xmean + quant99Val);
    else if(xmean - quant95Val >= lower_bound && xmean + quant95Val <= upper_bound)
        Log::console_msg(2, "[CI 95.0%] Convergence reached: {:e} in [{:e} , {:e}]\n", conv_vec.back().second, xmean - quant95Val, xmean + quant95Val);
        /*if(xmean - quant995Val <= conv_vec.back().second && xmean + quant995Val >= conv_vec.back().second)
            Log::console_msg(2, "[CI 99.5%] Convergence reached: {:e} +- {:e} / {:e} = {:e} in [{:e} , {:e}]\n",
                             conv_vec.back().second, varn, std::sqrt(varn), quant995Val, xmean - quant99Val, xmean + quant995Val);
        else if(xmean - quant99Val <= conv_vec.back().second && xmean + quant99Val >= conv_vec.back().second)
            Log::console_msg(2, "[CI 99.0%] Convergence reached: {:e} in [{:e} , {:e}]\n", conv_vec.back().second, xmean - quant99Val, xmean + quant99Val);
        else if(xmean - quant95Val <= conv_vec.back().second && xmean + quant95Val >= conv_vec.back().second)
            Log::console_msg(2, "[CI 95.0%] Convergence reached: {:e} in [{:e} , {:e}]\n", conv_vec.back().second, xmean - quant95Val, xmean + quant95Val);
        */
    else {
        double dist = std::min(lower_bound - xmean + quant95Val, upper_bound - xmean - quant95Val);
        Log::console_msg(2, "[4] Convergence distance to a=0.95: {:e} --> {:e} close to [{:e} , {:e}]\n", dist, conv_vec.back().second, xmean - quant95Val, xmean + quant95Val);
        //Log::console_msg(2, "[4] Abs to a95: {:e} --> Rel to a95: {:e} / {:e}\n", quant95Val, (conv_vec.back().second - xmean) / xmean, (quant95Val / conv_vec.back().second));
    }
    return quant95Val * sqrt(varn*varn/conv_vec.size());
}

/**
* \brief Restart values for ASCBR in case method parameters changed
*/
void Formulas::RestartASCBR(int formulaId){
    auto& convData = convergenceValues[formulaId];
    convData.chain_length = 0;
    convData.upper_bound = convData.lower_bound = 0.0;

    freq_accum[formulaId].clear();
    freq_accum[formulaId].resize(cb_length);
}

/**
* \brief Check stopping criterium based on *acceptable shifting convergence band rule* (ASCBR)
 * \param formulaId formula whose conv values should be checked for
 * \return true = end condition met: run algorithm until (chain length == cb_length)
*/
bool Formulas::CheckASCBR(int formulaId) {
    // Initialize
    auto& convData = convergenceValues[formulaId];
    if(convData.conv_vec.empty()) return false;

    const size_t cb_len = cb_length; // convergence band length

    // Step 1: increment step, +1 MC event
    // done outside

    // Step 2: Increment value for variable total T
    // done outside in FetchNewConvValue()

    // & calculate new mean
    //double conv_mean_local = convData.conv_total / convData.n_samples;
    const double conv_mean_local = convData.conv_vec.back().second;

    // Step 3: Check if mean still within bounds
    bool withinBounds = (conv_mean_local >= convData.lower_bound && conv_mean_local <= convData.upper_bound);

    if(!withinBounds){
        // half-width of convergence band
        const double eps = useAbsEps ? 0.5 * std::pow(10.0, -1.0*epsilon) : 0.5 * conv_mean_local * std::pow(10.0, -1.0*epsilon);
        // absolute: exponent d = significant digits after decimal point

        // Step 4: Update bounds
        convData.upper_bound = conv_mean_local + eps;
        convData.lower_bound = conv_mean_local - eps;
        // Step 5: Update length of inbound "chain"
        ++freq_accum[formulaId][(convData.chain_length>=cb_length)?0:convData.chain_length];
        convData.chain_length = 0;
    }
    else{
        // Step 5: Update length of inbound "chain"
        ++convData.chain_length;
    }

    return convData.chain_length >= cb_len;
}

double Formulas::ApproxShapeParameter(int formulaId, int index_from = 1) {
    // Initialize
    double shape_param = 0.0;
    double den = 0.0;
    for(int i = 1; i < cb_length; ++i){
        den += (double) i * freq_accum[formulaId][i];
    }
    if(den <= 1e-8) den = 1.0;
    shape_param = freq_accum.size() > index_from ? 1.0 - freq_accum[formulaId][index_from] / den : 0.0;

    return shape_param;
}