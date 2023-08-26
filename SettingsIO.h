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

#ifndef MOLFLOW_PROJ_SETTINGSIO_H
#define MOLFLOW_PROJ_SETTINGSIO_H

#include <string>
#include <list>
#include <vector>
#include <map>
#include "GeometryTypes.h" // selecitons

// Settings for the CLI application for IO related things
// Also adds some functions to handle IO pathes/files
namespace SettingsIO {
    extern std::string workFile; //! Uncompressed file that is used for the simulation
    extern std::string workPath; //! Output path for various simulation files (e.g. autosave), can be different from outputPath e.g. for cluster simulations
    extern std::string inputFile; //! Uncompressed or compressed file that is used for the simulation
    extern std::string inputPath; //! Path to input file
    extern std::string outputFile; //! Output file name
    extern std::string outputPath; //! Outputpath for output file
    extern std::vector<std::string> extraFiles; //! deprecated
    extern std::map<std::string, std::vector<std::string>> cachedLines; //! deprecated
    extern std::vector<std::vector<std::string>> formulas; //! cached formula values to keep them for a valid output file
    extern std::vector<SelectionGroup> selections; //! cached selection groups to keep them for a valid output file

    extern bool overwrite; //! whether overwriting the inputfile is allowed
    extern bool isArchive; //! whether the input file is an archive
    extern bool outputFacetDetails; //! whether output for all facet details is wanted
    extern bool outputFacetQuantities; //! whether output for derived facet quantities is wanted
    extern double autogenerateTest; //! whether an automatically generated test case should be used as input

    int prepareIO();
    int initDirectories();
    int initFromZip();
    void cleanup_files();
}


#endif //MOLFLOW_PROJ_SETTINGSIO_H
