//
// Created by pbahr on 18/06/2020.
//

#ifndef MOLFLOW_PROJ_STRINGHELPER_H
#define MOLFLOW_PROJ_STRINGHELPER_H

#include <string>
#include <vector>

void splitFacetList(std::vector<size_t>& outputFacetIds, std::string inputString, size_t nbFacets);
std::string AbbreviateString(const std::string& input, size_t maxLength);
std::vector<std::string> SplitString(std::string const& input);
std::vector<std::string> SplitString(std::string const& input, const char& delimiter);

bool endsWith(std::string const& fullString, std::string const& ending);
bool beginsWith(std::string const& fullString, std::string const& beginning);
std::string space2underscore(std::string text);
bool iequals(std::string a, std::string b);


#endif //MOLFLOW_PROJ_STRINGHELPER_H
