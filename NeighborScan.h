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

#ifndef MOLFLOW_PROJ_NEIGHBORSCAN_H
#define MOLFLOW_PROJ_NEIGHBORSCAN_H

#include "Geometry_shared.h"
#include "FacetData.h"

struct NeighborFacet;
struct CommonEdge {
    CommonEdge(int id, int v1, int v2) : facetId(), v1(v1), v2(v2){
        facetId.emplace_back(id);
        angle = 0.0;
        swapped = false;
    };
    CommonEdge(size_t id, size_t v1, size_t v2) : facetId(), v1(v1), v2(v2){
        facetId.emplace_back(static_cast<int>(id));
        angle = 0.0;
        swapped = false;
    };
    CommonEdge(int id1, int id2, int v1, int v2) : facetId(), v1(v1), v2(v2){
        facetId.emplace_back(id1);
        facetId.emplace_back(id2);
        angle = 0.0;
        swapped = false;
    };
    void Merge(CommonEdge& e2){
        for(auto id : e2.facetId) {
            if (std::find(facetId.begin(), facetId.end(), id) == facetId.end()) {
                // someName not in name, add it
                facetId.push_back(id);
            }
        }
        std::sort(facetId.begin(),facetId.end());
    }

    std::vector<int> facetId; // which facet owns edge
    int v1; // index vector 1
    int v2; // index vector 2
    bool swapped;
    double angle;
};

struct OverlappingEdge {
    OverlappingEdge(int id, int v1, int v2) : facetId1(id), f1v1(v1), f1v2(v2){
        angle = 0.0;
    };
    OverlappingEdge(size_t id, size_t v1, size_t v2) : facetId1(id), f1v1(v1), f1v2(v2){
        angle = 0.0;
    };
    OverlappingEdge(int id1, int id2, int f1v1, int f1v2, int f2v1, int f2v2)
    : facetId1(id1), f1v1(f1v1), f1v2(f1v2), facetId2(id2), f2v1(f2v1), f2v2(f2v2){
        angle = 0.0;
    };

    int facetId1; // which facet owns edge
    int facetId2; // which facet owns edge
    int f1v1; // index vector 1
    int f1v2; // index vector 2
    int f2v1; // index vector 1
    int f2v2; // index vector 2

    double angle;
};

class NeighborScan {
public:

    static int GetAnalysedCommonEdges(const std::vector<Facet *>& facets, std::vector<CommonEdge> &commonEdges);
    static int GetAnalysedUnorientedCommonEdges(const std::vector<Facet *>& facets, std::vector<CommonEdge> &commonEdges);
    static int GetAnalysedOverlappingEdges(const std::vector<Facet *> &facets, const std::vector<Vector3d> &vectors, std::vector<OverlappingEdge> &commonEdges);
    static int GetCommonEdgesList(const std::vector<Facet *> &facets, std::vector<CommonEdge> &commonEdges);
    static int GetUnorientedCommonEdgesList(const std::vector<Facet *>& facets, std::vector<CommonEdge> &commonEdges);
    static int GetOverlappingEdges(const std::vector<Facet *> &facets, const std::vector<Vector3d> &vectors,
                                   std::vector<OverlappingEdge> &overlappingEdges);

    static int CompareOverlap(std::vector<OverlappingEdge>& edges_overlap, std::vector<CommonEdge>& edges_un);
};

int RemoveDuplicates(std::vector<CommonEdge>& edge_v);
std::optional<double> AngleBetween2Vertices(Vector3d& v1, Vector3d& v2);

#endif //MOLFLOW_PROJ_NEIGHBORSCAN_H