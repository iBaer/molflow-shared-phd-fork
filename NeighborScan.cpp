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

#include <Helper/MathTools.h>
#include <Helper/ConsoleLogger.h>
#include <Helper/Chronometer.h>
#include "NeighborScan.h"

#include <cmath>

std::optional<double> AngleBetween2Vertices(Vector3d &v1, Vector3d &v2) {
    double denum = (v1.Length() * v2.Length());

    double numerator = std::clamp(Dot(v1, v2), -1.0, 1.0);

    if (denum == 0.0) {
        //Log::console_error("[NeighborAnalysis] Neighbors found with invalid angle: %lf/%lf\n", numerator, denum);
        return std::nullopt;
    }
    double angle_calc = numerator / (denum);
    if (angle_calc - 1e-7 <= 1.0 && angle_calc + 1e-7 >= -1.0) {
        angle_calc = std::clamp(angle_calc, -1.0, 1.0);
    } else {
        //Log::console_error("[NeighborAnalysis] Neighbors (%d , %d) found with invalid angle: %lf/%lf = %lf\n", numerator, denum, angle_calc);
        return std::nullopt;
    }

    return std::acos(angle_calc);
}

template<template<typename, typename> class Container,
        typename CommonEdge,
        typename Allocator=std::allocator<CommonEdge> >
void CalculateNeighborAngles(Container<CommonEdge, Allocator> &edges, std::vector<Facet *> facets) {
    for (auto iter_o = edges.begin(); iter_o != edges.end();) {
        auto f = facets[iter_o->facetId[0]];
        auto g = facets[iter_o->facetId[1]];
        auto angle_opt = AngleBetween2Vertices(f->sh.N, g->sh.N);
        if (!angle_opt.has_value()) {
            Log::console_error("[NeighborAnalysis] Neighbors found with invalid angle: %d , %d\n", iter_o->facetId[0],
                               iter_o->facetId[1]);
            iter_o = edges.erase(iter_o);
            continue;
        }

        iter_o->angle = angle_opt.value();
        ++iter_o;
    }
}

template<template<typename, typename> class Container,
        typename OverlappingEdge,
        typename Allocator=std::allocator<OverlappingEdge> >
void CalculateNeighborAnglesOverlapped(Container<OverlappingEdge, Allocator> &edges, std::vector<Facet *> facets) {
    for (auto iter_o = edges.begin(); iter_o != edges.end();) {
        auto f = facets[iter_o->facetId1];
        auto g = facets[iter_o->facetId2];
        auto angle_opt = AngleBetween2Vertices(f->sh.N, g->sh.N);
        if (!angle_opt.has_value()) {
            Log::console_error("[NeighborAnalysis] Neighbors found with invalid angle: %d , %d\n", iter_o->facetId1,
                               iter_o->facetId2);
            iter_o = edges.erase(iter_o);
            continue;
        }

        iter_o->angle = angle_opt.value();
        ++iter_o;
    }
}

/*!
 * @brief Filters common edges that are oriented differently to each other
 * Checks for orientation with `swapped`, merges two edges together to combine facet id information
 * Requires: Edges in order of (v1 < v2) and secondary prime facet id in order a.facetId[0] < b.facetId[0]
 * @tparam Container template container type, works for
 * @tparam CommonEdge specific data type we are accessing from within an arbitrary container
 * @tparam Allocator memory allocator for specific data type @CommonEdge
 * @param edges e.g. list or vector of common edges
 * @todo in this version depending on the facet number order there can be duplicates, leading to an endless loop
 */
/*template <template <typename, typename> class Container,
        typename CommonEdge,
        typename Allocator=std::allocator<CommonEdge> >
void CombineEdges(Container<CommonEdge, Allocator>& edges){
    int nextAdd = 1;
    auto stop_size = std::pow(edges.size(), 1.3);
    Container<CommonEdge, Allocator> edge_cpy;
    for(auto iter_o = edges.begin(); iter_o != edges.end(); ){
        auto iter_next = std::next(iter_o,nextAdd);
        // limit to prevent smart selection getting stucked
        if(edges.size() > stop_size) {
            edges.clear();
            return;
        }
        if(23432 == iter_o->v1 && 23432 == iter_next->v1
        && (iter_o->facetId[0] == 30599 || iter_next->facetId[0] == 30599)){
            fmt::print("Debug check\n");
        }
        if(iter_next == edges.end()) {
            if(iter_o->facetId.size()<=1) {
                iter_o = edges.erase(iter_o);
            }
            else{
                ++iter_o;
            }
            nextAdd = 1;
            continue;
        }
        if(iter_o->v1 == iter_next->v1 && iter_o->v2 == iter_next->v2){
            if(iter_o->facetId[0] == iter_next->facetId[0] || iter_o->facetId[1] == iter_next->facetId[0]) {
                nextAdd++;
                continue;
            }
            else if(iter_o->swapped != iter_next->swapped) {
                // common edge found
                if(iter_o->facetId.size() <= 1)
                    iter_o->Merge(*iter_next);
                else{
                    //create new edge
                    CommonEdge cpy(*iter_o);
                    cpy.facetId.resize(1);
                    cpy.Merge(*iter_next);
                    edge_cpy.push_back(cpy);
                }
                nextAdd++;
                continue;
            }
            else {
                nextAdd++;
                continue;
            }
        }
        *//*else if(iter_o->facetId[0] == iter_next->facetId[0]) {
            nextAdd++;
            continue;
        }*//*
        else {
            if(iter_o->facetId.size()<=1)
                iter_o = edges.erase(iter_o);
            else if(!edge_cpy.empty()){
                if(std::next(iter_o,nextAdd+1) == edges.end()) {
                    for (auto &cpy: edge_cpy)
                        edges.push_back(cpy);
                    break;
                }
                else {
                    iter_o = edges.insert(std::next(iter_o, 1), edge_cpy.begin(), edge_cpy.end());
                    ++iter_o;
                }
                edge_cpy.clear();
            }
            else {
                ++iter_o;
            }
        }
        nextAdd = 1;
    }
}*/

int RemoveDuplicates(std::vector<CommonEdge> &edge_v) {
#if defined(DEBUG)
    for(auto& e : edge_v) {
            if (e.facetId.size() < 1) {
                fmt::print("Debug check\n");
            }
        }
#endif
    std::sort(edge_v.begin(), edge_v.end(),
              [](const CommonEdge &e0, const CommonEdge &e1) -> bool {
                  if (e0.facetId[0] == e1.facetId[0])
                      return e0.facetId[1] < e1.facetId[1];
                  else
                      return e0.facetId[0] < e1.facetId[0];
              });

    auto last = std::unique(edge_v.begin(), edge_v.end(), [](const CommonEdge &e0, const CommonEdge &e1) -> bool {
        if (e0.facetId[0] == e1.facetId[0])
            return e0.facetId[1] == e1.facetId[1];
        else
            return false;
    });

    int nRemoved = std::distance(last, edge_v.end());
    edge_v.erase(last, edge_v.end());

    return nRemoved;
}


int RemoveDuplicates(std::vector<OverlappingEdge> &edge_v) {
    std::sort(edge_v.begin(), edge_v.end(),
              [](const OverlappingEdge &e0, const OverlappingEdge &e1) -> bool {
                  if (e0.facetId1 == e1.facetId1)
                      return e0.facetId2 < e1.facetId2;
                  else
                      return e0.facetId1 < e1.facetId1;
              });

    auto last = std::unique(edge_v.begin(), edge_v.end(),
                            [](const OverlappingEdge &e0, const OverlappingEdge &e1) -> bool {
                                if (e0.facetId1 == e1.facetId1)
                                    return e0.facetId2 == e1.facetId2;
                                else
                                    return false;
                            });

    int nRemoved = std::distance(last, edge_v.end());
    edge_v.erase(last, edge_v.end());

    return nRemoved;
}

template<template<typename, typename> class Container,
        typename CommonEdge,
        typename Allocator=std::allocator<CommonEdge> >
void CombineEdges(Container<CommonEdge, Allocator> &edges) {
    int nextAdd = 1;
    auto stop_size = std::pow(edges.size(), 1.3);
    Container<CommonEdge, Allocator> edge_cpy;
    // loop over all edges, combine when shared edge is found, if multiple shared facets are found, create an external extra facet for later merging into main list
    // alternatively, allow to merge 3>= facets into list, and save order info
    for (auto iter_o = edges.begin(); iter_o != edges.end();) {
        auto iter_next = std::next(iter_o, nextAdd);
        // limit to prevent smart selection getting stucked
        if (edges.size() > (unsigned int) stop_size) {
            edges.clear();
            return;
        }
        if (iter_next == edges.end()) {
            if (iter_o->facetId.size() <= 1) {
                iter_o = edges.erase(iter_o);
            } else {
                ++iter_o;
            }
            nextAdd = 1;
            continue;
        }
        if (iter_o->v1 == iter_next->v1 && iter_o->v2 == iter_next->v2) {
            if (iter_o->facetId[0] == iter_next->facetId[0] ||
                (iter_o->facetId.size() >= 2 && iter_o->facetId[1] == iter_next->facetId[0])) {
                nextAdd++;
                continue;
            } else if (iter_o->swapped != iter_next->swapped) {
                // common edge found
                if (iter_o->facetId.size() <= 1)
                    iter_o->Merge(*iter_next);
                else {
                    //create new edge
                    CommonEdge cpy(*iter_o);
                    cpy.facetId.resize(1);
                    cpy.Merge(*iter_next);
                    edge_cpy.push_back(cpy);
                }
                nextAdd++;
                continue;
            } else {
                nextAdd++;
                continue;
            }
        }
            /*else if(iter_o->facetId[0] == iter_next->facetId[0]) {
                nextAdd++;
                continue;
            }*/
        else {
            if (iter_o->facetId.size() <= 1)
                iter_o = edges.erase(iter_o);
            else {
                ++iter_o;
            }
        }
        nextAdd = 1;
    }

    // Remove potential duplicates in extra edge cpy list
    if (!edge_cpy.empty()) {
        // Find duplicates
        auto last = std::unique(edge_cpy.begin(), edge_cpy.end(),
                                [](const CommonEdge &e0, const CommonEdge &e1) -> bool {
                                    if (e0.facetId[0] == e1.facetId[0])
                                        return e0.facetId[1] == e1.facetId[1];
                                    else
                                        return false;
                                });
        // Delete duplicate elements and insert new edges into list
        edge_cpy.erase(last, edge_cpy.end());
        edges.insert(edges.end(), edge_cpy.begin(), edge_cpy.end());
    }
}

template<template<typename, typename> class Container,
        typename CommonEdge,
        typename Allocator=std::allocator<CommonEdge> >
void CombineUnorientedEdges(Container<CommonEdge, Allocator> &edges) {
    int nextAdd = 1;
    auto stop_size = std::pow(edges.size(), 1.3);
    Container<CommonEdge, Allocator> edge_cpy;
    // loop over all edges, combine when shared edge is found, if multiple shared facets are found, create an external extra facet for later merging into main list
    // alternatively, allow to merge 3>= facets into list, and save order info
    for (auto iter_o = edges.begin(); iter_o != edges.end();) {
        auto iter_next = std::next(iter_o, nextAdd);
        // limit to prevent smart selection getting stucked
        if (edges.size() > (unsigned int) stop_size) {
            edges.clear();
            return;
        }
        if (iter_next == edges.end()) {
            if (iter_o->facetId.size() <= 1) {
                iter_o = edges.erase(iter_o);
            } else {
                ++iter_o;
            }
            nextAdd = 1;
            continue;
        }
        if (iter_o->v1 == iter_next->v1 && iter_o->v2 == iter_next->v2) {
            if (iter_o->facetId[0] == iter_next->facetId[0] ||
                (iter_o->facetId.size() >= 2 && iter_o->facetId[1] == iter_next->facetId[0])) {
                nextAdd++;
                continue;
            } else if (iter_o->swapped != iter_next->swapped) {
                // common edge found
                if (iter_o->facetId.size() <= 1)
                    iter_o->Merge(*iter_next);
                else {
                    //create new edge
                    CommonEdge cpy(*iter_o);
                    cpy.facetId.resize(1);
                    cpy.Merge(*iter_next);
                    edge_cpy.push_back(cpy);
                }
                nextAdd++;
                continue;
            } else {
                nextAdd++;
                continue;
            }
        } else {
            if (iter_o->facetId.size() <= 1)
                iter_o = edges.erase(iter_o);
            else {
                ++iter_o;
            }
        }
        nextAdd = 1;
    }

    // Remove potential duplicates in extra edge cpy list
    if (!edge_cpy.empty()) {
        // Find duplicates
        auto last = std::unique(edge_cpy.begin(), edge_cpy.end(),
                                [](const CommonEdge &e0, const CommonEdge &e1) -> bool {
                                    if (e0.facetId[0] == e1.facetId[0])
                                        return e0.facetId[1] == e1.facetId[1];
                                    else
                                        return false;
                                });
        // Delete duplicate elements and insert new edges into list
        edge_cpy.erase(last, edge_cpy.end());
        edges.insert(edges.end(), edge_cpy.begin(), edge_cpy.end());
    }
}


// TODO: Could return a list of neighbors directly tuple{id1, id2, angle}
int NeighborScan::GetAnalysedCommonEdges(const std::vector<Facet *> &facets, std::vector<CommonEdge> &commonEdges) {
    int res = GetCommonEdgesList(facets, commonEdges);
    CalculateNeighborAngles(commonEdges, facets);
    return res;
}

int NeighborScan::GetAnalysedUnorientedCommonEdges(const std::vector<Facet *> &facets,
                                                   std::vector<CommonEdge> &commonEdges) {
    int res = GetUnorientedCommonEdgesList(facets, commonEdges);
    CalculateNeighborAngles(commonEdges, facets);
    return res;
}

int NeighborScan::GetAnalysedOverlappingEdges(const std::vector<Facet *> &facets, const std::vector<Vector3d> &vectors,
                                              std::vector<OverlappingEdge> &commonEdges) {
    int res = GetOverlappingEdges(facets, vectors, commonEdges);
    CalculateNeighborAnglesOverlapped(commonEdges, facets);
    return res;
}

int NeighborScan::GetCommonEdgesList(const std::vector<Facet *> &facets, std::vector<CommonEdge> &commonEdges) {

    // Detect common edge between facet
    size_t p11, p12;

    // 1. Create list of edges
    std::list<CommonEdge> edges;
    for (size_t facetId = 0; facetId < facets.size(); facetId++) {
        auto facet = facets[facetId];
        for (size_t i = 0; i < facet->sh.nbIndex; i++) {
            if (i == facet->sh.nbIndex - 1) {
                // GetIndex will turn last (i+1) to 0
                p11 = facet->indices[i];
                p12 = facet->indices[0];
            } else {
                p11 = facet->indices[i];
                p12 = facet->indices[i + 1];
            }
            bool swapped = false;
            if (p11 > p12) {// keep order for easier sort
                std::swap(p11, p12);
                swapped = true;
            }
            auto &edge = edges.emplace_back(facetId, p11, p12);
            edge.swapped = swapped;
        }
    }

    // 2. Sort edges
    edges.sort(
            [](const CommonEdge &e0, const CommonEdge &e1) -> bool {
                if (e0.v1 == e1.v1)
                    if (e0.v2 == e1.v2)
                        return (int) e0.facetId[0] < (int) e1.facetId[0];
                    else
                        return (int) e0.v2 < (int) e1.v2;
                else
                    return e0.v1 < e1.v1;
            });

    // 3. Get pairs
    CombineEdges(edges);

    commonEdges.clear();
    for (auto &edge: edges) {
        commonEdges.emplace_back(edge);
    }
    RemoveDuplicates(commonEdges);
    //HandleLoneEdge(commonEdges);

    return commonEdges.size();

}

// Same as GetCommonEdges, but orientation does not matter (vertex order)
int
NeighborScan::GetUnorientedCommonEdgesList(const std::vector<Facet *> &facets, std::vector<CommonEdge> &commonEdges) {

    // Detect common edge between facet
    size_t p11, p12;

    // 1. Create list of edges
    std::list<CommonEdge> edges;
    for (size_t facetId = 0; facetId < facets.size(); facetId++) {
        auto facet = facets[facetId];
        for (size_t i = 0; i < facet->sh.nbIndex; i++) {
            p11 = facet->indices[i];
            if (i == facet->sh.nbIndex - 1) {
                // GetIndex will turn last (i+1) to 0
                p12 = facet->indices[0];
            } else {
                p12 = facet->indices[i + 1];
            }
            bool swapped = false;
            if (p11 > p12) {// keep order for easier sort
                std::swap(p11, p12);
                swapped = true;
            }
            auto &edge = edges.emplace_back(facetId, p11, p12);
            edge.swapped = swapped;
        }
    }

    // 2. Sort edges
    edges.sort(
            [](const CommonEdge &e0, const CommonEdge &e1) -> bool {
                if (e0.v1 == e1.v1)
                    if (e0.v2 == e1.v2)
                        return (int) e0.facetId[0] < (int) e1.facetId[0];
                    else
                        return (int) e0.v2 < (int) e1.v2;
                else
                    return e0.v1 < e1.v1;
            });

    // 3. Get pairs
    CombineUnorientedEdges(edges);

    commonEdges.clear();
    for (auto &edge: edges) {
        commonEdges.emplace_back(edge);
    }
    RemoveDuplicates(commonEdges);
    //HandleLoneEdge(commonEdges);

    return commonEdges.size();

}


struct FacetVectors {
    FacetVectors(size_t id, Vector3d d, Vector3d v1, Vector3d v2, int vid1, int vid2) : facetId(id), direction(d),
                                                                                        v1(v1), v2(v2), vid1(vid1),
                                                                                        vid2(vid2) {};
    size_t facetId;
    Vector3d direction;
    Vector3d v1, v2;
    int vid1, vid2;
};

double CosineSimilarity(const Vector3d &a, const Vector3d &b) {
    double numerator = a.x * b.x + a.y * b.y + a.z * b.z;
    double denominator = a.Length() * b.Length();

    if (denominator < 1e-6) return 1.0f;
    return numerator / denominator;
}

double DistancePointToLine(const Vector3d &a, const Vector3d &b, const Vector3d &p) {
    Vector3d dir_ba = b - a;

    double ba_len = dir_ba.Length();
    auto crossbaap = CrossProduct(b - a, a - p);
    if (ba_len != 0.0) {
        return crossbaap.Length() / (b - a).Length();
    } else {
        return 1e42;
    }
    /*Vector3d dir_ba = b - a;
    double ba_len = dir_ba.Length();
    if(ba_len != 0.0){
        return CrossProduct(b-a, a-p).Length() / (b-a).Length();
    }
    else{
        return 1e42;
    }*/
}

int NeighborScan::GetOverlappingEdges(const std::vector<Facet *> &facets, const std::vector<Vector3d> &vectors,
                                      std::vector<OverlappingEdge> &overlappingEdges) {

    Chronometer timer;
    timer.Start();

    std::list<FacetVectors> vectors_xsort;

    // Detect overlapping edges to find geometric neighbors
    size_t p11, p12;
    Vector3d v1, v2;
    // 1. Define edge/line by direction + starting point
    for (size_t facetId = 0; facetId < facets.size(); facetId++) {
        auto facet = facets[facetId];
        for (size_t i = 0; i < facet->sh.nbIndex; i++) { // GetIndex will turn last (i+1) to 0
            p11 = facet->indices[i];
            p12 = (i + 1 >= facet->sh.nbIndex) ? facet->indices[0] : facet->indices[i + 1];
            v1 = vectors[p11];
            v2 = vectors[p12];

            Vector3d dir((v2 - v1).Normalized());
            if (dir.x < 0) {
                // Keep x direction positive for sortability
                // For overlapping lines, the orientation does not matter
                dir = (-1.0) * dir;
            }
            vectors_xsort.emplace_back(facetId, dir, v1, v2, p11, p12);
        }
    }

    Log::console_msg(3, "#ModelReader: Defined edges {}s\n", timer.Elapsed());
    timer.ReInit(); timer.Start();

    // Sort _edges_ for _axis_
    vectors_xsort.sort(
            [](const FacetVectors &lhs, const FacetVectors &rhs) -> bool {
                auto &dir0 = lhs.direction;
                auto &dir1 = rhs.direction;

                if (dir0.x < dir1.x)
                    return true;
                else if (dir0.x == dir1.x) {
                    if (dir0.y < dir1.y) {
                        return true;
                    } else if (dir0.y == dir1.y) {
                        if (dir0.z < dir1.z) {
                            return true;
                        }
                    }
                }

                return false;
            });

    Log::console_msg(3, "#ModelReader: Sorted edges {}s\n", timer.Elapsed());
    timer.ReInit(); timer.Start();

    long unsigned int t_c = 0;
    long unsigned int t_c2 = 0;
    long unsigned int s_c = 0;
    long unsigned int s_c2 = 0;

    int x_c = 0;
    for (auto x_iter = vectors_xsort.begin(); x_iter != vectors_xsort.end(); x_iter++) {
        x_c++;
        //int would_skip = 0;
        int o_c = 0;
        for (auto iter_o = std::next(x_iter); iter_o != vectors_xsort.end(); iter_o++) {
            o_c++;
            t_c++;
            if (x_iter->facetId == iter_o->facetId)
                continue;

            if (std::fabs(x_iter->direction.x - iter_o->direction.x) >= 1e-5) {
                break;
                //would_skip += std::fabs(x_iter->direction.x - iter_o->direction.x) >= 1e-4 ? 1 : 0;
                //Log::console_msg(0, "SKIP AT Sim {} for x {} ... {}\n", sim, x_iter->direction.x, iter_o->direction.x);
            }

            t_c2++;

            double sim = CosineSimilarity(x_iter->direction, iter_o->direction);

            //https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html>v1
            if (IsEqual(std::fabs(sim), 1, 1e-6)) {
                s_c++;

                int nOverlaps = 0;
                if (IsEqual(DistancePointToLine(x_iter->v1, x_iter->v2, iter_o->v1), 0.0)) nOverlaps++;
                if (IsEqual(DistancePointToLine(x_iter->v1, x_iter->v2, iter_o->v2), 0.0)) nOverlaps++;
                if (IsEqual(DistancePointToLine(iter_o->v1, iter_o->v2, x_iter->v1), 0.0)) nOverlaps++;
                if (IsEqual(DistancePointToLine(iter_o->v1, iter_o->v2, x_iter->v2), 0.0)) nOverlaps++;
                if (nOverlaps >= 2) {
                    /*fmt::print("Vectors are overlapping and similar [{}] {} , {}, {} -- {} , {} , {}\n", sim, x_iter->direction.x,
                               x_iter->direction.y, x_iter->direction.z,
                               iter_o->direction.x, iter_o->direction.y, iter_o->direction.z);*/
                    // Insert by lower facet ID first
                    if (x_iter->facetId < iter_o->facetId)
                        overlappingEdges.emplace_back(x_iter->facetId, iter_o->facetId, x_iter->vid1, x_iter->vid2,
                                                      iter_o->vid1, iter_o->vid2);
                    else {
                        overlappingEdges.emplace_back(iter_o->facetId, x_iter->facetId, iter_o->vid1, iter_o->vid2,
                                                      x_iter->vid1, x_iter->vid2);
                    }

                    /*if (would_skip >= 1)
                    Log::console_msg(0, "[{}] SSim {} for x {} ... {}\n", would_skip, sim, x_iter->direction.x,
                                     iter_o->direction.x);*/
                }
                /*else if(nOverlaps == 1)
                {
                    fmt::print("Distances are {} , {} x {} , {}\n", DistancePointToLine(x_iter->v1, x_iter->v2, iter_o->v1)
                            ,DistancePointToLine(x_iter->v1, x_iter->v2, iter_o->v2)
                            ,DistancePointToLine(iter_o->v1, iter_o->v2, x_iter->v1)
                            ,DistancePointToLine(iter_o->v1, iter_o->v2, x_iter->v2));
                    fmt::print("Vectors are not overlapping but similar [{}] {} , {}, {} -- {} , {} , {} x {} , {}, {} -- {} , {} , {}\n", sim,
                               x_iter->v1.x, x_iter->v1.y, x_iter->v1.z,
                               x_iter->v2.x, x_iter->v2.y, x_iter->v2.z,
                               iter_o->v1.x, iter_o->v1.y, iter_o->v1.z,
                               iter_o->v2.x, iter_o->v2.y, iter_o->v2.z);
                }*/
            }
            else {
                s_c2++;
            }
        }
    }
    Log::console_msg(3, "#ModelReader: c{} , {} s{} , {}: Found {}/{} overlapping edges {}s\n", t_c, t_c2, s_c, s_c2, overlappingEdges.size(), vectors_xsort.size(), timer.Elapsed());
    timer.ReInit(); timer.Start();

    //int nRemoved = RemoveDuplicates(overlappingEdges);


    /*
    RemoveDuplicates(overlappingEdges);
    //HandleLoneEdge(overlappingEdges);
*/
    return overlappingEdges.size();

}

int NeighborScan::CompareOverlap(std::vector<OverlappingEdge>& edges_overlap, std::vector<CommonEdge>& edges_un) {
    Chronometer timer;
    timer.Start();

    // Sort _edges_ for facet id
    std::sort(edges_overlap.begin(), edges_overlap.end(),
              [](const OverlappingEdge &lhs, const OverlappingEdge &rhs) -> bool {
                  auto &dir0 = lhs.facetId1;
                  auto &dir1 = rhs.facetId1;

                  if (lhs.facetId1 < rhs.facetId1)
                      return true;
                  else if (lhs.facetId1 == rhs.facetId1) {
                      if (lhs.facetId2 < rhs.facetId2) {
                          return true;
                      }
                  }

                  return false;
              });

    std::sort(edges_un.begin(), edges_un.end(),
              [](const CommonEdge &lhs, const CommonEdge &rhs) -> bool {

                  if (lhs.facetId[0] < rhs.facetId[0])
                      return true;
                  else if (lhs.facetId[0] == rhs.facetId[0]) {
                      if (lhs.facetId[1] < rhs.facetId[1]) {
                          return true;
                      }
                  }

                  return false;
              });

    /*for (auto edgeo_it = edges_overlap.begin(); edgeo_it != edges_overlap.end();) {
        for (auto edge_it = edges_un.begin(); edge_it != edges_un.end(); edge_it++) {
            if (edge_it->facetId[0] == edgeo_it->facetId1 && edge_it->facetId[1] == edgeo_it->facetId2) {
                edgeo_it = edges_overlap.erase(edgeo_it);
                //edge_it = edges_un.begin();
                if (edgeo_it == edges_overlap.end()) {
                    break;
                }
            }
        }
        if (edgeo_it != edges_overlap.end())
            edgeo_it++;
    }*/

    auto edgeo_it = edges_overlap.begin();
    auto edge_it = edges_un.begin();

    unsigned long int n_loops = 0;
    for (; edgeo_it != edges_overlap.end() && edge_it != edges_un.end();) {
        if (edge_it->facetId[0] == edgeo_it->facetId1 && edge_it->facetId[1] == edgeo_it->facetId2) {
            edgeo_it = edges_overlap.erase(edgeo_it);
            //edge_it = edges_un.begin();
            continue;
        }

        if(edge_it->facetId[0] <= edgeo_it->facetId1) {
            edge_it++;
        } else {
            edgeo_it++;
        }

        n_loops++;
    }

    Log::console_msg(3, "#ModelReader: Compared {} overlapping edges {}s\n", n_loops, timer.Elapsed());

    return edges_overlap.size();
}