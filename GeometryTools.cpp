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
#include "GeometryTools.h"
#include "Facet_shared.h"
#include "NeighborScan.h"
#include "SimulationModel.h"

#if defined(USE_CGAL)
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds> Delaunay;
typedef K::Point_2 Point;
#endif

std::vector<InterfaceFacet*> GeometryTools::GetTriangulatedGeometry(Geometry* geometry, std::vector<size_t> facetIndices, GLProgress* prg)
{
    std::vector<InterfaceFacet*> triangleFacets;
    for (size_t i = 0; i < facetIndices.size(); i++) {
        if (prg) prg->SetProgress((double)i/(double(facetIndices.size())));
        size_t nb = geometry->GetFacet(facetIndices[i])->sh.nbIndex;
        if (nb > 3) {
            // Create new triangle facets (does not invalidate old ones, you have to manually delete them)
            std::vector<InterfaceFacet*> newTriangles = Triangulate(geometry->GetFacet(facetIndices[i]));
            triangleFacets.insert(std::end(triangleFacets), std::begin(newTriangles), std::end(newTriangles));
        }
        else {
            //Copy
            auto* newFacet = new InterfaceFacet(nb);
            newFacet->indices = geometry->GetFacet(facetIndices[i])->indices;
            newFacet->CopyFacetProperties(geometry->GetFacet(facetIndices[i]), false);
            triangleFacets.push_back(newFacet);
        }
    }
    return triangleFacets;
}

#ifndef USE_CGAL
std::vector<InterfaceFacet*> GeometryTools::Triangulate(InterfaceFacet *f) {

    // Triangulate a facet (rendering purpose)
    // The facet must have at least 3 points
    // Use the very simple "Two-Ears" theorem. It computes in O(n^2).
    std::vector<InterfaceFacet*> triangleFacets;
    if (f->nonSimple) {
        // Not a simple polygon
        // Abort triangulation
        return triangleFacets;
    }

    // Build a Polygon
    GLAppPolygon p;
    p.pts = f->vertices2;
    //p.sign = f->sign;

    std::unique_ptr<InterfaceFacet> facetCopy(new InterfaceFacet(f->sh.nbIndex)); //Create a copy and don't touch original
    facetCopy->CopyFacetProperties(f);
    facetCopy->indices = f->indices;

    // Perform triangulation
    while (p.pts.size() > 3) {
        int e = FindEar(p);
        //DrawEar(f, p, e, addTextureCoord);

        // Create new triangle facet and copy polygon parameters, but change indices
        InterfaceFacet* triangle = GetTriangleFromEar(facetCopy.get(), p, e);
        triangleFacets.push_back(triangle);

        // Remove the ear
        p.pts.erase(p.pts.begin() + e);
    }

    // Draw the last ear
    InterfaceFacet* triangle = GetTriangleFromEar(facetCopy.get(), p, 0);
    triangleFacets.push_back(triangle);
    //DrawEar(f, p, 0, addTextureCoord);

    return triangleFacets;
}
#else
std::vector<InterfaceFacet*> GeometryTools::Triangulate(InterfaceFacet *f) {

    // Triangulate a facet (rendering purpose)
    // The facet must have at least 3 points
    // Use the very simple "Two-Ears" theorem. It computes in O(n^2).
    std::vector<InterfaceFacet*> triangleFacets;
    if (f->nonSimple) {
        // Not a simple polygon
        // Abort triangulation
        return triangleFacets;
    }

    // Build a Polygon
    GLAppPolygon p;
    p.pts = f->vertices2;
    //p.sign = f->sign;

    std::unique_ptr<InterfaceFacet> facetCopy(new InterfaceFacet(f->sh.nbIndex)); //Create a copy and don't touch original
    facetCopy->CopyFacetProperties(f);
    facetCopy->indices = f->indices;

    std::vector<std::pair<CGAL::Point_2<K>, unsigned>> indexed_vertices;
    for (unsigned i = 0; i < p.pts.size(); ++i) {
        indexed_vertices.emplace_back(CGAL::Point_2<K>(p.pts[i].u, p.pts[i].v), facetCopy->indices[i]);
    }


    // Perform Delaunay triangulation
    Delaunay triangulation;
    triangulation.insert(indexed_vertices.begin(), indexed_vertices.end());

    // Convert each face into a triangle and add to output list
    for (Delaunay::Finite_faces_iterator fit = triangulation.finite_faces_begin(); fit != triangulation.finite_faces_end(); ++fit) {
        Delaunay::Face_handle face = fit;
        int triangleIndices[3];

        triangleIndices[0] = face->vertex(0)->info();
        triangleIndices[1] = face->vertex(1)->info();
        triangleIndices[2] = face->vertex(2)->info();

    InterfaceFacet* triangle = new InterfaceFacet(3);
    triangle->CopyFacetProperties(facetCopy.get(),0);
    for (size_t i = 0; i < 3; i++) {
        triangle->indices[i] = triangleIndices[i];
    }
    triangleFacets.emplace_back(triangle);
    }

    return triangleFacets;
}
#endif
int  GeometryTools::FindEar(const GLAppPolygon& p){

    int i = 0;
    bool earFound = false;
    while (i < p.pts.size() && !earFound) {
        if (IsConvex(p, i))
            earFound = !ContainsConcave(p, i - 1, i, i + 1);
        if (!earFound) i++;
    }

    // REM: Theoretically, it should always find an ear (2-Ears theorem).
    // However on degenerated geometry (flat poly) it may not find one.
    // Returns first point in case of failure.
    if (earFound)
        return i;
    else
        return 0;
}

// Return indices to vertices3d for new triangle facet
InterfaceFacet* GeometryTools::GetTriangleFromEar(InterfaceFacet *f, const GLAppPolygon& p, int ear) {

    //Commented out sections: theoretically in a right-handed system the vertex order is inverse
    //However we'll solve it simpler by inverting the geometry viewer Front/back culling mode setting

    auto* indices = new size_t[3];

    indices[0] = f->indices[Previous(ear, p.pts.size())];
    indices[1] = f->indices[IDX(ear, p.pts.size())];
    indices[2] = f->indices[Next(ear, p.pts.size())];

    // Create new triangle facet and copy polygon parameters, but change indices
    InterfaceFacet* triangle = new InterfaceFacet(3);
    triangle->CopyFacetProperties(f,0);
    for (size_t i = 0; i < 3; i++) {
        triangle->indices[i] = indices[i];
    }
    f->indices.erase(f->indices.begin()+IDX(ear, p.pts.size()));
    delete[] indices;

/*    const Vector2d* p1 = &(p.pts[Previous(ear, p.pts.size())]);
    const Vector2d* p2 = &(p.pts[IDX(ear, p.pts.size())]);
    const Vector2d* p3 = &(p.pts[Next(ear, p.pts.size())]);*/

    return triangle;
}

/*std::optional<double> AngleBetween2Vertices(Vector3d& v1, Vector3d& v2){
    double denum = (v1.Norme() * v2.Length());

    double numerator = std::clamp(Dot(v1, v2), -1.0, 1.0);

    if(denum == 0.0){
        //Log::console_error("[NeighborAnalysis] Neighbors found with invalid angle: {}/{}\n", numerator, denum);
        return std::nullopt;
    }
    double angle_calc = numerator / (denum);
    if(angle_calc - 1e-7 <= 1.0 && angle_calc + 1e-7 >= -1.0) {
        angle_calc = std::clamp(angle_calc, -1.0, 1.0);
    }
    else{
        //Log::console_error("[NeighborAnalysis] Neighbors ({} , {}) found with invalid angle: {}/{} = {}\n", numerator, denum, angle_calc);
        return std::nullopt;
    }

    return std::acos(angle_calc);
}*/

template <template <typename, typename> class Container,
        typename CommonEdge,
        typename Allocator=std::allocator<CommonEdge> >
void CalculateNeighborAngles(Container<CommonEdge, Allocator>& edges, std::vector<Facet *> facets){
    for(auto iter_o = edges.begin(); iter_o != edges.end(); ){
        auto f = facets[iter_o->facetId[0]];
        auto g = facets[iter_o->facetId[1]];
        auto angle_opt =  AngleBetween2Vertices(f->sh.N, g->sh.N);
        if(!angle_opt.has_value()) {
            Log::console_error("[NeighborAnalysis] Neighbors found with invalid angle: {} , {}\n", iter_o->facetId[0], iter_o->facetId[1]);
            iter_o = edges.erase(iter_o);
            continue;
        }

        iter_o->angle = angle_opt.value();
        ++iter_o;
    }
}
/*

int RemoveDuplicates(std::vector<CommonEdge>& edge_v){
    for(auto& e : edge_v) {
        if (e.facetId.size() < 1) {
            fmt::print("Debug check\n");
        }
    }
    std::sort(&edge_v[0], &edge_v[edge_v.size()],
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

    int nRemoved = std::distance(last,edge_v.end());
    edge_v.erase(last, edge_v.end());

    return nRemoved;
}
*/

int HandleLoneEdge(std::vector<CommonEdge>& edge_v){
    int nRemoved = 0;
    for(auto iter_o = edge_v.begin(); iter_o != edge_v.end(); ){
        if((*iter_o).facetId.size() <= 1){
            Log::console_msg(2, "[{}] Lone edge : {} -- {}\n",(*iter_o).facetId[0],(*iter_o).v1, (*iter_o).v2);
            iter_o = edge_v.erase(iter_o);
            ++nRemoved;
        }
        else
            ++iter_o;
    }
    return nRemoved;
}

int HandleTransparent(std::vector<CommonEdge>& edge_v, Geometry* geometry){
    int nRemoved = 0;
    for(auto iter_o = edge_v.begin(); iter_o != edge_v.end(); ){
        auto f = geometry->GetFacet((*iter_o).facetId[0]);
        auto g = geometry->GetFacet((*iter_o).facetId[1]);

        if(f->sh.is2sided || g->sh.is2sided) {
            size_t c1, c2, l;
            if (!Geometry::GetCommonEdges(f, g, &c1, &c2, &l)) {
                iter_o = edge_v.erase(iter_o);
                ++nRemoved;
                continue;
            }
        }
        ++iter_o;
    }

    return nRemoved;
}

static std::vector<std::vector<CommonEdge>> edges_algo(6); // first index = facet id, second index = neighbor id

void GeometryTools::AnalyseGeometry(const std::shared_ptr<SimulationModel>& model) {
    edges_algo.resize(6);
    Chronometer stop;
    int dupli = 0, lone = 0, trans = 0;

    // original, but slow to compare
    /*stop.Start();
    AnalyseNeighbors(nullptr);
    printf("[{}] Neighbor T1: {}\n", edges_algo[0].size(), stop.Elapsed());
    //int dupli = RemoveDuplicates(edges_algo[0]);
    //int lone = HandleLoneEdge(edges_algo[0]);
    //int trans = HandleTransparent(edges_algo[0], geometry);
    printf("[{} - {} - {} - {}] Neighbor T1 Dupli: {}\n", edges_algo[0].size(), dupli, lone, trans, stop.Elapsed());
*/
    //stop.ReInit();
    stop.Start();
    GetCommonEdgesVec(nullptr, edges_algo[1]);
    Log::console_msg(2, "[{}] Neighbor T2 Vertex: {}\n", edges_algo[1].size(), stop.Elapsed());
    //dupli = RemoveDuplicates(edges_algo[1]);
    //lone = HandleLoneEdge(edges_algo[1]);

    auto& facets = model->facets;
    std::vector<Facet*> facet_ptr;
    facet_ptr.resize(facets.size());
    std::transform(facets.begin(), facets.end(), facet_ptr.begin(),
                   [](std::shared_ptr<SimulationFacet>& f){return (Facet*)(f.get());}
               );

    CalculateNeighborAngles(edges_algo[1], facet_ptr);
    //trans = HandleTransparent(edges_algo[1], geometry);
    Log::console_msg(2, "[{} - {} - {} - {}] Neighbor T2 Dupli: {}\n", edges_algo[1].size(), dupli, lone, trans, stop.Elapsed());
    edges_algo[0] = std::move(edges_algo[1]);

    stop.ReInit(); stop.Start();
    GetCommonEdgesList(facet_ptr, edges_algo[2]);
    Log::console_msg(2, "[{}] Neighbor T3 List: {}\n", edges_algo[2].size(), stop.Elapsed());
    //dupli = RemoveDuplicates(edges_algo[2]);
    //lone = HandleLoneEdge(edges_algo[2]);
    CalculateNeighborAngles(edges_algo[2], facet_ptr);
    //trans = HandleTransparent(edges_algo[2], geometry);
    Log::console_msg(2, "[{} - {} - {} - {}] Neighbor T3 Dupli: {}\n", edges_algo[2].size(), dupli, lone, trans, stop.Elapsed());

    stop.ReInit(); stop.Start();
    GetCommonEdgesHash(nullptr, edges_algo[3]);
    Log::console_msg(2, "[{}] Neighbor T4 Hashmap: {}\n", edges_algo[3].size(), stop.Elapsed());
    //dupli = RemoveDuplicates(edges_algo[3]);
    //lone = HandleLoneEdge(edges_algo[3]);
    //trans = HandleTransparent(edges_algo[3], geometry);
    CalculateNeighborAngles(edges_algo[3], facet_ptr);
    Log::console_msg(2, "[{} - {} - {} - {}] Neighbor T4 Dupli: {}\n", edges_algo[3].size(), dupli, lone, trans, stop.Elapsed());

    stop.ReInit(); stop.Start();
    GetCommonEdgesMap(model.get(), edges_algo[4]);
    Log::console_msg(2, "[{}] Neighbor T5 Map: {}\n", edges_algo[4].size(), stop.Elapsed());
    //dupli = RemoveDuplicates(edges_algo[3]);
    //lone = HandleLoneEdge(edges_algo[3]);
    //trans = HandleTransparent(edges_algo[3], geometry);
    CalculateNeighborAngles(edges_algo[4], facet_ptr);
    Log::console_msg(2, "[{} - {} - {} - {}] Neighbor T5 Dupli: {}\n", edges_algo[4].size(), dupli, lone, trans, stop.Elapsed());

    stop.ReInit(); stop.Start();
    GetCommonEdgesSingleVertex(nullptr, edges_algo[5]);
    Log::console_msg(2, "[{}] Neighbor T6 Map: {}\n", edges_algo[5].size(), stop.Elapsed());
    //dupli = RemoveDuplicates(edges_algo[3]);
    //lone = HandleLoneEdge(edges_algo[3]);
    //trans = HandleTransparent(edges_algo[3], geometry);
    //CalculateNeighborAngles(edges_algo[5], geometry);
    //printf("[%lu - %d - %d - %d] Neighbor T6 Dupli: %lf\n", edges_algo[5].size(), dupli, lone, trans, stop.Elapsed());


    // Skip and only compare with newer algorithms
    //CompareAlgorithm(nullptr, 1);
    CompareAlgorithm(nullptr, 2);
    CompareAlgorithm(nullptr, 3);
    CompareAlgorithm(model.get(), 4);
}

void GeometryTools::CompareAlgorithm(SimulationModel *model, size_t index) {
    // Put in order for comparism
    auto& compEdge = edges_algo[index];

    Log::console_msg(2, "[T1 vs. T{}] Neighbors found: {} - {}\n", index+1, edges_algo[0].size(), compEdge.size());

    int j = 0;
    size_t minCompSize = compEdge.size();
    for(int i = 0; i < edges_algo[0].size() && j < minCompSize; i++){
        bool diffSize = edges_algo[0][i].facetId.size() != compEdge[j].facetId.size() || edges_algo[0][i].facetId.size() != compEdge[j].facetId.size();
        bool diffVert = edges_algo[0][i].v1 != compEdge[j].v1 || edges_algo[0][i].v2 != compEdge[j].v2;
        bool diffFac = edges_algo[0][i].facetId[0] != compEdge[j].facetId[0] || edges_algo[0][i].facetId[1] != compEdge[j].facetId[1];
        bool diffAngle = edges_algo[0][i].angle != compEdge[j].angle;
        if(diffFac)
            Log::console_msg(2, "[{} , {}] {} -- {} <{} , {}> c({}) x {} -- {} <{} , {}> c({})\n", i, j,
                   edges_algo[0][i].facetId[0], edges_algo[0][i].facetId[1], edges_algo[0][i].v1, edges_algo[0][i].v2, edges_algo[0][i].angle * (180 / M_PI),
                   compEdge[j].facetId[0], compEdge[j].facetId[1], compEdge[j].v1, compEdge[j].v2, compEdge[j].angle * (180 / M_PI));
        if(diffSize){
            Log::console_msg(2, "[{} , {}] {} vs {} <{} , {}>: \n", i, j, edges_algo[0][i].facetId.size(), compEdge[j].facetId.size(), edges_algo[0][i].v1, edges_algo[0][i].v2);
            for(auto id : edges_algo[0][i].facetId) Log::console_msg(2, "{} , ", id);
            Log::console_msg(2, " --- ");
            for(auto id : compEdge[j].facetId) Log::console_msg(2, "{} , ", id);
            Log::console_msg(2, "\n");
        }
        if(diffAngle){
            Log::console_msg(2, "[{} , {}] {} -- {} <{} , {}> c({}) x {} -- {} <{} , {}> c({})\n", i, j,
                   edges_algo[0][i].facetId[0], edges_algo[0][i].facetId[1], edges_algo[0][i].v1, edges_algo[0][i].v2, edges_algo[0][i].angle * (180 / M_PI),
                   compEdge[j].facetId[0], compEdge[j].facetId[1], compEdge[j].v1, compEdge[j].v2, compEdge[j].angle * (180 / M_PI));
        }
        if(diffFac){
            if(edges_algo[0][i].facetId[0] > compEdge[j].facetId[0]) {
                j++;
                i--;
            }
        }
        else{
            j++;
        }
    }
}

// Update facet list of geometry by removing polygon facets and replacing them with triangular facets with the same properties
void GeometryTools::PolygonsToTriangles(Geometry* geometry) {
    auto allIndices = geometry->GetAllFacetIndices();
    std::vector<InterfaceFacet*> triangleFacets = GetTriangulatedGeometry(geometry , allIndices);
    geometry->RemoveFacets(allIndices);
    geometry->AddFacets(triangleFacets);
}

void GeometryTools::PolygonsToTriangles(Geometry* geometry, std::vector<size_t> selectedIndices) {
    std::vector<InterfaceFacet*> triangleFacets = GetTriangulatedGeometry(geometry, selectedIndices);
    geometry->RemoveFacets(selectedIndices);
    geometry->AddFacets(triangleFacets);
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

template <template <typename, typename> class Container,
        typename CommonEdge,
        typename Allocator=std::allocator<CommonEdge> >
void CombineEdges(Container<CommonEdge, Allocator>& edges){
    int nextAdd = 1;
    auto stop_size = std::pow(edges.size(), 1.3);
    Container<CommonEdge, Allocator> edge_cpy;
    // loop over all edges, combine when shared edge is found, if multiple shared facets are found, create an external extra facet for later merging into main list
    // alternatively, allow to merge 3>= facets into list, and save order info
    for(auto iter_o = edges.begin(); iter_o != edges.end(); ){
        auto iter_next = std::next(iter_o,nextAdd);
        // limit to prevent smart selection getting stucked
        if(edges.size() > (unsigned int) stop_size) {
            edges.clear();
            return;
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
            if (iter_o->facetId[0] == iter_next->facetId[0] || (iter_o->facetId.size() >= 2 && iter_o->facetId[1] == iter_next->facetId[0])) {
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
            /*else if(iter_o->facetId[0] == iter_next->facetId[0]) {
                nextAdd++;
                continue;
            }*/
        else {
            if(iter_o->facetId.size()<=1)
                iter_o = edges.erase(iter_o);
            else {
                ++iter_o;
            }
        }
        nextAdd = 1;
    }

    // Remove potential duplicates in extra edge cpy list
    if(!edge_cpy.empty()) {
        // Find duplicates
        auto last = std::unique(edge_cpy.begin(), edge_cpy.end(), [](const CommonEdge &e0, const CommonEdge &e1) -> bool {
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
int GeometryTools::GetAnalysedCommonEdges(std::vector<Facet *> facets, std::vector<CommonEdge> &commonEdges){
    int res = GetCommonEdgesList(facets, commonEdges);
    CalculateNeighborAngles(commonEdges, facets);
    return res;
}

std::vector<std::vector<NeighborFacet>> GeometryTools::AnalyseNeighbors(Geometry* geometry){

    // 1. First find neighbors
    // search for common edges
    size_t count = 0;
    std::vector<std::vector<NeighborFacet>> neighbors(geometry->GetNbFacet()); // first index = facet id, second index = neighbor id

    std::vector<CommonEdge> edges; // first index = facet id, second index = neighbor id
    for (int i = 0; i < geometry->GetNbFacet()-1; i++) {
        InterfaceFacet *f = geometry->GetFacet(i);
        for (int j = i+1; j < geometry->GetNbFacet(); j++) {
            InterfaceFacet *g = geometry->GetFacet(j);
            size_t  c1,c2,l;
            if(Geometry::GetCommonEdges(f,g, &c1,&c2,&l)){
                auto angle_opt = AngleBetween2Vertices(f->sh.N, g->sh.N);
                if(!angle_opt.has_value()) {
                    Log::console_error("[NeighborAnalysis] Neighbors found with invalid angle: {} , {}\n", i, j);
                    continue;
                }

                double angle = angle_opt.value();

                // neighbor found, add to both lists
                NeighborFacet n1{}, n2{};
                n1.angleDiff = angle;
                n1.id = i;
                n2.angleDiff = angle;
                n2.id = j;

                neighbors[i].push_back(n1);
                neighbors[j].push_back(n2);
                count++;

                CommonEdge e((size_t)i, f->GetIndex(c1), g->GetIndex(c2));
                e.facetId.emplace_back(j);
                if(e.v1 > e.v2)
                    std::swap(e.v1, e.v2);
                e.angle = angle;
                edges.emplace_back(e);
            }
        }
    }

    edges_algo[0] = edges;
    return neighbors;
}

int GeometryTools::GetCommonEdgesVec(Geometry *geometry, std::vector<CommonEdge> &commonEdges) {

    // Detect common edge between facet
    size_t p11, p12;

    // 1. Create list of edges
    commonEdges.clear();
    std::vector<CommonEdge>& edges = commonEdges;

    for (size_t facetId = 0; facetId < geometry->GetNbFacet(); facetId++) {
        auto facet = geometry->GetFacet(facetId);
        for (size_t i = 0; i < facet->sh.nbIndex; i++) { // GetIndex will turn last (i+1) to 0

            p11 = facet->GetIndex(i);
            p12 = facet->GetIndex(i + 1);

            bool swapped = false;
            if(p11 > p12) {// keep order for easier sort
                std::swap(p11, p12);
                swapped = true;
            }
            auto& edge = edges.emplace_back(facetId, p11, p12);
            edge.swapped = swapped;
        }
    }

    // 2. Sort edges
    std::sort(&edges[0], &edges[edges.size()],
              [](const CommonEdge &e0, const CommonEdge &e1) -> bool {
                  if (e0.v1 == e1.v1)
                      if(e0.v2 == e1.v2)
                          return (int)e0.facetId[0] < (int)e1.facetId[0];
                      else
                          return (int)e0.v2 < (int)e1.v2;
                  else
                      return e0.v1 < e1.v1;
              });

    // 3. Get pairs
    CombineEdges(edges);
    RemoveDuplicates(commonEdges);
    //HandleLoneEdge(commonEdges);

    return edges.size();
}

int GeometryTools::GetCommonEdgesList(std::vector<Facet *> facets, std::vector<CommonEdge> &commonEdges) {

    // Detect common edge between facet
    size_t p11, p12;

    // 1. Create list of edges
    std::list<CommonEdge> edges;
    for (size_t facetId = 0; facetId < facets.size(); facetId++) {
        auto facet = facets[facetId];
        for (size_t i = 0; i < facet->sh.nbIndex; i++) { // GetIndex will turn last (i+1) to 0
            p11 = facet->indices[i];
            p12 = facet->indices[i + 1];

            bool swapped = false;
            if(p11 > p12) {// keep order for easier sort
                std::swap(p11, p12);
                swapped = true;
            }
            auto& edge = edges.emplace_back(facetId, p11, p12);
            edge.swapped = swapped;
        }
    }

    // 2. Sort edges
    edges.sort(
              [](const CommonEdge &e0, const CommonEdge &e1) -> bool {
                  if (e0.v1 == e1.v1)
                      if(e0.v2 == e1.v2)
                          return (int)e0.facetId[0] < (int)e1.facetId[0];
                      else
                          return (int)e0.v2 < (int)e1.v2;
                  else
                      return e0.v1 < e1.v1;
              });

    // 3. Get pairs
    CombineEdges(edges);

    commonEdges.clear();
    for(auto& edge : edges){
        commonEdges.emplace_back(edge);
    }
    RemoveDuplicates(commonEdges);
    //HandleLoneEdge(commonEdges);

    return commonEdges.size();

}

int GeometryTools::GetCommonEdgesHash(Geometry *geometry, std::vector<CommonEdge> &commonEdges) {

    // Detect common edge between facet
    size_t p11, p12;

    // 1. Create list of edges
    commonEdges.clear();
    std::vector<CommonEdge>& edges = commonEdges;
    std::unordered_map<size_t,std::vector<int>> hashMap;

    for (size_t facetId = 0; facetId < geometry->GetNbFacet(); facetId++) {
        auto facet = geometry->GetFacet(facetId);
        for (size_t i = 0; i < facet->sh.nbIndex; i++) { // GetIndex will turn last (i+1) to 0

            p11 = facet->GetIndex(i);
            p12 = facet->GetIndex(i + 1);

            size_t hash = ((p11 << 16) | (p12));
            auto edge = hashMap.emplace(std::make_pair(hash, std::vector<int>{(int)facetId}));
            if(!edge.second) {// edge already exists, append facetid
                edge.first->second.push_back(facetId);
            }
        }
    }

    for(auto& entry : hashMap){
        p11 = entry.first >> 16;
        p12 = entry.first % (int)65536; // % 2^16

        if(p11 > p12) // only check one direction
            continue;
        auto rEdge = hashMap.find((p12 << 16) | (p11)); // look for reverse edge
        if(rEdge != hashMap.end()) {
            for(auto fac1 : entry.second) {
                for (auto fac2: rEdge->second) {
                    if(fac1 > fac2)
                        edges.emplace_back(CommonEdge(fac2, fac1, p11, p12));
                    else if(fac1 < fac2)
                        edges.emplace_back(CommonEdge(fac1, fac2, p11, p12));
                }
            }
        }
    }

    RemoveDuplicates(edges);
    return commonEdges.size();

}

int GeometryTools::GetCommonEdgesMap(SimulationModel *model, std::vector<CommonEdge> &commonEdges) {

    // Detect common edge between facet
    size_t p11, p12;

    // 1. Create list of edges
    commonEdges.clear();
    std::vector<CommonEdge>& edges = commonEdges;
    std::map<size_t,std::vector<int>> hashMap;

    for (size_t facetId = 0; facetId < model->facets.size(); facetId++) {
        auto facet = model->facets[facetId];
        for (size_t i = 0; i < facet->sh.nbIndex; i++) { // GetIndex will turn last (i+1) to 0

            p11 = facet->indices[i];
            p12 = facet->indices[i + 1];

            size_t hash = ((p11 << 16) | (p12));
            auto edge = hashMap.emplace(std::make_pair(hash, std::vector<int>{(int)facetId}));
            if(!edge.second) {// edge already exists, append facetid
                edge.first->second.push_back(facetId);
            }
        }
    }

    for(auto& entry : hashMap){
        p11 = entry.first >> 16;
        p12 = entry.first % (int)65536; // % 2^16

        if(p11 > p12) // only check one direction
            continue;
        auto rEdge = hashMap.find((p12 << 16) | (p11)); // look for reverse edge
        if(rEdge != hashMap.end()) {
            for(auto fac1 : entry.second) {
                for (auto fac2: rEdge->second) {
                    if(fac1 > fac2)
                        edges.emplace_back(CommonEdge(fac2, fac1, p11, p12));
                    else if(fac1 < fac2)
                        edges.emplace_back(CommonEdge(fac1, fac2, p11, p12));
                }
            }
        }
    }

    RemoveDuplicates(edges);
    return commonEdges.size();

}

int GeometryTools::GetCommonEdgesSingleVertex(Geometry *geometry, std::vector<CommonEdge> &commonEdges) {

    // 1. The first step is also to form the pairs of integers (i.e., two arrays of integers)
    std::vector<std::pair<int, int>> vertex_elements; // pairs of vertex IDs, element IDs

    for (size_t facetId = 0; facetId < geometry->GetNbFacet(); facetId++) {
        auto facet = geometry->GetFacet(facetId);
        for (size_t i = 0; i < facet->sh.nbIndex; i++) {
            vertex_elements.emplace_back(std::make_pair(facet->GetIndex(i), facetId));
        }
    }

    // 2. sort according to the first array of integers
    std::sort(&vertex_elements[0], &vertex_elements[vertex_elements.size()]);

    // 3. use the parallel segmented reduction and scan
    // to further determine both the total number
    // and the detailed indices of the neighboring elements

    commonEdges.clear();
    std::vector<CommonEdge>& edges = commonEdges;
    std::vector<int> positions;
    positions.resize(geometry->GetNbVertex());
    int pos = 0;
    positions[0] = 0;
    for(auto& entry : vertex_elements){
        auto vertex = entry.first;
        pos++;
        positions[vertex] = pos;
    }

    for(auto iter_o = vertex_elements.begin(); iter_o != vertex_elements.end(); ++iter_o){
        if((iter_o+1) == vertex_elements.end())
            break;
        auto vertex1 = iter_o->first;
        auto element1 = iter_o->second;

        for(auto iter_i = iter_o+1; iter_i != vertex_elements.end() && iter_o->first != iter_i->first && iter_o->second == iter_i->second; ++iter_i){
            auto vertex2 = iter_i->first;
            auto element2 = iter_i->second;

            auto& edge = edges.emplace_back(CommonEdge(element1, vertex1, vertex2));
            auto edge2 = CommonEdge(element2, vertex1, vertex2);
            edge.Merge(edge2);
        }
    }

    RemoveDuplicates(edges);
    return commonEdges.size();

}

/*

int NeighborScan::GetOverlappingEdges(std::vector<Facet *> facets, const std::vector<Vector3d>& vectors, std::vector<OverlappingEdge> &overlappingEdges) {

    std::list<FacetVectors> vectors_xsort;

    // Detect overlapping edges to find geometric neighbors
    size_t p11, p12;
    Vector3d v1, v2;
    // 1. Define edge/line by direction + starting point
    for (size_t facetId = 0; facetId < facets.size(); facetId++) {
        auto facet = facets[facetId];
        for (size_t i = 0; i < facet->sh.nbIndex; i++) { // GetIndex will turn last (i+1) to 0
            p11 = facet->indices[i];
            p12 = (i+1 >= facet->sh.nbIndex) ? facet->indices[0] : facet->indices[i + 1];
            v1 = vectors[p11];
            v2 = vectors[p12];

            Vector3d dir((v2-v1).Normalized());
            if(dir.x < 0){
                // Keep x direction positive for sortability
                // For overlapping lines, the orientation does not matter
                dir = (-1.0) * dir;
            }
            vectors_xsort.emplace_back(facetId, dir, v1, v2, p11, p12);
        }
    }

    // Sort _edges_ for _axis_
    vectors_xsort.sort(
            [](const FacetVectors &v0, const FacetVectors &v1) -> bool {
                return v0.direction.x < v1.direction.x;
            });

    for(auto x_iter = vectors_xsort.begin(); x_iter != vectors_xsort.end(); x_iter++){
        for(auto iter_o = std::next(x_iter); iter_o != vectors_xsort.end(); iter_o++){
            if(x_iter->facetId == iter_o->facetId)
                continue;
            double sim = CosineSimilarity(x_iter->direction, iter_o->direction);

            //https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html>v1
            if(IsEqual(std::fabs(sim), 1, 1e-6)){
                int nOverlaps = 0;
                if(IsEqual(DistancePointToLine(x_iter->v1, x_iter->v2, iter_o->v1), 0.0)) nOverlaps++;
                if(IsEqual(DistancePointToLine(x_iter->v1, x_iter->v2, iter_o->v2), 0.0)) nOverlaps++;
                if(IsEqual(DistancePointToLine(iter_o->v1, iter_o->v2, x_iter->v1), 0.0)) nOverlaps++;
                if(IsEqual(DistancePointToLine(iter_o->v1, iter_o->v2, x_iter->v2), 0.0)) nOverlaps++;
                if(nOverlaps >= 2)
                {
                    */
/*fmt::print("Vectors are overlapping and similar [{}] {} , {}, {} -- {} , {} , {}\n", sim, x_iter->direction.x,
                               x_iter->direction.y, x_iter->direction.z,
                               iter_o->direction.x, iter_o->direction.y, iter_o->direction.z);*//*

                    if(x_iter->facetId < iter_o->facetId)
                        overlappingEdges.emplace_back(x_iter->facetId, iter_o->facetId, x_iter->vid1, x_iter->vid2, iter_o->vid1, iter_o->vid2);
                    else {
                        overlappingEdges.emplace_back(iter_o->facetId, x_iter->facetId, iter_o->vid1, iter_o->vid2, x_iter->vid1, x_iter->vid2);
                    }
                }
                */
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
                }*//*

            }
        }
    }

    int nRemoved = RemoveDuplicates(overlappingEdges);


    */
/*
    RemoveDuplicates(overlappingEdges);
    //HandleLoneEdge(overlappingEdges);
*//*

    return overlappingEdges.size();

}
*/
