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

#include "SimulationModel.h"
#include "FacetData.h"
#include "SimulationFacet.h"
#include "Helper/MathTools.h"
#include <cmath>

size_t SimulationModel::size() {
    size_t modelSize = 0;
    modelSize += facets.capacity();
    for (auto &fac : facets)
        modelSize += fac->GetMemSize();
    modelSize += structures.capacity();
    for (auto &struc : structures)
        modelSize += struc.GetMemSize();
    modelSize += sizeof(std::vector<Vector3d>) + sizeof(Vector3d) * vertices3.capacity();
    //modelSize += tdParams.GetMemSize();
    modelSize += sizeof(otfParams);
    modelSize += sizeof(wp);
    modelSize += sizeof(sh);
    modelSize += sizeof(m);
    modelSize += sizeof(initialized);

    return modelSize;
}

/**
* \brief Initialises geometry properties that haven't been loaded from file
* \return error code: 0=no error, 1=error
*/
int SimulationModel::InitialiseFacets() {
    if (!m.try_lock()) {
        return 1;
    }

    for (const auto& f : facets) {
        auto& facet = *f;
        // Main facet params
        // Current facet
        //SubprocessFacet *f = model->facets[i];
        CalculateFacetParams(&facet);

        // Set some texture parameters
        // bool Facet::SetTexture(double width, double height, bool useMesh)
        if (facet.sh.texWidth_precise * facet.sh.texHeight_precise > 0.0000001) {
            const double ceilCutoff = 0.9999999;
            facet.sh.texWidth = (int) std::ceil(facet.sh.texWidth_precise *
                                                ceilCutoff); //0.9999999: cut the last few digits (convert rounding error 1.00000001 to 1, not 2)
            facet.sh.texHeight = (int) std::ceil(facet.sh.texHeight_precise * ceilCutoff);
        } else {
            facet.sh.texWidth = 0;
            facet.sh.texHeight = 0;
            facet.sh.texWidth_precise = 0.0;
            facet.sh.texHeight_precise = 0.0;
        }
    }

    m.unlock();
    return 0;
}

/*!
 * @brief Calculates various facet parameters without sanity checking @see Geometry::CalculateFacetParams(Facet* f)
 * @param f individual subprocess facet
 */
void SimulationModel::CalculateFacetParams(Facet* f) {
    // Calculate facet normal
    Vector3d p0 = vertices3[f->indices[0]];
    Vector3d v1;
    Vector3d v2;
    bool consecutive = true;
    size_t ind = 2;

    // TODO: Handle possible collinear consequtive vectors
    size_t i0 = f->indices[0];
    size_t i1 = f->indices[1];
    while (ind < f->sh.nbIndex && consecutive) {
        size_t i2 = f->indices[ind++];

        v1 = vertices3[i1] - vertices3[i0]; // v1 = P0P1
        v2 = vertices3[i2] - vertices3[i1]; // v2 = P1P2
        f->sh.N = CrossProduct(v1, v2);              // Cross product
        consecutive = (f->sh.N.Length() < 1e-3);
    }
    f->sh.N = f->sh.N.Normalized();                  // Normalize

    // Calculate Axis Aligned Bounding Box
    f->sh.bb.min = Vector3d(1e100, 1e100, 1e100);
    f->sh.bb.max = Vector3d(-1e100, -1e100, -1e100);

    for (const auto& i : f->indices) {
        const Vector3d& p = vertices3[i];
        f->sh.bb.min.x = std::min(f->sh.bb.min.x,p.x);
        f->sh.bb.min.y = std::min(f->sh.bb.min.y, p.y);
        f->sh.bb.min.z = std::min(f->sh.bb.min.z, p.z);
        f->sh.bb.max.x = std::max(f->sh.bb.max.x, p.x);
        f->sh.bb.max.y = std::max(f->sh.bb.max.y, p.y);
        f->sh.bb.max.z = std::max(f->sh.bb.max.z, p.z);
    }

    // Facet center (AxisAlignedBoundingBox center)
    f->sh.center = 0.5 * (f->sh.bb.max + f->sh.bb.min);

    // Plane equation
    //double A = f->sh.N.x;
    //double B = f->sh.N.y;
    //double C = f->sh.N.z;
    //double D = -Dot(f->sh.N, p0);

    Vector3d p1 = vertices3[f->indices[1]];

    Vector3d U, V;

    U = (p1 - p0).Normalized(); //First side

    // Construct a normal vector V:
    V = CrossProduct(f->sh.N, U); // |U|=1 and |N|=1 => |V|=1

    // u,v vertices (we start with p0 at 0,0)
    f->vertices2[0].u = 0.0;
    f->vertices2[0].v = 0.0;
    Vector2d BBmin; BBmin.u = 0.0; BBmin.v = 0.0;
    Vector2d BBmax; BBmax.u = 0.0; BBmax.v = 0.0;

    for (size_t j = 1; j < f->sh.nbIndex; j++) {
        Vector3d p = vertices3[f->indices[j]];
        Vector3d v = p - p0;
        f->vertices2[j].u = Dot(U, v);  // Project p on U along the V direction
        f->vertices2[j].v = Dot(V, v);  // Project p on V along the U direction

        // Bounds
        BBmax.u  = std::max(BBmax.u , f->vertices2[j].u);
        BBmax.v = std::max(BBmax.v, f->vertices2[j].v);
        BBmin.u = std::min(BBmin.u, f->vertices2[j].u);
        BBmin.v = std::min(BBmin.v, f->vertices2[j].v);
    }

    // Calculate facet area (Meister/Gauss formula)
    double area = 0.0;
    for (size_t j = 0; j < f->sh.nbIndex; j++) {
        size_t j_next = Next(j,f->sh.nbIndex);
        area += f->vertices2[j].u*f->vertices2[j_next].v - f->vertices2[j_next].u*f->vertices2[j].v; //Equal to Z-component of vectorial product
    }
    if (area > 0.0) {

    }
    else if (area < 0.0) {
        //This is a case where a concave facet doesn't obey the right-hand rule:
        //it happens when the first rotation (usually around the second index) is the opposite as the general outline rotation

        //Do a flip
        f->sh.N = -1.0 * f->sh.N;
        V = -1.0 * V;
        BBmin.v = BBmax.v = 0.0;
        for (auto& v : f->vertices2) {
            v.v = -1.0 * v.v;
            BBmax.v = std::max(BBmax.v, v.v);
            BBmin.v = std::min(BBmin.v, v.v);
        }
    }

    f->sh.area = std::abs(0.5 * area);

    // Compute the 2D basis (O,U,V)
    double uD = (BBmax.u - BBmin.u);
    double vD = (BBmax.v - BBmin.v);

    // Origin
    f->sh.O = p0 + BBmin.u * U + BBmin.v * V;

    // Rescale U and V vector
    f->sh.nU = U;
    f->sh.U = U * uD;

    f->sh.nV = V;
    f->sh.V = V * vD;

    f->sh.Nuv = CrossProduct(f->sh.U,f->sh.V); //Not normalized normal vector

    // Rescale u,v coordinates
    for (auto& p : f->vertices2) {
        p.u = (p.u - BBmin.u) / uD;
        p.v = (p.v - BBmin.v) / vD;
    }

#if defined(MOLFLOW)
    f->sh.maxSpeed = 4.0 * std::sqrt(2.0*8.31*f->sh.temperature / 0.001 / wp.gasMass);
#endif
}

#include <Vector.h>
#include <random>

struct Triangle {
    Vector3d v0, v1, v2;
    Triangle() : v0(), v1(), v2(){};
    Triangle(const Vector3d& v0, const Vector3d& v1, const Vector3d& v2) : v0(v0), v1(v1), v2(v2) {}
};

// Calculate the area of a triangle
float triangleArea(const Triangle& tri) {
    Vector3d edge1 = tri.v1 - tri.v0;
    Vector3d edge2 = tri.v2 - tri.v0;
    Vector3d crossProduct = CrossProduct(edge1, edge2);
    return 0.5f * sqrtf(crossProduct.x * crossProduct.x + crossProduct.y * crossProduct.y + crossProduct.z * crossProduct.z);
}

// Calculate the centroid of a triangle
Vector3d calculateCentroid(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3) {
    return (v1 + v2 + v3) / 3.0;
}

float randomFloat(float min, float max) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dist(min, max);
    return dist(gen);
}

int SimulationModel::AnalyzeGeom() {
    std::vector<Triangle> triangles;

    // Convert geometry (vertices and facets) to triangles
    for (const auto& facet : facets) {
        for (size_t i = 1; i < facet->indices.size() - 1; ++i) {
            Triangle triangle;
            triangle.v0 = vertices3[facet->indices[0]];
            triangle.v1 = vertices3[facet->indices[i]];
            triangle.v2 = vertices3[facet->indices[i + 1]];
            triangles.push_back(triangle);
        }
    }

    // Compute the AABB and total surface area
    Vector3d minCorner(std::numeric_limits<float>::max());
    Vector3d maxCorner(std::numeric_limits<float>::lowest());
    float totalArea = 0.0f;

    for (const auto& tri : triangles) {
        Vector3d vertices[3] = {tri.v0, tri.v1, tri.v2};

        for (const auto& vertex : vertices) {
            minCorner.x = std::min(minCorner.x, vertex.x);
            minCorner.y = std::min(minCorner.y, vertex.y);
            minCorner.z = std::min(minCorner.z, vertex.z);

            maxCorner.x = std::max(maxCorner.x, vertex.x);
            maxCorner.y = std::max(maxCorner.y, vertex.y);
            maxCorner.z = std::max(maxCorner.z, vertex.z);
        }

        totalArea += triangleArea(tri);
    }

    Vector3d sceneExtent = maxCorner - minCorner;
    float sceneVolume = sceneExtent.x * sceneExtent.y * sceneExtent.z;
    sceneVolume = isinf(sceneVolume) ? 0 : sceneVolume;

    // Compute spatial distribution and depth complexity
    int numSamplePairs = 1000;
    float avgDistance = 0.0f;
    int sceneDepthComplexity = 0;


    Vector3d rayDirection(0, 0, -1); // Change the direction to test different viewpoints

    for (int i = 0; i < numSamplePairs && triangles.size(); ++i) {
        // Randomly select two triangles and compute the distance between their centroids
        int idx1 = static_cast<int>(randomFloat(0, triangles.size()));
        int idx2 = static_cast<int>(randomFloat(0, triangles.size()));

        Vector3d centroid1 = (triangles[idx1].v0 + triangles[idx1].v1 + triangles[idx1].v2) / 3.0f;
        Vector3d centroid2 = (triangles[idx2].v0 + triangles[idx2].v1 + triangles[idx2].v2) / 3.0f;

        avgDistance += Distance(centroid1, centroid2) / numSamplePairs;

        // Compute scene depth complexity
        Vector3d startPos = centroid1 + rayDirection * sceneExtent.z;
        for (const auto& tri : triangles) {
            Vector3d edge1 = tri.v1 - tri.v0;
            Vector3d edge2 = tri.v2 - tri.v0;
            Vector3d normal = CrossProduct(edge1, edge2);

            // Check if the ray intersects the triangle using the dot product with the normal
            if (normal.x * rayDirection.x + normal.y * rayDirection.y + normal.z * rayDirection.z < 0) {
                ++sceneDepthComplexity;
            }
        }
    }

    // Calculate the average polygon vertex count
    size_t totalVertexCount = 0;
    for (const auto& facet : facets) {
        totalVertexCount += facet->indices.size();
    }
    double averageVertexCount = static_cast<double>(totalVertexCount) / facets.size();
    std::vector<size_t> vertexCounts;
    vertexCounts.reserve(facets.size());
    for (const auto& facet : facets) {
        vertexCounts.push_back(facet->indices.size());
    }
    std::sort(vertexCounts.begin(), vertexCounts.end());

    size_t medianVertexCount;
    if (vertexCounts.size() % 2 == 0) {
        // Even number of polygons
        medianVertexCount = (vertexCounts[vertexCounts.size() / 2 - 1] + vertexCounts[vertexCounts.size() / 2]) / 2;
    } else {
        // Odd number of polygons
        medianVertexCount = vertexCounts[vertexCounts.size() / 2];
    }

    // Calculate the median polygon area
    std::vector<double> polygonAreas;
    polygonAreas.reserve(triangles.size());
    for (const Triangle& triangle : triangles) {
        polygonAreas.push_back(triangleArea(triangle));
    }
    std::sort(polygonAreas.begin(), polygonAreas.end());
    double medianPolygonArea = !polygonAreas.empty() ? polygonAreas[polygonAreas.size() / 2] : 0;

    // Calculate the bounding sphere radius
    Vector3d bboxCenter = (minCorner + maxCorner) / 2.0;
    double boundingSphereRadius = 0.0;
    for (const Vector3d& vertex : vertices3) {
        double distance = Distance(vertex,bboxCenter);
        boundingSphereRadius = std::max(boundingSphereRadius, distance);
    }


    // Calculate the average and maximum triangle edge length
    double totalEdgeLength = 0.0;
    double maxEdgeLength = 0.0;
    for (const Triangle& triangle : triangles) {
        double edge1Length = (triangle.v0 - triangle.v1).Length();
        double edge2Length = (triangle.v1 - triangle.v2).Length();
        double edge3Length = (triangle.v2 - triangle.v1).Length();

        totalEdgeLength += edge1Length + edge2Length + edge3Length;

        maxEdgeLength = std::max(maxEdgeLength, std::max(edge1Length, std::max(edge2Length, edge3Length)));
    }
    double averageEdgeLength = totalEdgeLength / (3.0 * triangles.size());

    float aspectRatio = std::max(sceneExtent.x, std::max(sceneExtent.y, sceneExtent.z)) / std::min(sceneExtent.x, std::min(sceneExtent.y, sceneExtent.z));

    int numTriangles = static_cast<int>(triangles.size());
    float avgArea = numTriangles > 0 ? (totalArea /  numTriangles) : 0;
    float primitiveDensity = sceneVolume > 0 ? (numTriangles / sceneVolume) : 0;

    // Calculate the average polygon vertex count
    size_t nbRealTri = 0;
    size_t nbRealQuad = 0;
    size_t nbRealPol = 0;
    for (const auto& facet : facets) {
        nbRealTri += facet->indices.size() == 3 ? 1 : 0;
        nbRealQuad += facet->indices.size() == 4 ? 1 : 0;
        nbRealPol += facet->indices.size() > 4 ? 1 : 0;
    }
    // Output results for scene
    std::cout << "Number of facets: " << facets.size() << std::endl;
    std::cout << "Number of triangles (real): " << nbRealTri << std::endl;
    std::cout << "Number of real rectangles: " << nbRealQuad << std::endl;
    std::cout << "Number of real n.polygons: " << nbRealPol << std::endl;

    std::cout << "Average polygon vertex count: " << averageVertexCount << std::endl;
    std::cout << "Median polygon vertex count: " << medianVertexCount << std::endl;
    std::cout << "Median polygon area: " << medianPolygonArea << std::endl;
    std::cout << "Bounding sphere radius: " << boundingSphereRadius << std::endl;

    std::cout << "Number of triangles (tri-mesh): " << numTriangles << std::endl;
    std::cout << "Avg Edge of triangles: " << averageEdgeLength << std::endl;
    std::cout << "Max Edge of triangles: " << maxEdgeLength << std::endl;

    //std::cout << "AABB min corner: (" << minCorner.x << ", " << minCorner.y << ", " << minCorner.z << ")" << std::endl;
    //std::cout << "AABB max corner: (" << maxCorner.x << ", " << maxCorner.y << ", " << maxCorner.z << ")" << std::endl;
    std::cout << "Average triangle area: " << avgArea << std::endl;
    std::cout << "Scene volume: " << sceneVolume << std::endl;
    std::cout << "Average distance between triangle centroids: " << avgDistance << std::endl;
    std::cout << "Aspect ratio: " << aspectRatio << std::endl;
    std::cout << "Scene depth complexity: " << sceneDepthComplexity << std::endl;
    std::cout << "Primitive density: " << primitiveDensity << std::endl;

    return 0;
}

/*#include <iostream>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include <random>*/
