//
// Created by pbahr on 08/02/2021.
//

#ifndef MOLFLOW_PROJ_FACETDATA_H
#define MOLFLOW_PROJ_FACETDATA_H

#include <vector>
#include <RayTracing/Primitive.h>
#include "Vector.h"
#include "Buffer_shared.h"
#include <Random.h>

class Surface {
public:
    virtual ~Surface() = default;
    virtual bool IsHardHit(const Ray &r) { return true; };
};

class TransparentSurface : public Surface{
public:
    bool IsHardHit(const Ray &r) override {
        return false;
    };
};

class AlphaSurface : public Surface{
    double opacity{1.0};
public:
    AlphaSurface(double opacity) : opacity(opacity){};
    bool IsHardHit(const Ray &r) override {
        return (r.rng->rnd() < opacity);
    };
};

struct GeomPrimitive : public RTPrimitive {
    GeomPrimitive() : RTPrimitive(), sh(0){ surf = nullptr; };
    explicit GeomPrimitive(size_t nbIndex) : RTPrimitive(), sh(nbIndex){ surf = nullptr; };
    ~GeomPrimitive() override{
        if (surf) {
            //delete surf;
            // don' t delete, origin is an unreferenced shared ptr
            surf = nullptr;
        }
    }
    size_t globalId; //Global index (to identify when superstructures are present)
    FacetProperties sh;
    std::vector<size_t>      indices;          // Indices (Reference to geometry vertex)
    std::vector<Vector2d> vertices2;        // Vertices (2D plane space, UV coordinates)
    Surface* surf;
};

struct Facet : public GeomPrimitive {
    Facet() : GeomPrimitive(){ };
    Facet(size_t nbIndex) : GeomPrimitive(nbIndex) { };
    virtual ~Facet(){ }

    //size_t iSCount{0};

    void ComputeBB() override { bb = sh.bb;};
    bool Intersect(Ray &r) override;

};

struct TriangleFacet : public GeomPrimitive {
    TriangleFacet() : GeomPrimitive(3){ };
    virtual ~TriangleFacet() { }
    std::vector<Vector3d>* vertices3;        // Vertices (2D plane space, UV coordinates)
    Vector2d texCoord[3];

    //size_t iSCount{0};

    void ComputeBB() override { bb = sh.bb;};
    bool Intersect(Ray &r) override;
};

#endif //MOLFLOW_PROJ_FACETDATA_H
