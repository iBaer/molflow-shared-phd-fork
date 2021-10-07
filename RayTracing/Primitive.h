//
// Created by pascal on 4/22/21.
//

#ifndef MOLFLOW_PROJ_RTPRIMITIVE_H
#define MOLFLOW_PROJ_RTPRIMITIVE_H

#include "BoundingBox.h"
#include <atomic>

class Ray;

struct IntersectCount{
    IntersectCount() : nbChecks(0), nbIntersects(0){};
    IntersectCount(const IntersectCount& cpy) : nbChecks(cpy.nbChecks.load()), nbIntersects(cpy.nbIntersects.load()){
        nbPrim = cpy.nbPrim;
        level = cpy.level;
    };
    IntersectCount(IntersectCount&& mvr)  noexcept : nbChecks(mvr.nbChecks.load()), nbIntersects(mvr.nbIntersects.load()){
        nbPrim = mvr.nbPrim;
        level = mvr.level;
    };
    IntersectCount& operator=(const IntersectCount& cpy){
        nbChecks = cpy.nbChecks.load();
        nbIntersects = cpy.nbIntersects.load();
        nbPrim = cpy.nbPrim;
        level = cpy.level;

        return *this;
    }

    IntersectCount& operator=(IntersectCount&& cpy) noexcept {
        nbChecks = cpy.nbChecks.load();
        nbIntersects = cpy.nbIntersects.load();
        nbPrim = cpy.nbPrim;
        level = cpy.level;

        return *this;
    }

    void Reset() {
        nbChecks = 0;
        nbIntersects = 0;
    }


    std::atomic<size_t> nbChecks{0};
    std::atomic<size_t> nbIntersects{0};
    size_t nbPrim{0};
    size_t level{0};
};

// General ray tracing interface
class RTPrimitive {
public:
    virtual ~RTPrimitive() = default;
    virtual void ComputeBB() = 0;
    virtual bool Intersect(Ray &r) = 0;
    virtual bool IntersectStat(RayStat &r) = 0;
    AxisAlignedBoundingBox bb;
};
#endif //MOLFLOW_PROJ_RTPRIMITIVE_H
