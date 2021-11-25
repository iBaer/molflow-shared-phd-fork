//
// Created by pascal on 4/22/21.
//

#ifndef MOLFLOW_PROJ_RAY_H
#define MOLFLOW_PROJ_RAY_H

#include "Vector.h"
#include "RTHelper.h"

//struct SubProcessFacetTempVar;
class MersenneTwister;

struct HitChain {
    size_t hitId;
    SubProcessFacetTempVar *hit;
    HitChain *next;
};

struct HitLink {
    HitLink() : hitId(9999999999), hit(SubProcessFacetTempVar()) {};
    HitLink(size_t id, SubProcessFacetTempVar h) : hitId(id), hit(h) {};

    // Move constructor called on resize, prevent from deleting SubProcessFacetTempVar
    HitLink(const HitLink &rhs) = default;

    HitLink(HitLink &&rhs) noexcept:
            hitId(rhs.hitId),
            hit(rhs.hit) {};

    HitLink &operator=(const HitLink &src) {
        hitId = src.hitId;
        hit = src.hit;
        return *this;
    };

    HitLink &operator=(HitLink &&src) {
        hitId = src.hitId;
        hit = src.hit;
        return *this;
    };

    ~HitLink();

    size_t hitId;
    SubProcessFacetTempVar hit;
};

struct Payload {
};
constexpr double inf_d = 1.0e99;

class Ray {
public:
    Ray() : tMax(inf_d), time(0.f), structure(-1), lastIntersected(-1), hitChain(nullptr), rng(nullptr), pay(nullptr) {}

    Ray(const Vector3d &o, const Vector3d &d, Payload *payload, double tMax = inf_d,
        double time = 0.f, int structure = -1)
            : origin(o), direction(d), tMax(tMax), time(time), structure(structure), hitChain(nullptr), rng(nullptr),
              pay(payload) {}

    ~Ray() {
        if (pay) {
            delete pay;
            pay=nullptr;
        }
    }

    Vector3d operator()(double t) const { return origin + direction * t; }

    Vector3d origin;
    Vector3d direction;

    // To keep track of shortest intersection
    double tMax;

    double time; // Only for td simulations in Molflow
    int lastIntersected; //

    int structure; //

    //const Medium *medium;
    Payload *pay;

    HitChain *hitChain;
    std::vector<HitLink> hits;
    HitLink hardHit;
    std::vector<HitLink> transparentHits; // TODO: Remove Extra debug structure
    MersenneTwister *rng;
};

struct RayStatistics {
    size_t traversalSteps{0};
    size_t nbIntersectionTests{0};
    size_t nbBoxIntersectionTests{0};
    size_t nbDownTest{0};
    size_t nbUpTest{0};

    RayStatistics& operator+=(const RayStatistics& src) {
        this->traversalSteps += src.traversalSteps;
        this->nbIntersectionTests += src.nbIntersectionTests;
        this->nbBoxIntersectionTests += src.nbBoxIntersectionTests;
        this->nbDownTest += src.nbDownTest;
        this->nbUpTest += src.nbUpTest;

        return *this;
    }

    inline void Reset(){
        traversalSteps = 0;
        nbIntersectionTests = 0;
        nbBoxIntersectionTests = 0;
        nbDownTest = 0;
        nbUpTest = 0;
    }
};

class RayStat : public Ray {
public:
    RayStat() : Ray() {}

    explicit RayStat(const Ray &src) : Ray(src) {}

    RayStat(const Vector3d &o, const Vector3d &d, Payload *payload = nullptr, double tMax = inf_d,
        double time = 0.f, int structure = 0)
            : Ray(o, d, payload, tMax, time, structure) {}

    virtual ~RayStat() = default;
    RayStat &operator=(const RayStat &src) noexcept {
        origin = src.origin;
        direction = src.direction;
        lastIntersected = src.lastIntersected;
        structure = src.structure;
        tMax = src.tMax;
        time = src.time;

        return *this;
    }

    void ResetStats(){
        stats.Reset();
    }

    //Statistics
    RayStatistics stats;
    std::vector<size_t> traversedNodes;
};

#endif //MOLFLOW_PROJ_RAY_H
