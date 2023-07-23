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

/*
 * KDTree code based on pbrt-v3
 * Copyright (c) 1998-2015, Matt Pharr, Greg Humphreys, and Wenzel Jakob.
 */


#include <BoundingBox.h>
#include <cmath>
#include <memory>
#include <set>
#include <cassert>
#include <random>
#include <Helper/ConsoleLogger.h>

#include "KDTree.h"
#include "Ray.h"

#include <Helper/Chronometer.h>

#include <fstream>

namespace Profiling {
    CummulativeBenchmark intersectStatsKD{};
    CummulativeBenchmark intersectStatsKDRope{};
    CummulativeBenchmark intersectStatsKDRopeRestart{};
}

static Chronometer time_trav;
static Chronometer time_rope;
static Chronometer time_roperestart;

namespace STATS {
    //STAT_MEMORY_COUNTER("Memory/BVH tree", treeBytes);
    //STAT_RATIO("BVH/Primitives per leaf node", totalPrimitives, totalLeafNodes);

    static int totalPrimitives = 0;
    static int totalLeafNodes = 0;
    static int interiorNodes = 0;
    static int leafNodes = 0;
}

inline int Log2Int(uint32_t v) {
#if defined(_MSC_VER)
    unsigned long lz = 0;
    if (_BitScanReverse(&lz, v)) return lz;
    return 0;
#else
    return 31 - __builtin_clz(v);
#endif
}

inline int Log2Int(int32_t v) { return Log2Int((uint32_t) v); }

inline int Log2Int(uint64_t v) {
#if defined(_MSC_VER)
    unsigned long lz = 0;
#if defined(_WIN64)
    _BitScanReverse64(&lz, v);
#else
    if  (_BitScanReverse(&lz, v >> 32))
        lz += 32;
    else
        _BitScanReverse(&lz, v & 0xffffffff);
#endif // _WIN64
    return lz;
#else  // _MSC_VER
    return 63 - __builtin_clzll(v);
#endif
}

inline int Log2Int(int64_t v) { return Log2Int((uint64_t) v); }

namespace STATS_KD {
    //STAT_MEMORY_COUNTER("Memory/BVH tree", treeBytes);
    //STAT_RATIO("BVH/Primitives per leaf node", totalPrimitives, totalLeafNodes);

    static int totalPrimitives = 0;
    static int totalLeafNodes = 0;
    static int interiorNodes = 0;
    static int leafNodes = 0;
    static int leafBadRefine = 0;
    static int leafHigherCost = 0;
    static int splitRay = 0;
    static int splitSAH = 0;

    void _reset() {
        totalPrimitives = 0;
        totalLeafNodes = 0;
        interiorNodes = 0;
        leafNodes = 0;
        leafBadRefine = 0;
        leafHigherCost = 0;
        splitRay = 0;
        splitSAH = 0;
    }
}

void KdAccelNode::InitInterior(int axis, int ac, double s) {
    split = s;
    flags = axis;
    aboveChild |= (ac << 2);
    STATS_KD::interiorNodes++;
}

enum class EdgeType {
    Start, End
};

struct BoundEdge {
    // BoundEdge Public Methods
    BoundEdge() = default;
    BoundEdge(double t, int primNum, bool starting) : t(t), primNum(primNum) {
        type = starting ? EdgeType::Start : EdgeType::End;
    }

    double t{};
    int primNum{};
    EdgeType type;
};

void KdTreeAccel::PrintTreeInfo() {
    std::string splitName;
    switch (splitMethod) {
        case (SplitMethod::SAH):
            splitName = "SAH";
            break;
        case (SplitMethod::ProbSplit):
            splitName = "Prob";
            break;
        case (SplitMethod::ProbHybrid):
            splitName = "ProbHybrid";
            break;
        case (SplitMethod::TestSplit):
            splitName = "Test";
            break;
        case (SplitMethod::HybridSplit):
            splitName = "Hybrid";
            break;
        case (SplitMethod::HybridBin):
            splitName = "HybridBin";
            break;
        default:
            splitName = "Unknown split";
    }
    Log::console_msg_master(4, "--- KD STATS ---\n");
    Log::console_msg_master(4, "- Split: {} -\n", splitName.c_str());
    Log::console_msg_master(4, " Total Primitives: {}\n", STATS_KD::totalPrimitives);
    Log::console_msg_master(4, " Total Leaf Nodes: {}\n", STATS_KD::totalLeafNodes);
    Log::console_msg_master(4, " Interior Nodes:   {}\n", STATS_KD::interiorNodes);
    Log::console_msg_master(4, " Leaf Nodes:       {}\n", STATS_KD::leafNodes);
    Log::console_msg_master(4, " Leaf 4 BadRefine: {} [{}]\n", STATS_KD::leafBadRefine,
           (double) STATS_KD::leafBadRefine / STATS_KD::leafNodes);
    Log::console_msg_master(4, " Leaf 4 BadCost:   {} [{}]\n", STATS_KD::leafHigherCost,
           (double) STATS_KD::leafHigherCost / STATS_KD::leafNodes);
    Log::console_msg_master(4, " Split SAH:        {} [{}]\n", STATS_KD::splitSAH,
           (double) STATS_KD::splitSAH / (STATS_KD::splitRay + STATS_KD::splitSAH));
    Log::console_msg_master(4, " Split Ray:        {} [{}]\n", STATS_KD::splitRay,
           (double) STATS_KD::splitRay / (STATS_KD::splitRay + STATS_KD::splitSAH));
}

void KdTreeAccel::ComputeBB() {
    bb = AxisAlignedBoundingBox();
    if (nodes) {
        for (const auto &prim: primitives) {
            bb = AxisAlignedBoundingBox::Union(bb, prim->bb);
        }
    }
}

// KdTreeAccel Method Definitions
KdTreeAccel::KdTreeAccel(SplitMethod splitMethod, std::vector<std::shared_ptr<Primitive>> p,
                         const std::vector<double> &probabilities, const std::vector<TestRay> &battery,
                         double hybridWeight,
                         int isectCost, int traversalCost, double emptyBonus, int maxPrims, int maxDepth)
        : splitMethod(splitMethod),
          isectCost(isectCost),
          traversalCost(traversalCost),
          maxPrims(maxPrims),
          emptyBonus(emptyBonus),
          primitives(std::move(p)) {

    nodes = nullptr;

    if (primitives.empty())
        return;
    if ((splitMethod == KdTreeAccel::SplitMethod::ProbSplit || splitMethod == KdTreeAccel::SplitMethod::ProbHybrid) && probabilities.empty())
        return;
    STATS_KD::_reset();

    // Build kd-tree for accelerator
    nextFreeNode = nAllocedNodes = 0;
    if (maxDepth <= 0)
        maxDepth = std::round(8.0f + 1.3f * Log2Int(int64_t(primitives.size())));

    // Compute bounds for kd-tree construction
    std::vector<AxisAlignedBoundingBox> primBounds;
    primBounds.reserve(primitives.size());
    for (const std::shared_ptr<Primitive> &prim: primitives) {
        AxisAlignedBoundingBox b = prim->sh.bb;
        bounds = AxisAlignedBoundingBox::Union(bounds, b);
        primBounds.push_back(b);
    }

    // Allocate working memory for kd-tree construction
    std::unique_ptr<BoundEdge[]> edges[3];
    for (int i = 0; i < 3; ++i)
        edges[i] = std::make_unique<BoundEdge[]>(2 * primitives.size());
    std::unique_ptr<int[]> prims0(new int[primitives.size()]);
    std::unique_ptr<int[]> prims1(new int[(maxDepth + 1) * primitives.size()]);

    // Initialize _primNums_ for kd-tree construction
    std::unique_ptr<int[]> primNums(new int[primitives.size()]);
    for (size_t i = 0; i < primitives.size(); ++i) {
        primNums[i] = i;
    }

    std::vector<IntersectCount>(0).swap(ints);

    bool withBattery = (splitMethod == KdTreeAccel::SplitMethod::TestSplit)
                       || (splitMethod == KdTreeAccel::SplitMethod::HybridSplit)
                       || (splitMethod == KdTreeAccel::SplitMethod::HybridBin);
    bool withSmartBattery = withBattery;

    if (true && withBattery) {
        Log::console_msg_master(4, "--- KD with HitBattery ---\n");
        Log::console_msg_master(4, " Battery size: {}\n", battery.size());
        // if battery is too large, only use a random sample

        auto *batter_ptr = const_cast<std::vector<TestRay> *>(&battery);
        if (battery.size() > HITCACHELIMIT) {
            batter_ptr = new std::vector<TestRay>;
            std::sample(battery.begin(), battery.end(), std::back_inserter(*batter_ptr),
                        HITCACHELIMIT, std::mt19937{std::random_device{}()});
            Log::console_msg_master(4, " Sample size: {}\n", batter_ptr->size());
        }
        std::vector<TestRayLoc> indices;
        indices.reserve(batter_ptr->size());
        for (int i = 0; i < batter_ptr->size(); ++i)
            indices.emplace_back(i, 0, 1.0e99);

        double costLimit = 1.0e99;
        // Start recursive construction of kd-tree
        buildTree(0, bounds, primBounds, primNums.get(), primitives.size(),
                  maxDepth, edges, prims0.get(), prims1.get(), 0, hybridWeight, *batter_ptr, indices, -1, 0.0, 1.0e99,
                  probabilities, costLimit);

        if (battery.size() > HITCACHELIMIT)
            delete batter_ptr;
    } else if (!battery.empty() && withSmartBattery) {
        std::unique_ptr<std::vector<RaySegment>> rays(new std::vector<RaySegment>(battery.size()));
        for (int r = 0; r < battery.size(); ++r) {
            rays->at(r).ray = &(battery[r]);
        }

        RayBoundary *rb_stack = new RayBoundary[battery.size() * 2 * 3];
        //RayBoundary* rb_arr[3];
        RayBoundary **rb_arr_ptr = new RayBoundary *[3];
        auto rb_arr = rb_arr_ptr;
        rb_arr[0] = &rb_stack[battery.size() * 2 * 0];
        rb_arr[1] = &rb_stack[battery.size() * 2 * 1];
        rb_arr[2] = &rb_stack[battery.size() * 2 * 2];
        auto *rb_ptr = rb_stack;

        for (int i = 0; i < battery.size(); ++i) {
            auto *ray = &rays->at(i);
            // left bound
            assert(i + battery.size() * 5 < battery.size() * 2 * 3);
            rb_ptr[i].rs = ray;
            rb_ptr[i].flag = 0;
            rb_ptr[i + battery.size() * 2].rs = ray;
            rb_ptr[i + battery.size() * 2].flag = 0;
            rb_ptr[i + battery.size() * 4].rs = ray;
            rb_ptr[i + battery.size() * 4].flag = 0;
            // right bound
            rb_ptr[i + battery.size()].rs = ray;
            rb_ptr[i + battery.size()].flag = 2;
            rb_ptr[i + battery.size() * 3].rs = ray;
            rb_ptr[i + battery.size() * 3].flag = 2;
            rb_ptr[i + battery.size() * 5].rs = ray;
            rb_ptr[i + battery.size() * 5].flag = 2;

            assert(rb_ptr[i + battery.size() * 5].flag == rb_arr[2][i + battery.size()].flag);
        }

        BoundaryStack boundaryStack{0, static_cast<int>(battery.size() - 1), 0, static_cast<int>(battery.size() - 1)};
        double costLimit = 1.0e99;
        // Start recursive construction of kd-tree
        buildTreeRDH(0, bounds, primBounds, primNums.get(), primitives.size(),
                     maxDepth, edges, prims0.get(), prims1.get(), 0, hybridWeight, rays, rb_arr_ptr,
                     -1, 0.0, 1.0e99, probabilities, costLimit, boundaryStack);

        delete[] rb_stack;
    } else {
        buildTree(0, bounds, primBounds, primNums.get(), primitives.size(),
                  maxDepth, edges, prims0.get(), prims1.get(), 0, probabilities, hybridWeight, -1);
    }

    for (int nodeNum = 0; nodeNum < STATS_KD::interiorNodes + STATS_KD::leafNodes; ++nodeNum) {
        ints[nodeNum].level = maxDepth - ints[nodeNum].level;
    }

    Log::console_msg_master(4, "--- KD STATS ---\n");
    Log::console_msg_master(4, " Total Primitives: {}\n", STATS_KD::totalPrimitives);
    Log::console_msg_master(4, " Total Leaf Nodes: {}\n", STATS_KD::totalLeafNodes);
    Log::console_msg_master(4, " Interior Nodes:   {}\n", STATS_KD::interiorNodes);
    Log::console_msg_master(4, " Leaf Nodes:       {}\n", STATS_KD::leafNodes);

    hasRopes = false;
    addParent(&nodes[0], nullptr);
    PrintTreeInfo();
}

void KdTreeAccel::AddRopes(bool optimized) {
    if(!nodes)
        return;
    ResetStats();

    hasRopes = true;
    KdAccelNode* ropes_loc[6] = { NULL }; // 2 ropes for each dimension
    attachRopes(&this->nodes[0], ropes_loc, optimized);
}

void KdTreeAccel::RemoveRopes(){
    hasRopes = false;
    ResetStats();

    for(int n = 0; n < nAllocedNodes; ++n){
        for (auto & rope : nodes[n].ropes) {
            rope = nullptr;
        }
    }
}

void KdAccelNode::InitLeaf(int *primNums, int np,
                           std::vector<int> *primitiveIndices) {
    flags = 3;
    nPrims |= (np << 2);
    // Store primitive ids for leaf node
    if (np == 0)
        onePrimitive = 0;
    else if (np == 1)
        onePrimitive = primNums[0];
    else {
        primitiveIndicesOffset = primitiveIndices->size();
        for (int i = 0; i < np; ++i) primitiveIndices->push_back(primNums[i]);
    }
    STATS_KD::leafNodes++;
    STATS_KD::totalPrimitives += np;
}

KdTreeAccel::~KdTreeAccel() {
    if (nodes) {
        delete[] nodes;
        nodes = nullptr;
    }
}

std::tuple<double, int, int> KdTreeAccel::SplitSAH(int axis, const AxisAlignedBoundingBox &nodeBounds,
                                                   const std::vector<AxisAlignedBoundingBox> &allPrimBounds,
                                                   int *primNums,
                                                   int nPrimitives, const std::unique_ptr<BoundEdge[]> edges[3]) {

    // Choose split axis position for interior node
    int bestAxis = -1, bestOffset = -1;
    double bestCost = 1.0e99;
    double oldCost = isectCost * double(nPrimitives);
    double totalSA = nodeBounds.SurfaceArea();
    double invTotalSA = 1.0 / totalSA;
    Vector3d d = nodeBounds.max - nodeBounds.min;

    // Initialize edges for _axis_
    for (int i = 0; i < nPrimitives; ++i) {
        int pn = primNums[i];
        const AxisAlignedBoundingBox &apBounds = allPrimBounds[pn];
        edges[axis][2 * i] = BoundEdge(apBounds.min[axis], pn, true);
        edges[axis][2 * i + 1] = BoundEdge(apBounds.max[axis], pn, false);
    }

    // Sort _edges_ for _axis_
    std::sort(&edges[axis][0], &edges[axis][2 * nPrimitives],
              [](const BoundEdge &e0, const BoundEdge &e1) -> bool {
                  if (e0.t == e1.t)
                      return (int) e0.type < (int) e1.type;
                  else
                      return e0.t < e1.t;
              });

    // Compute cost of all splits for _axis_ to find best
    int nBelow = 0, nAbove = nPrimitives;

    for (int i = 0; i < 2 * nPrimitives; ++i) {
        if (edges[axis][i].type == EdgeType::End) {
            --nAbove;
        }
        double edgeT = edges[axis][i].t;
        if (edgeT > nodeBounds.min[axis] && edgeT < nodeBounds.max[axis]) {
            // Compute cost for split at _i_th edge

            // Compute child surface areas for split at _edgeT_
            int otherAxis0 = (axis + 1) % 3, otherAxis1 = (axis + 2) % 3;
            double belowSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                                  (edgeT - nodeBounds.min[axis]) *
                                  (d[otherAxis0] + d[otherAxis1]));
            double aboveSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                                  (nodeBounds.max[axis] - edgeT) *
                                  (d[otherAxis0] + d[otherAxis1]));
            double pBelow = belowSA * invTotalSA;
            double pAbove = aboveSA * invTotalSA;
            double eb = (nAbove == 0 || nBelow == 0) ? emptyBonus : 0;

            double cost =
                    traversalCost +
                    isectCost * (1 - eb) * (pBelow * nBelow + pAbove * nAbove);

            // Update best split if this is lowest cost so far
            if (cost < bestCost) {
                bestCost = cost;
                bestAxis = axis;
                bestOffset = i;
            }
        }
        if (edges[axis][i].type == EdgeType::Start) {
            ++nBelow;
        }
    }

    //CHECK(nBelow == nPrimitives && nAbove == 0);
    return {bestCost, bestAxis, bestOffset};
};

std::tuple<double, int, int> KdTreeAccel::SplitProb(int axis, const AxisAlignedBoundingBox &nodeBounds,
                                                    const std::vector<AxisAlignedBoundingBox> &allPrimBounds,
                                                    int *primNums, int nPrimitives,
                                                    const std::unique_ptr<BoundEdge[]> edges[3],
                                                    const std::vector<double> &probabilities, const double weight) const {

    // Choose split axis position for interior node
    int bestAxis = -1, bestOffset = -1;
    double bestCost = 1.0e99;
    double oldCost = isectCost * double(nPrimitives);
    double totalSA = nodeBounds.SurfaceArea();
    double invTotalSA = 1.0 / totalSA;
    Vector3d d = nodeBounds.max - nodeBounds.min;

    // Initialize edges for _axis_
    for (int i = 0; i < nPrimitives; ++i) {
        int pn = primNums[i];
        const AxisAlignedBoundingBox &apBounds = allPrimBounds[pn];
        edges[axis][2 * i] = BoundEdge(apBounds.min[axis], pn, true);
        edges[axis][2 * i + 1] = BoundEdge(apBounds.max[axis], pn, false);
    }

    // Sort _edges_ for _axis_
    std::sort(&edges[axis][0], &edges[axis][2 * nPrimitives],
              [](const BoundEdge &e0, const BoundEdge &e1) -> bool {
                  if (e0.t == e1.t)
                      return (int) e0.type < (int) e1.type;
                  else
                      return e0.t < e1.t;
              });

    // Compute cost of all splits for _axis_ to find best
    int nBelow = 0, nAbove = nPrimitives;
    // w/ probability
    double probBelow = 0.0;
    double probAbove = 0.0;

    // add probability for all edges going in
    std::set<int> edgeIn;
    if (!probabilities.empty()) {
        for (int i = 0; i < 2 * nPrimitives; ++i) {
            auto inside = edgeIn.find(edges[axis][i].primNum);
            if (inside != edgeIn.end()) {
                //probAbove += edges[axis][i].type == EdgeType::End ? probabilities[edges[axis][i].primNum] : 0.0;
                probAbove += probabilities[edges[axis][i].primNum];
                edgeIn.erase(inside);
            } else {
                edgeIn.insert(edges[axis][i].primNum);
            }
        }
    }

    for (int i = 0; i < 2 * nPrimitives; ++i) {
        if (edges[axis][i].type == EdgeType::End) {
            --nAbove;
            if (!probabilities.empty()) {
                auto inside = edgeIn.find(edges[axis][i].primNum);
                if (inside == edgeIn.end())
                    probAbove -= probabilities[edges[axis][i].primNum];
            }
        }
        double edgeT = edges[axis][i].t;
        if (edgeT > nodeBounds.min[axis] && edgeT < nodeBounds.max[axis]) {
            // Compute cost for split at _i_th edge

            // Compute child surface areas for split at _edgeT_
            int otherAxis0 = (axis + 1) % 3, otherAxis1 = (axis + 2) % 3;
            double belowSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                                  (edgeT - nodeBounds.min[axis]) *
                                  (d[otherAxis0] + d[otherAxis1]));
            double aboveSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                                  (nodeBounds.max[axis] - edgeT) *
                                  (d[otherAxis0] + d[otherAxis1]));
            double pBelow = belowSA * invTotalSA;
            double pAbove = aboveSA * invTotalSA;
            double eb = (nAbove == 0 || nBelow == 0) ? emptyBonus : 0;

            const double weightSA = 0.0;
            double pca = (probAbove + probBelow > 0) ? weightSA + (2.0 - weightSA) * probAbove / (probAbove + probBelow)
                                                     : 1.0;
            double pcb = (probAbove + probBelow > 0) ? weightSA + (2.0 - weightSA) * probBelow / (probAbove + probBelow)
                                                     : 1.0;

            double cost =
                    traversalCost +
                    isectCost * (1 - eb) * (pcb * pBelow * nBelow + pca * pAbove * nAbove);

            // Update best split if this is lowest cost so far
            if (cost < bestCost) {
                bestCost = cost;
                bestAxis = axis;
                bestOffset = i;
            }
        }
        if (edges[axis][i].type == EdgeType::Start) {
            ++nBelow;
            if (!probabilities.empty()) {
                auto inside = edgeIn.find(edges[axis][i].primNum);
                if (inside == edgeIn.end())
                    probBelow += probabilities[edges[axis][i].primNum];
            }
        }
    }

    //CHECK(nBelow == nPrimitives && nAbove == 0);
    return {bestCost, bestAxis, bestOffset};
};

std::tuple<double, int, int> KdTreeAccel::SplitTest(int axis, const AxisAlignedBoundingBox &nodeBounds,
                                                    const std::vector<AxisAlignedBoundingBox> &allPrimBounds,
                                                    int *primNums, int nPrimitives,
                                                    const std::unique_ptr<BoundEdge[]> edges[3],
                                                    const std::unique_ptr<std::vector<RaySegment>> &battery,
                                                    RayBoundary **rb_stack,
                                                    const std::vector<double> &primChance, double tMax,
                                                    const BoundaryStack &boundary_stack) {

    // Choose split axis position for interior node
    int bestAxis = -1, bestOffset = -1;
    double bestCost = 1.0e99;
    double oldCost = isectCost * double(nPrimitives);

    // Initialize edges for _axis_
    for (int i = 0; i < nPrimitives; ++i) {
        int pn = primNums[i];
        const AxisAlignedBoundingBox &apBounds = allPrimBounds[pn];
        edges[axis][2 * i] = BoundEdge(apBounds.min[axis], pn, true);
        edges[axis][2 * i + 1] = BoundEdge(apBounds.max[axis], pn, false);
    }

    // Sort _edges_ for _axis_
    std::sort(&edges[axis][0], &edges[axis][2 * nPrimitives],
              [](const BoundEdge &e0, const BoundEdge &e1) -> bool {
                  if (e0.t == e1.t)
                      return (int) e0.type < (int) e1.type;
                  else
                      return e0.t < e1.t;
              });

    // Compute cost of all splits for _axis_ to find best
    int nBelow = 0, nAbove = nPrimitives;

    for (int i = 0; i < 2 * nPrimitives; ++i) {
        if (edges[axis][i].type == EdgeType::End) {
            --nAbove;
        }
        double edgeT = edges[axis][i].t;

        auto boundary = rb_stack[axis];
        BoundaryStack bound{boundary_stack}; // dummy todo: replace with parameter
        size_t countAbove = bound.leftMax - bound.leftMin;
        size_t countBelow = bound.rightMax - bound.rightMin;
        int splitIndex = -1;
        for (int b = bound.leftMin; b < bound.rightMax; ++b) {
            auto tVal = (boundary[b].flag == 0) ? boundary[b].rs->tMin : boundary[b].rs->tMax;
            if (edgeT < tVal) {
                continue;
            } else {
                // found
                splitIndex = b;
            }
        }

        countAbove = bound.leftMin + std::max(splitIndex, 0);
        countBelow = bound.leftMax - countAbove;

        if (edgeT > nodeBounds.min[axis] && edgeT < nodeBounds.max[axis]) {
            // Compute cost for split at _i_th edge
            double eb = (nAbove == 0 || nBelow == 0) ? emptyBonus : 0;

            {

                // count boundaries

                double cost = traversalCost
                              + isectCost * (1 - eb) *
                                ((double) countBelow / (countAbove + countBelow) * nBelow +
                                 ((double) countAbove / (countAbove + countBelow) * nAbove));
                if (cost < bestCost) {
                    bestCost = cost;
                    bestAxis = axis;
                    bestOffset = i;
                }
            }
        }
        if (edges[axis][i].type == EdgeType::Start) {
            ++nBelow;
        }
    }
    //CHECK(nBelow == nPrimitives && nAbove == 0);

    return {bestCost, bestAxis, bestOffset};
};

std::tuple<double, int, int> KdTreeAccel::SplitTest(int axis, const AxisAlignedBoundingBox &nodeBounds,
                                                    const std::vector<AxisAlignedBoundingBox> &allPrimBounds,
                                                    int *primNums, int nPrimitives,
                                                    const std::unique_ptr<BoundEdge[]> edges[3],
                                                    const std::vector<TestRay> &battery,
                                                    const std::vector<TestRayLoc> &local_battery,
                                                    const std::vector<double> &primChance, double tMax) {

    // Choose split axis position for interior node
    int bestAxis = -1, bestOffset = -1;
    double bestCost = 1.0e99;
    double oldCost = isectCost * double(nPrimitives);

    // Initialize edges for _axis_
    for (int i = 0; i < nPrimitives; ++i) {
        int pn = primNums[i];
        const AxisAlignedBoundingBox &apBounds = allPrimBounds[pn];
        edges[axis][2 * i] = BoundEdge(apBounds.min[axis], pn, true);
        edges[axis][2 * i + 1] = BoundEdge(apBounds.max[axis], pn, false);
    }

    double probBelow = 0.0;
    double probAbove = 0.0;
    for (int i = 0; i < 2 * nPrimitives; ++i) {
        /*if (edges[axis][i].type == EdgeType::Start)
            probAbove = edges[axis][i].primNum;*/
        if (edges[axis][i].type == EdgeType::End)
            probAbove += primChance[edges[axis][i].primNum];
    }
    // normalise to [0..1]
    double probMax = probAbove;

    // Sort _edges_ for _axis_
    std::sort(&edges[axis][0], &edges[axis][2 * nPrimitives],
              [](const BoundEdge &e0, const BoundEdge &e1) -> bool {
                  if (e0.t == e1.t)
                      return (int) e0.type < (int) e1.type;
                  else
                      return e0.t < e1.t;
              });

    // Compute cost of all splits for _axis_ to find best
    int nBelow = 0, nAbove = nPrimitives;
    double inv_denum = 1.0 / (double) local_battery.size();

    for (int i = 0; i < 2 * nPrimitives; ++i) {
        if (edges[axis][i].type == EdgeType::End) {
            --nAbove;
            probAbove -= primChance[edges[axis][i].primNum];
        }
        double edgeT = edges[axis][i].t;
        if (edgeT > nodeBounds.min[axis] && edgeT < nodeBounds.max[axis]) {
            // Compute cost for split at _i_th edge
            double eb = (nAbove == 0 || nBelow == 0) ? emptyBonus : 0;

            {//if(local_battery.size() >= HITCACHEMIN) {
                int hitCountB = 0;
                int hitCountA = 0;
                int hitCountBoth = 0;

                auto ray = Ray();
                double locTMax = tMax;
#pragma omp parallel default(none) firstprivate(ray,locTMax) shared(local_battery)
                {
                    int hitCountB_local = 0;
                    int hitCountA_local = 0;
                    int hitCountBoth_local = 0;
#pragma omp for
                    for (int sample_id = 0; sample_id < local_battery.size(); sample_id++) {
                        auto &ind = local_battery[sample_id];
                        ray.origin = battery[ind.index].pos;
                        ray.direction = battery[ind.index].dir;
                        Vector3d invDir(1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z);
                        int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};
                        /*hitcount1 += b1.IntersectBox(ray, invDir, dirIsNeg) ? 1 : 0;
                        hitcount0 += b0.IntersectBox(ray, invDir, dirIsNeg) ? 1 : 0;*/

                        // Compute parametric distance along ray to split plane
                        double tSplit = edgeT;
                        //int axis = node->SplitAxis();
                        double tPlane = (tSplit - ray.origin[axis]) * invDir[axis];

                        // Get node children pointers for ray
                        const KdAccelNode *firstChild, *secondChild;
                        int belowFirst =
                                (ray.origin[axis] < tSplit) ||
                                (ray.origin[axis] == tSplit && ray.direction[axis] <= 0);

                        // Advance to next child node, possibly enqueue other child
                        if (tPlane > ind.tMax || tPlane <= 0)
                            hitCountA_local++;// test first child
                        else if (tPlane < ind.tMin)
                            hitCountB_local++;// test second child
                        else {
                            hitCountBoth_local++;
                        }
                        // combine local
#pragma omp critical
                        {
                            hitCountA += hitCountA_local;
                            hitCountB += hitCountB_local;
                            hitCountBoth += hitCountBoth_local;

                        }
                    }
                }
                double pAbove = (double) (hitCountA + hitCountBoth) * inv_denum;
                double pBelow = (double) (hitCountB + hitCountBoth) * inv_denum;
                double pBoth = (double) (hitCountBoth) * inv_denum;

                /*if(hitCountB-hitCountA < 0) {
                    pAbove = (double) (std::abs(hitCountB - hitCountA)) / (double) (hitCountA + hitCountB);
                    pBelow = 1.0 - pAbove;
                }
                else {
                    pBelow = (double) (std::abs(hitCountA - hitCountB)) / (double) (hitCountA + hitCountB);
                    pAbove = 1.0 - pBelow;
                }*/

                double cost =
                        traversalCost
                        + isectCost * (1 - eb) *
                          (/*pBelow **/ probBelow / probMax * nBelow + /*pAbove **/ probAbove / probMax * nAbove);
                if (cost < bestCost) {
                    bestCost = cost;
                    bestAxis = axis;
                    bestOffset = i;
                }
            }
        }
        if (edges[axis][i].type == EdgeType::Start) {
            ++nBelow;
            probBelow += primChance[edges[axis][i].primNum];
        }
    }
    //CHECK(nBelow == nPrimitives && nAbove == 0);

    return {bestCost, bestAxis, bestOffset};
};

std::tuple<double, int, int> KdTreeAccel::SplitHybridProb(int axis, const AxisAlignedBoundingBox &nodeBounds,
                                                          const std::vector<AxisAlignedBoundingBox> &allPrimBounds,
                                                          int *primNums, int nPrimitives,
                                                          const std::unique_ptr<BoundEdge[]> edges[3],
                                                          const std::vector<double> &probabilities,
                                                          const double hybridWeight) const {

    // Choose split axis position for interior node
    int bestAxis = -1, bestOffset = -1;
    double bestCost = 1.0e99;
    double oldCost = isectCost * double(nPrimitives);
    double totalSA = nodeBounds.SurfaceArea();
    double invTotalSA = 1.0 / totalSA;
    Vector3d d = nodeBounds.max - nodeBounds.min;

    // Initialize edges for _axis_
    for (int i = 0; i < nPrimitives; ++i) {
        int pn = primNums[i];
        const AxisAlignedBoundingBox &apBounds = allPrimBounds[pn];
        edges[axis][2 * i] = BoundEdge(apBounds.min[axis], pn, true);
        edges[axis][2 * i + 1] = BoundEdge(apBounds.max[axis], pn, false);
    }

    // Sort _edges_ for _axis_
    std::sort(&edges[axis][0], &edges[axis][2 * nPrimitives],
              [](const BoundEdge &e0, const BoundEdge &e1) -> bool {
                  if (e0.t == e1.t)
                      return (int) e0.type < (int) e1.type;
                  else
                      return e0.t < e1.t;
              });

    double probBelow = 0.0;
    double probAbove = 0.0;

    // add probability for all edges going in
    std::set<int> edgeIn;
    if (!probabilities.empty()) {
        for (int i = 0; i < 2 * nPrimitives; ++i) {
            probAbove += probabilities[edges[axis][i].primNum];

            /*auto inside = edgeIn.find(edges[axis][i].primNum);
            if (inside != edgeIn.end()) {
                //probAbove += edges[axis][i].type == EdgeType::End ? probabilities[edges[axis][i].primNum] : 0.0;
                probAbove += probabilities[edges[axis][i].primNum];
                edgeIn.erase(inside);
            } else {
                edgeIn.insert(edges[axis][i].primNum);
            }*/
        }
    }

    // normalise to [0..1]
    double probMax = probAbove;
    //probAbove = 1.0;

    // Compute cost of all splits for _axis_ to find best
    int nBelow = 0, nAbove = nPrimitives;

#ifdef WRITE_KD_COSTFILE
    std::ofstream cost_writer;
    if(nPrimitives > 800 && nPrimitives < 850){
        cost_writer.open(fmt::format("costfile_{}.txt",nextFreeNode-1));
    }
#endif
    for (int i = 0; i < 2 * nPrimitives; ++i) {
        if (edges[axis][i].type == EdgeType::End) {
            --nAbove;
            if (!probabilities.empty()) {
                //auto inside = edgeIn.find(edges[axis][i].primNum);
                //if (inside == edgeIn.end())
                    probAbove -= probabilities[edges[axis][i].primNum];
            }
        }
        double edgeT = edges[axis][i].t;
        if (edgeT > nodeBounds.min[axis] && edgeT < nodeBounds.max[axis]) {
            // Compute cost for split at _i_th edge

            // Compute child surface areas for split at _edgeT_
            int otherAxis0 = (axis + 1) % 3, otherAxis1 = (axis + 2) % 3;
            double belowSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                                  (edgeT - nodeBounds.min[axis]) *
                                  (d[otherAxis0] + d[otherAxis1]));
            double aboveSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                                  (nodeBounds.max[axis] - edgeT) *
                                  (d[otherAxis0] + d[otherAxis1]));
            double pBelow = belowSA * invTotalSA;
            double pAbove = aboveSA * invTotalSA;
            double eb = (nAbove == 0 || nBelow == 0) ? emptyBonus : 0;

            double cost_sah =
                    traversalCost +
                    isectCost * (1 - eb) * (pBelow * nBelow + pAbove * nAbove);

            double un_cost = traversalCost
                             + isectCost * (1 - eb) *
                               (/*pBelow **/ probBelow / (probMax) * nBelow + /*pAbove **/ probAbove / (probMax) * nAbove);

            {
                double cost = 0.0;
                double alpha = 0.0;
                double beta = 1.0 - alpha;

                if(hybridWeight >= 0.0){
                    alpha = hybridWeight;
                    beta = 1.0 - alpha;
                    cost = traversalCost
                           + isectCost * (1 - eb) *
                             (/*pBelow **/ probBelow / (probAbove+probBelow) * nBelow + /*pAbove **/ probAbove / (probAbove+probBelow) * nAbove);

                }
                else {
                    if (1.0f / (float) this->primitives.size() <= probMax) {
                        alpha = 1.0;
                        beta = 1.0 - alpha;
                        cost = traversalCost
                               + isectCost * (1 - eb) *
                                 (/*pBelow **/ probBelow / (probAbove+probBelow) * nBelow +
                                               /*pAbove **/ probAbove / (probAbove+probBelow) * nAbove);
                    }
                }

                //TODO: Maybe use nB Hits used for prob calculation as smoothing factor (beta * smooth)
                //double weight = alpha * (1.0 - (1.0 / (1.0 + beta)));
                double linCost = alpha * cost + beta * cost_sah;

                if (linCost < bestCost) {
                    bestCost = linCost;
                    bestAxis = axis;
                    bestOffset = i;
                }

#ifdef WRITE_KD_COSTFILE
                if(cost_writer.is_open()){
                    cost_writer << edges[axis][i].t << " " << probabilities[edges[axis][i].primNum] << " " << cost_sah << " " << cost << " " << un_cost << "\n";
                }
#endif
            }
        }
        if (edges[axis][i].type == EdgeType::Start) {
            ++nBelow;
            if (!probabilities.empty()) {
                //auto inside = edgeIn.find(edges[axis][i].primNum);
                //if (inside == edgeIn.end())
                    probBelow += probabilities[edges[axis][i].primNum];
            }
        }
    }

#ifdef WRITE_KD_COSTFILE
    if(cost_writer.is_open()){
        cost_writer.close();
    }
#endif
    //CHECK(nBelow == nPrimitives && nAbove == 0);

    return {bestCost, bestAxis, bestOffset};
};

std::tuple<double, int, int> KdTreeAccel::SplitHybrid(int axis, const AxisAlignedBoundingBox &nodeBounds,
                                                      const std::vector<AxisAlignedBoundingBox> &allPrimBounds,
                                                      int *primNums, int nPrimitives,
                                                      const std::unique_ptr<BoundEdge[]> edges[3],
                                                      const std::vector<TestRay> &battery,
                                                      const std::vector<TestRayLoc> &local_battery,
                                                      const std::vector<double> &primChance, double tMax,
                                                      double hybridWeight) const {

    // Choose split axis position for interior node
    int bestAxis = -1, bestOffset = -1;
    double bestCost = 1.0e99;
    double oldCost = isectCost * double(nPrimitives);
    double totalSA = nodeBounds.SurfaceArea();
    double invTotalSA = 1.0 / totalSA;
    Vector3d d = nodeBounds.max - nodeBounds.min;

    // Initialize edges for _axis_
    for (int i = 0; i < nPrimitives; ++i) {
        int pn = primNums[i];
        const AxisAlignedBoundingBox &apBounds = allPrimBounds[pn];
        edges[axis][2 * i] = BoundEdge(apBounds.min[axis], pn, true);
        edges[axis][2 * i + 1] = BoundEdge(apBounds.max[axis], pn, false);
    }

    double probBelow = 0.0;
    double probAbove = 0.0;
    for (int i = 0; i < 2 * nPrimitives; ++i) {
        /*if (edges[axis][i].type == EdgeType::Start)
            probAbove = edges[axis][i].primNum;*/
        if (edges[axis][i].type == EdgeType::End)
            probAbove += primChance[edges[axis][i].primNum];
    }
    // normalise to [0..1]
    double probMax = probAbove;
    //probAbove = 1.0;

    // Sort _edges_ for _axis_
    std::sort(&edges[axis][0], &edges[axis][2 * nPrimitives],
              [](const BoundEdge &e0, const BoundEdge &e1) -> bool {
                  if (e0.t == e1.t)
                      return (int) e0.type < (int) e1.type;
                  else
                      return e0.t < e1.t;
              });

    // Compute cost of all splits for _axis_ to find best
    int nBelow = 0, nAbove = nPrimitives;
    double inv_denum = 1.0 / local_battery.size();

#ifdef WRITE_KD_COSTFILE
    std::ofstream cost_writer;
    if(local_battery.size() > 1000 && local_battery.size() < 1200){
        cost_writer.open(fmt::format("costfile_rdh_{}.txt",nextFreeNode-1));
    }
#endif

    int nbEdgeInsideNode = 0;
    for (int i = 0; i < 2 * nPrimitives; ++i) {
        if (edges[axis][i].type == EdgeType::End) {
            --nAbove;
            probAbove -= primChance[edges[axis][i].primNum];
        }
        double edgeT = edges[axis][i].t;
        if (edgeT > nodeBounds.min[axis] && edgeT < nodeBounds.max[axis]) {
            nbEdgeInsideNode++;
            // Compute cost for split at _i_th edge

            // Compute child surface areas for split at _edgeT_
            int otherAxis0 = (axis + 1) % 3, otherAxis1 = (axis + 2) % 3;
            double belowSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                                  (edgeT - nodeBounds.min[axis]) *
                                  (d[otherAxis0] + d[otherAxis1]));
            double aboveSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                                  (nodeBounds.max[axis] - edgeT) *
                                  (d[otherAxis0] + d[otherAxis1]));
            double pBelow = belowSA * invTotalSA;
            double pAbove = aboveSA * invTotalSA;
            double eb = (nAbove == 0 || nBelow == 0) ? emptyBonus : 0;

            double cost_sah =
                    traversalCost +
                    isectCost * (1 - eb) * (pBelow * nBelow + pAbove * nAbove);

            {//if(local_battery.size() >= HITCACHEMIN) {
                int hitCountB = 0;
                int hitCountA = 0;
                int hitCountBoth = 0;
                double pBoth = 0.0;

                auto ray = Ray();
                double locTMax = tMax;
#pragma omp parallel default(none) firstprivate(ray) shared(battery, axis, hitCountA, hitCountB, hitCountBoth, local_battery, edges, bestAxis, bestOffset, edgeT)
                {
                    int hitCountB_loc = 0;
                    int hitCountA_loc = 0;
                    int hitCountBoth_loc = 0;

#pragma omp for
                    for (int sample_id = 0; sample_id < local_battery.size(); sample_id++) {
                        auto &ind = local_battery[sample_id];
                        ray.origin = battery[ind.index].pos;
                        ray.direction = battery[ind.index].dir;
                        Vector3d invDir(1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z);
                        int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};
                        /*hitcount1 += b1.IntersectBox(ray, invDir, dirIsNeg) ? 1 : 0;
                        hitcount0 += b0.IntersectBox(ray, invDir, dirIsNeg) ? 1 : 0;*/

                        // Compute parametric distance along ray to split plane
                        double tSplit = edgeT;
                        //int axis = node->SplitAxis();
                        double tPlane = (tSplit - ray.origin[axis]) * invDir[axis];

                        // Get node children pointers for ray
                        const KdAccelNode *firstChild, *secondChild;
                        int belowFirst =
                                (ray.origin[axis] < tSplit) ||
                                (ray.origin[axis] == tSplit && ray.direction[axis] <= 0);

                        // Advance to next child node, possibly enqueue other child
                        if (tPlane > ind.tMax || tPlane <= 0)
                            hitCountA_loc++;// test first child
                        else if (tPlane < ind.tMin)
                            hitCountB_loc++;// test second child
                        else {
                            hitCountBoth_loc++;
                        }
                    }

                    // parallel reduce
#pragma omp critical
                    {
                        hitCountA += hitCountA_loc;
                        hitCountB += hitCountB_loc;
                        hitCountBoth += hitCountBoth_loc;
                    }
                }
                pAbove = (double) (hitCountA + hitCountBoth) * inv_denum;
                pBelow = (double) (hitCountB + hitCountBoth) * inv_denum;
                pBoth = (double) (hitCountBoth) * inv_denum;

                /*if(hitCountB-hitCountA < 0) {
                    pAbove = (double) (std::abs(hitCountB - hitCountA)) / (double) (hitCountA + hitCountB);
                    pBelow = 1.0 - pAbove;
                }
                else {
                    pBelow = (double) (std::abs(hitCountA - hitCountB)) / (double) (hitCountA + hitCountB);
                    pAbove = 1.0 - pBelow;
                }*/

                double cost = traversalCost
                              + isectCost * (1.0 - eb) *
                                (pBelow * nBelow + pAbove * nAbove);

                /*double alpha = 0.1;
                double beta = 0.9;
                double weight = alpha * (1.0 - (1.0 / (1.0 + beta * (double) local_battery.size())));
                */
                double weight = hybridWeight;
                double linCost = weight * cost + (1.0 - weight) * cost_sah;

                if (linCost < bestCost) {
                    bestCost = linCost;
                    bestAxis = axis;
                    bestOffset = i;
                }

#ifdef WRITE_KD_COSTFILE
                if(cost_writer.is_open()){
                    cost_writer << edges[axis][i].t << " " << cost_sah << " " << cost << " " << linCost << "\n";
                }
#endif
            }
        }
        if (edges[axis][i].type == EdgeType::Start) {
            ++nBelow;
            probBelow += primChance[edges[axis][i].primNum];
        }
    }

    /*if(nbEdgeInsideNode == 0){
        fmt::print("Warning, no edges were inside node, could result into malformed node\n");
    }*/

#ifdef WRITE_KD_COSTFILE
    if(cost_writer.is_open()){
        cost_writer.close();
    }
#endif
    //CHECK(nBelow == nPrimitives && nAbove == 0);

    return {bestCost, bestAxis, bestOffset};
};

struct RDHBucketInfo {
    int count = 0;
    int indPos = 0;
    double pos = 0.0;
    double probBelow = 0.0;
    double probAbove = 0.0;

    int hitCountB = 0;
    int hitCountA = 0;
    int hitCountBoth = 0;
    double pBoth = 0.0;

    int nAbove = 0;
    int nBelow = 0;
    int countl = 0;
    int countr = 0;
    int ncloserl = 0;
    int ncloserr = 0;
    AxisAlignedBoundingBox bounds0;
    AxisAlignedBoundingBox bounds1;
};

int IntersectBins(Ray& ray, std::vector<RDHBucketInfo>& buckets, int axis){

    for (auto &bin : buckets) {
        double distl = 1e99;
        double distr = 1e99;

        if(bin.count==0)
            continue; // is count the correct attribute to skip ?
        double t0;
        double t1;
        if(bin.bounds0.IntersectP(ray, &t0, &t1, axis)){
            bin.countl++;
            distl = t0;
        }

        if(bin.bounds1.IntersectP(ray, &t0, &t1, axis)){
            bin.countr++;
            distr = t0;
        }

        if(distl < distr){
            bin.ncloserl++;
        }
        else{
            bin.ncloserr++;
        }
    }
    return 0;
}

std::tuple<double, int, int, bool, bool> KdTreeAccel::SplitHybridBin(int axis, const AxisAlignedBoundingBox &nodeBounds,
                                                                     const std::vector<AxisAlignedBoundingBox> &allPrimBounds,
                                                                     int *primNums, int nPrimitives,
                                                                     const std::unique_ptr<BoundEdge[]> edges[3],
                                                                     const std::vector<TestRay> &battery,
                                                                     const std::vector<TestRayLoc> &local_battery,
                                                                     const std::vector<double> &primChance, double tMax) const {

    // Choose split axis position for interior node
    int bestAxis = -1, bestOffset = -1;
    double bestCost = 1.0e99;
    double oldCost = isectCost * double(nPrimitives);
    double totalSA = nodeBounds.SurfaceArea();
    double invTotalSA = 1.0 / totalSA;
    Vector3d d = nodeBounds.max - nodeBounds.min;

    // Initialize edges for _axis_
    for (int i = 0; i < nPrimitives; ++i) {
        int pn = primNums[i];
        const AxisAlignedBoundingBox &apBounds = allPrimBounds[pn];
        edges[axis][2 * i] = BoundEdge(apBounds.min[axis], pn, true);
        edges[axis][2 * i + 1] = BoundEdge(apBounds.max[axis], pn, false);
    }

    double probBelow = 0.0;
    double probAbove = 0.0;
    for (int i = 0; i < 2 * nPrimitives; ++i) {
        /*if (edges[axis][i].type == EdgeType::Start)
            probAbove = edges[axis][i].primNum;*/
        if (edges[axis][i].type == EdgeType::End)
            probAbove += primChance[edges[axis][i].primNum];
    }
    // normalise to [0..1]
    double probMax = probAbove;
    //probAbove = 1.0;

    // Sort _edges_ for _axis_
    std::sort(&edges[axis][0], &edges[axis][2 * nPrimitives],
              [](const BoundEdge &e0, const BoundEdge &e1) -> bool {
                  if (e0.t == e1.t)
                      return (int) e0.type < (int) e1.type;
                  else
                      return e0.t < e1.t;
              });

    // subdivide into n bins
    constexpr int nBuckets = 128;
    std::vector<RDHBucketInfo> buckets(nBuckets);

    //<<Initialize BucketInfo for SAH partition buckets>>

    auto posMin = bounds.min;
    auto posMax = bounds.max;
    for (auto &b: buckets) {
        b.nAbove = nPrimitives;
        b.probAbove = probAbove;
        b.bounds0 = nodeBounds;
        b.bounds1 = nodeBounds;
    }

    int start = 0;
    int end = 2 * nPrimitives;
    for (int i = 0; i < 2 * nPrimitives; ++i) {
        auto edge = edges[axis][i].t;
        start = i;
        if (edge > nodeBounds.min[axis] && edge < nodeBounds.max[axis]) {
            break;
        }
    }
    for (int i = 2 * nPrimitives - 1; i >= start; --i) {
        auto edge = edges[axis][i].t;
        end = i;
        if (edge > nodeBounds.min[axis] && edge < nodeBounds.max[axis]) {
            break;
        }
    }

    double stride = (nodeBounds.max[axis] - nodeBounds.min[axis]) / nBuckets;
    //for (int i = 0; i < 2 * nPrimitives; ++i) {
    for (int i = start; i <= end; ++i) {
        auto& edge = edges[axis][i];
        double edgeT = edge.t;
        if (edgeT <= nodeBounds.min[axis] || edgeT >= nodeBounds.max[axis]) {
            continue;
        }
        int b = (double) (edgeT - nodeBounds.min[axis]) / stride;
    /*double stride = (double)(end - start) / nBuckets;
    if(stride == 0)
        stride = 1;
    if(stride < nBuckets)
        buckets.resize((int)stride);

    for (int i = start; i <= end; ++i) {
        auto& edge = edges[axis][i];
        double edgeT = edge.t;
        int b = (double) (i - start) / stride;*/

        //int b = nBuckets * edges[axis][i].t;
        if (b == nBuckets) b = nBuckets - 1;
        assert(b >= 0);
        assert(b < nBuckets);
        for(int bb = b; bb < buckets.size(); ++bb){
            buckets[bb].indPos = i;
            if (edge.type == EdgeType::End) {
                buckets[bb].count--;
                buckets[bb].probAbove -= primChance[edge.primNum];
                --buckets[bb].nAbove;
            }
            if (edge.type == EdgeType::Start) {
                buckets[bb].count++;
                buckets[bb].probBelow += primChance[edge.primNum];
                ++buckets[bb].nBelow;
            }
        }
        buckets[b].pos = edgeT;
        //buckets[b].bounds = AxisAlignedBoundingBox::Union(buckets[b].bounds, edges[axis][i]);
    }

    for (auto &bucket: buckets) {
        bucket.bounds0.max[bestAxis] = bucket.bounds1.min[bestAxis] = bucket.pos;
    }

    // ray intersection computation
    // ray sampling
    bool skipL = false;
    bool skipR = false;
    std::vector<TestRayLoc> battery_below, battery_above;
    if(!skipL || !skipR)
    {
        auto ray = Ray();


        double locTMax = tMax;
#pragma omp parallel default(none) firstprivate(ray, locTMax) shared(local_battery, buckets, battery, axis, edges, bestAxis, bestOffset)
        {
            std::vector<RDHBucketInfo> buckets_local = buckets;
#pragma omp for
            for (int sample_id = 0; sample_id < local_battery.size(); sample_id++) {
                auto &ind = local_battery[sample_id];
                ray.origin = battery[ind.index].pos;
                ray.direction = battery[ind.index].dir;
                /*Vector3d invDir(1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z);
                int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};
                */
                IntersectBins(ray, buckets_local, axis);
            }

            // combine local
#pragma omp critical
            {
                for(int bi = 0; bi < buckets.size(); bi++){
                    buckets[bi].countl += buckets_local[bi].countl;
                    buckets[bi].countr += buckets_local[bi].countr;
                    buckets[bi].ncloserl += buckets_local[bi].ncloserl;
                    buckets[bi].ncloserr += buckets_local[bi].ncloserr;
                }
            }
        }
    }

    // Compute cost of all splits for _axis_ to find best
    //int nBelow = 0, nAbove = nPrimitives;
    double inv_denum = 1.0 / (double)local_battery.size();

    //<<Compute costs for splitting after each bucket>>
    std::vector<double> cost(buckets.size());
    std::vector<double> cost_min(buckets.size());
    std::vector<double> cost_mid(buckets.size());
    std::vector<double> cost_sah(buckets.size());
    std::vector<double> cost_rdh(buckets.size());

    for (int i = 0; i < buckets.size(); ++i) {
        cost[i] = 9e99;
        cost_min[i] = 9e99;
        cost_mid[i] = 9e99;
        cost_sah[i] = 9e99;
        cost_rdh[i] = 9e99;
    }

    /*bool*/ skipL = false;
    /*bool*/ skipR = false;

    for (int i = 0; i < buckets.size(); ++i) {
        if(buckets[i].count == 0)
            continue;
        double edgeT = buckets[i].pos;
        if (edgeT > nodeBounds.min[axis] && edgeT < nodeBounds.max[axis]) {
            // Compute cost for split at _i_th edge
            double eb = (buckets[i].nAbove == 0 || buckets[i].nBelow == 0) ? emptyBonus : 0;

            // SAH - Compute child surface areas for split at _edgeT_
            {
                int otherAxis0 = (axis + 1) % 3, otherAxis1 = (axis + 2) % 3;
                double belowSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                                      (edgeT - nodeBounds.min[axis]) *
                                      (d[otherAxis0] + d[otherAxis1]));
                double aboveSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                                      (nodeBounds.max[axis] - edgeT) *
                                      (d[otherAxis0] + d[otherAxis1]));
                double pBelow = belowSA * invTotalSA;
                double pAbove = aboveSA * invTotalSA;
                double eb = (buckets[i].nAbove == 0 || buckets[i].nBelow == 0) ? emptyBonus : 0;

                cost_sah[i] =
                        traversalCost +
                        isectCost * (1 - eb) * (pBelow * buckets[i].nBelow + pAbove * buckets[i].nAbove);
            }

            // RDH
            {
                double pAbove = (double) (buckets[i].hitCountA + buckets[i].hitCountBoth) * inv_denum;
                double pBelow = (double) (buckets[i].hitCountB + buckets[i].hitCountBoth) * inv_denum;
                buckets[i].pBoth = (double) (buckets[i].hitCountBoth) * inv_denum;


                cost_rdh[i] = traversalCost
                              + isectCost * (1 - eb) *
                                (/*pBelow **/ probBelow / probMax * buckets[i].nBelow + /*pAbove **/ probAbove / probMax * buckets[i].nAbove);
            }

            // Hybrid
            {
                double alpha = 0.1;
                double beta = 0.9;
                double weight = alpha * (1.0 - (1.0 / (1.0 + beta * (double) local_battery.size())));
                cost[i] = weight * cost_rdh[i] + (1.0 - weight) * cost_sah[i];
            }

            // DAC
            {
                double pAbove = (double) (buckets[i].countl) * inv_denum;
                double pBelow = (double) (buckets[i].countr) * inv_denum;

                cost[i] = traversalCost
                              + isectCost * (1 - eb) * (pBelow * buckets[i].nBelow + pAbove * buckets[i].nAbove);
            }


            if (cost[i] < bestCost) {
                bestCost = cost[i];
                bestAxis = axis;
                bestOffset = buckets[i].indPos;
                skipL = (double) (buckets[i].countl) * inv_denum > 0.5;
                skipR = (double) (buckets[i].countr) * inv_denum > 0.5;
            }
        }
        // -----
    }

    //<<Find bucket to split at that minimizes PROB metric>>
    /*bestCost = cost[0];
    for (int i = 1; i < nBuckets - 1; ++i) {
        if(buckets[i].count == 0)
            continue;
        if (cost[i] < bestCost) {
            bestCost = cost[i];
            bestAxis = axis;
            bestOffset = buckets[i].indPos;
        }
    }*/

    //CHECK(nBelow == nPrimitives && nAbove == 0);

    return {bestCost, bestAxis, bestOffset, skipL, skipR};
};

void KdTreeAccel::buildTree(int nodeNum, const AxisAlignedBoundingBox &nodeBounds,
                            const std::vector<AxisAlignedBoundingBox> &allPrimBounds, int *primNums, int nPrimitives,
                            int depth, const std::unique_ptr<BoundEdge[]> edges[3], int *prims0, int *prims1,
                            int badRefines, const std::vector<double> &probabilities, const double hybridWeight,
                            int prevSplitAxis) {
    //CHECK_EQ(nodeNum, nextFreeNode);
    assert(nodeNum == nextFreeNode);
    // Get next free node from _nodes_ array
    if (nextFreeNode == nAllocedNodes) {
        int nNewAllocNodes = std::max(2 * nAllocedNodes, 512);
        KdAccelNode *n = new KdAccelNode[nNewAllocNodes];
        if (nAllocedNodes > 0) {
            memcpy(n, nodes, nAllocedNodes * sizeof(KdAccelNode));
            delete[] nodes;
        }
        nodes = n;
        nAllocedNodes = nNewAllocNodes;
    }
    ++nextFreeNode;

    nodes[nodeNum].nodeId = nodeNum;
    nodes[nodeNum].parent = nullptr;
    for(auto& rope : nodes[nodeNum].ropes)
        rope = nullptr;

    // Excl for rope traversal
    nodes[nodeNum].bbox = nodeBounds;

    ints.emplace_back();
    auto &intstat = ints.back(); // ints[nodeNum]
    intstat.nbPrim = nPrimitives;
    intstat.level = depth;

    // Initialize leaf node if termination criteria met
    if (nPrimitives <= maxPrims || depth == 0) {
        nodes[nodeNum].InitLeaf(primNums, nPrimitives, &primitiveIndices);
        return;
    }

    // Initialize interior node and continue recursion

    // Choose split axis position for interior node
    int bestAxis = -1, bestOffset = -1;
    double bestCost = 1.0e99;
    /*double oldCost = isectCost * double(nPrimitives);
    double totalSA = nodeBounds.SurfaceArea();
    double invTotalSA = 1.0 / totalSA;
    Vector3d d = nodeBounds.max - nodeBounds.min;*/

    // Choose which axis to split along
    int axis = nodeBounds.MaximumExtent();
    int retries = 0;
    retrySplit:

    if (splitMethod == SplitMethod::ProbSplit && !probabilities.empty()) {
        std::tie(bestCost, bestAxis, bestOffset) = SplitProb(axis, nodeBounds, allPrimBounds, primNums, nPrimitives,
                                                             edges, probabilities, 0.5);
    }
    else if (splitMethod == SplitMethod::ProbHybrid && !probabilities.empty()) {
        std::tie(bestCost, bestAxis, bestOffset) = SplitHybridProb(axis, nodeBounds, allPrimBounds, primNums,
                                                                   nPrimitives,
                                                                   edges, probabilities, hybridWeight);
    }
    else /*if (splitMethod == SplitMethod::SAH)*/
        std::tie(bestCost, bestAxis, bestOffset) = SplitSAH(axis, nodeBounds, allPrimBounds, primNums, nPrimitives,
                                                            edges);

    // Create leaf if no good splits were found
    if (bestAxis == -1 && retries < 2) {
        ++retries;
        if (prevSplitAxis >= 0) {
            axis = (axis + 1) % 3;
            if (retries == 1 && axis == prevSplitAxis) {
                axis = (axis + 1) % 3;
            } else if (retries == 2) {
                axis = prevSplitAxis;
            }
        } else
            axis = (axis + 1) % 3;
        goto retrySplit;
    }

    double oldCost = isectCost * double(nPrimitives);
    if (bestCost > oldCost) ++badRefines;
    if ((bestCost > 4 * oldCost && nPrimitives < 16) || bestAxis == -1 ||
        badRefines == 3) {
        if (bestCost > 4 * oldCost) STATS_KD::leafHigherCost++;
        if (badRefines == 3) STATS_KD::leafBadRefine++;
        nodes[nodeNum].InitLeaf(primNums, nPrimitives, &primitiveIndices);
        return;
    }

    // Classify primitives with respect to split
    int n0 = 0, n1 = 0;
    for (int i = 0; i < bestOffset; ++i)
        if (edges[bestAxis][i].type == EdgeType::Start)
            prims0[n0++] = edges[bestAxis][i].primNum;
    for (int i = bestOffset + 1; i < 2 * nPrimitives; ++i)
        if (edges[bestAxis][i].type == EdgeType::End)
            prims1[n1++] = edges[bestAxis][i].primNum;

    STATS_KD::splitSAH++;

    // Recursively initialize children nodes
    double tSplit = edges[bestAxis][bestOffset].t;
    AxisAlignedBoundingBox bounds0 = nodeBounds, bounds1 = nodeBounds;
    bounds0.max[bestAxis] = bounds1.min[bestAxis] = tSplit;

    buildTree(nodeNum + 1, bounds0, allPrimBounds, prims0, n0, depth - 1, edges,
              prims0, prims1 + nPrimitives, badRefines, probabilities, hybridWeight, bestAxis);
    nodes[nodeNum + 1].parent = &nodes[nodeNum];

    int aboveChild = nextFreeNode;
    nodes[nodeNum].InitInterior(bestAxis, aboveChild, tSplit);
    buildTree(aboveChild, bounds1, allPrimBounds, prims1, n1, depth - 1, edges,
              prims0, prims1 + nPrimitives, badRefines, probabilities, hybridWeight, bestAxis);
    nodes[aboveChild].parent = &nodes[nodeNum];
}

// default construction algorithm for RDH
void KdTreeAccel::buildTree(int nodeNum, const AxisAlignedBoundingBox &nodeBounds,
                            const std::vector<AxisAlignedBoundingBox> &allPrimBounds, int *primNums, int nPrimitives,
                            int depth, const std::unique_ptr<BoundEdge[]> edges[3], int *prims0, int *prims1,
                            int badRefines, double hybridWeight, const std::vector<TestRay> &battery,
                            const std::vector<TestRayLoc> &local_battery, int prevSplitAxis, double tMin, double tMax,
                            const std::vector<double> &primChance, double oldCost) {
    //CHECK_EQ(nodeNum, nextFreeNode);
    assert(nodeNum == nextFreeNode);
    // Get next free node from _nodes_ array
    if (nextFreeNode == nAllocedNodes) {
        int nNewAllocNodes = std::max(2 * nAllocedNodes, 512);
        KdAccelNode *n = new KdAccelNode[nNewAllocNodes];
        if (nAllocedNodes > 0) {
            memcpy(n, nodes, nAllocedNodes * sizeof(KdAccelNode));
            delete[] nodes;
        }
        nodes = n;
        nAllocedNodes = nNewAllocNodes;
    }
    ++nextFreeNode;

    nodes[nodeNum].nodeId = nodeNum;
    nodes[nodeNum].parent = nullptr;
    ints.emplace_back();
    auto &intstat = ints.back(); // ints[nodeNum]
    intstat.nbPrim = nPrimitives;
    intstat.level = depth;

    // Initialize leaf node if termination criteria met
    if (nPrimitives <= maxPrims || depth == 0) {
        nodes[nodeNum].InitLeaf(primNums, nPrimitives, &primitiveIndices);
        return;
    }

    // Initialize interior node and continue recursion

    // Choose split axis position for interior node
    int bestAxis = -1, bestOffset = -1;
    double bestCost = 1.0e99;

    //double oldCost = isectCost * double(nPrimitives);
    /*double totalSA = nodeBounds.SurfaceArea();
    double invTotalSA = 1.0 / totalSA;
    Vector3d d = nodeBounds.max - nodeBounds.min;*/

    // Choose which axis to split along
    int axis = nodeBounds.MaximumExtent();
    int retries = 0;

    bool skipL = false;
    bool skipR = false;

    retrySplit:
    if (splitMethod == SplitMethod::TestSplit && !local_battery.empty())
        std::tie(bestCost, bestAxis, bestOffset) = SplitTest(axis, nodeBounds, allPrimBounds, primNums, nPrimitives,
                                                             edges, battery, local_battery, primChance, tMax);
    else if (splitMethod == SplitMethod::HybridSplit && !local_battery.empty())
        std::tie(bestCost, bestAxis, bestOffset) = SplitHybrid(axis, nodeBounds, allPrimBounds, primNums, nPrimitives,
                                                               edges, battery, local_battery, primChance, tMax,
                                                               hybridWeight);
    else if (splitMethod == SplitMethod::HybridBin && !local_battery.empty())
        std::tie(bestCost, bestAxis, bestOffset, skipL, skipR) = SplitHybridBin(axis, nodeBounds, allPrimBounds, primNums,
                                                                  nPrimitives,
                                                                  edges, battery, local_battery, primChance, tMax);
    else
        std::tie(bestCost, bestAxis, bestOffset) = SplitSAH(axis, nodeBounds, allPrimBounds, primNums, nPrimitives,
                                                            edges);

    /*if(local_battery.size() < HITCACHEMIN){
            bestCost = bestCostSAH;
            bestAxis = bestAxisSAH;
            bestOffset = bestOffsetSAH;
    }*/
    // Create leaf if no good splits were found
    if (bestAxis == -1 && retries < 2) {
        ++retries;
        if (prevSplitAxis >= 0) {
            axis = (axis + 1) % 3;
            if (retries == 1 && axis == prevSplitAxis) {
                axis = (axis + 1) % 3;
            } else if (retries == 2) {
                axis = prevSplitAxis;
            }
        } else
            axis = (axis + 1) % 3;
        goto retrySplit;
    }

    if (bestCost > 4.0 * oldCost) {
        ++badRefines;
    }

    if ((bestCost > 4.0 * oldCost && nPrimitives < 16) || bestAxis == -1 ||
        (badRefines >= 3 && nPrimitives < 128) || (badRefines >= 5)) {
        if (bestCost > 4.0 * oldCost) {
            ++STATS_KD::leafHigherCost;
            /*if(nPrimitives > 16)
                std::tie(bestCost, bestAxis, bestOffset) = SplitHybrid(axis, nodeBounds, allPrimBounds, primNums, nPrimitives,
                                                                   edges, battery, local_battery, primChance, tMax,
                                                                   hybridWeight);*/
        }
        if (badRefines >= 3)
            ++STATS_KD::leafBadRefine;
        nodes[nodeNum].InitLeaf(primNums, nPrimitives, &primitiveIndices);
        return;
    }

    // todo: determine traversal order

    std::vector<TestRayLoc> battery_below, battery_above;
    if(!skipL || !skipR)
    {
        auto ray = Ray();
        double locTMax = tMax;
#pragma omp parallel default(none) firstprivate(ray) shared(battery, local_battery, edges, bestAxis, bestOffset, battery_above, battery_below)
        {
            std::vector<TestRayLoc> battery_below_local, battery_above_local;
#pragma omp for
            for (int sample_id = 0; sample_id < local_battery.size(); sample_id++) {
                auto& ind = local_battery[sample_id];
                ray.origin = battery[ind.index].pos;
                ray.direction = battery[ind.index].dir;
                Vector3d invDir(1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z);
                int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};
                /*hitcount1 += b1.IntersectBox(ray, invDir, dirIsNeg) ? 1 : 0;
                hitcount0 += b0.IntersectBox(ray, invDir, dirIsNeg) ? 1 : 0;*/

                {
                    // Compute parametric distance along ray to split plane
                    double tSplit = edges[bestAxis][bestOffset].t;
                    //int axis = node->SplitAxis();
                    double tPlane = (tSplit - ray.origin[bestAxis]) * invDir[bestAxis];

                    // Get node children pointers for ray
                    int belowFirst =
                            (ray.origin[bestAxis] < tSplit) ||
                            (ray.origin[bestAxis] == tSplit && ray.direction[bestAxis] <= 0);

                    // Advance to next child node, possibly enqueue other child
                    {
                        if (tPlane > ind.tMax || tPlane <= 0) {
                            battery_above_local.emplace_back(ind.index, ind.tMin, ind.tMax);
                        } else if (tPlane < ind.tMin) {
                            battery_below_local.emplace_back(ind.index, ind.tMin, ind.tMax);
                        } else {
                            battery_above_local.emplace_back(ind.index, ind.tMin, tPlane);
                            battery_below_local.emplace_back(ind.index, tPlane, ind.tMax);
                        }
                    }
                }
            }

            // combine local
#pragma omp critical
            {
                battery_above.insert(battery_above.end(), battery_above_local.begin(), battery_above_local.end());
                battery_below.insert(battery_below.end(), battery_below_local.begin(), battery_below_local.end());
            }
        }
    }

    // Classify primitives with respect to split
    int n0 = 0, n1 = 0;
    for (int i = 0; i < bestOffset; ++i)
        if (edges[bestAxis][i].type == EdgeType::Start)
            prims0[n0++] = edges[bestAxis][i].primNum;
    for (int i = bestOffset + 1; i < 2 * nPrimitives; ++i)
        if (edges[bestAxis][i].type == EdgeType::End)
            prims1[n1++] = edges[bestAxis][i].primNum;

    STATS_KD::splitRay++;

    // Recursively initialize children nodes
    double tSplit = edges[bestAxis][bestOffset].t;
    AxisAlignedBoundingBox bounds0 = nodeBounds, bounds1 = nodeBounds;
    bounds0.max[bestAxis] = bounds1.min[bestAxis] = tSplit;
    if ((skipL && local_battery.size() < HITCACHEMIN) || (!skipL && battery_below.size() < HITCACHEMIN) ) {
        std::vector<double> dummy;
        buildTree(nodeNum + 1, bounds0, allPrimBounds, prims0, n0, depth - 1, edges,
                  prims0, prims1 + nPrimitives, badRefines, dummy, hybridWeight, bestAxis);
    } else if (!skipL || !skipR){
        buildTree(nodeNum + 1, bounds0, allPrimBounds, prims0, n0, depth - 1, edges,
                  prims0, prims1 + nPrimitives, badRefines, hybridWeight, battery, battery_below, bestAxis, tMin,
                  tSplit,
                  primChance, bestCost);
    }
    else{
        buildTree(nodeNum + 1, bounds0, allPrimBounds, prims0, n0, depth - 1, edges,
                  prims0, prims1 + nPrimitives, badRefines, hybridWeight, battery, local_battery, bestAxis, tMin,
                  tSplit,
                  primChance, bestCost);
    }
    nodes[nodeNum + 1].parent = &nodes[nodeNum];

    int aboveChild = nextFreeNode;
    nodes[nodeNum].InitInterior(bestAxis, aboveChild, tSplit);
    if ((skipL && local_battery.size() < HITCACHEMIN) || (!skipL && battery_above.size() < HITCACHEMIN)) {
        std::vector<double> dummy;
        buildTree(aboveChild, bounds1, allPrimBounds, prims1, n1, depth - 1, edges,
                  prims0, prims1 + nPrimitives, badRefines, dummy, hybridWeight, bestAxis);
    } else if (!skipL || !skipR) {
        buildTree(aboveChild, bounds1, allPrimBounds, prims1, n1, depth - 1, edges,
                  prims0, prims1 + nPrimitives, badRefines, hybridWeight, battery, battery_above, bestAxis, tSplit,
                  tMax, primChance, bestCost);
    }
    else{
        buildTree(aboveChild, bounds1, allPrimBounds, prims1, n1, depth - 1, edges,
                  prims0, prims1 + nPrimitives, badRefines, hybridWeight, battery, local_battery, bestAxis, tSplit,
                  tMax, primChance, bestCost);
    }
    nodes[aboveChild].parent = &nodes[nodeNum];
}

// default construction algorithm for RDH
void KdTreeAccel::buildTreeRDH(int nodeNum, const AxisAlignedBoundingBox &nodeBounds,
                               const std::vector<AxisAlignedBoundingBox> &allPrimBounds, int *primNums, int nPrimitives,
                               int depth, const std::unique_ptr<BoundEdge[]> edges[3], int *prims0, int *prims1,
                               int badRefines, double hybridWeight,
                               const std::unique_ptr<std::vector<RaySegment>> &battery, RayBoundary **rb_stack,
                               int prevSplitAxis, double tMin, double tMax, const std::vector<double> &primChance,
                               double oldCost, BoundaryStack bound) {
    //CHECK_EQ(nodeNum, nextFreeNode);
    assert(nodeNum == nextFreeNode);
    // Get next free node from _nodes_ array
    if (nextFreeNode == nAllocedNodes) {
        int nNewAllocNodes = std::max(2 * nAllocedNodes, 512);
        KdAccelNode *n = new KdAccelNode[nNewAllocNodes];
        if (nAllocedNodes > 0) {
            memcpy(n, nodes, nAllocedNodes * sizeof(KdAccelNode));
            delete[] nodes;
        }
        nodes = n;
        nAllocedNodes = nNewAllocNodes;
    }
    ++nextFreeNode;

    nodes[nodeNum].nodeId = nodeNum;
    ints.emplace_back();
    auto &intstat = ints.back(); // ints[nodeNum]
    intstat.nbPrim = nPrimitives;
    intstat.level = depth;

    // Initialize leaf node if termination criteria met
    if (nPrimitives <= maxPrims || depth == 0) {
        nodes[nodeNum].InitLeaf(primNums, nPrimitives, &primitiveIndices);
        return;
    }

    // Update rays for RDH
    {
        auto ray = Ray();
        double locTMax = tMax;
#pragma omp parallel for default(none) firstprivate(ray) shared(nodeNum, battery, nodeBounds, edges)
        for (int bat = 0; bat < battery->size(); ++bat) {
            auto &test = battery->at(bat).ray;
            ray.origin = test->pos;
            ray.direction = test->dir;
            Vector3d invDir(1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z);
            int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};

            double hitt0 = 0.0;
            double hitt1 = 0.0;

            bool hit = nodeBounds.IntersectP(ray, &hitt0, &hitt1);
            if (!hit) continue;

            // Advance to next child node, possibly enqueue other child
            if (hitt1 < battery->at(bat).tMax) {
                battery->at(bat).tMax = hitt1;
            }
            if (hitt0 > battery->at(bat).tMin) {
                battery->at(bat).tMin = hitt0;
            }
        } // for
    }


    // Choose split axis position for interior node
    int bestAxis = -1, bestOffset = -1;
    double bestCost = 1.0e99;

    //double oldCost = isectCost * double(nPrimitives);
    /*double totalSA = nodeBounds.SurfaceArea();
    double invTotalSA = 1.0 / totalSA;
    Vector3d d = nodeBounds.max - nodeBounds.min;*/

    // Choose which axis to split along
    int axis = nodeBounds.MaximumExtent();
    int retries = 0;

    retrySplit:
    if (splitMethod == SplitMethod::TestSplit && battery)
        std::tie(bestCost, bestAxis, bestOffset) = SplitTest(axis, nodeBounds, allPrimBounds, primNums, nPrimitives,
                                                             edges, battery, rb_stack, primChance, tMax, bound);
        /*else if(splitMethod == SplitMethod::HybridSplit && battery)
            std::tie(bestCost, bestAxis, bestOffset) = SplitHybrid(axis, nodeBounds, allPrimBounds, primNums, nPrimitives,
                                                                   edges, battery, local_battery, primChance, tMax);*/
    else
        std::tie(bestCost, bestAxis, bestOffset) = SplitSAH(axis, nodeBounds, allPrimBounds, primNums, nPrimitives, edges);

    /*if(local_battery.size() < HITCACHEMIN){
            bestCost = bestCostSAH;
            bestAxis = bestAxisSAH;
            bestOffset = bestOffsetSAH;
    }*/
    // Create leaf if no good splits were found
    if (bestAxis == -1 && retries < 2) {
        ++retries;
        if (prevSplitAxis >= 0) {
            axis = (axis + 1) % 3;
            if (retries == 1 && axis == prevSplitAxis) {
                axis = (axis + 1) % 3;
            } else if (retries == 2) {
                axis = prevSplitAxis;
            }
        } else
            axis = (axis + 1) % 3;
        goto retrySplit;
    }

    if (bestCost > oldCost) {
        //++badRefines;
    }

    if ((bestCost > 4.0 * oldCost && nPrimitives < 16) || bestAxis == -1 ||
        badRefines == 3) {
        if (bestCost > 4.0 * oldCost)
            ++STATS_KD::leafHigherCost;
        if (badRefines >= 3)
            ++STATS_KD::leafBadRefine;
        nodes[nodeNum].InitLeaf(primNums, nPrimitives, &primitiveIndices);
        return;
    }

    // Classify primitives with respect to split
    int n0 = 0, n1 = 0;
    for (int i = 0; i < bestOffset; ++i)
        if (edges[bestAxis][i].type == EdgeType::Start)
            prims0[n0++] = edges[bestAxis][i].primNum;
    for (int i = bestOffset + 1; i < 2 * nPrimitives; ++i)
        if (edges[bestAxis][i].type == EdgeType::End)
            prims1[n1++] = edges[bestAxis][i].primNum;

    STATS_KD::splitRay++;

    // Recursively initialize children nodes
    double tSplit = edges[bestAxis][bestOffset].t;
    AxisAlignedBoundingBox bounds0 = nodeBounds, bounds1 = nodeBounds;
    bounds0.max[bestAxis] = bounds1.min[bestAxis] = tSplit;

    /* --- RDH ---
     * When a splitting plane position is decided and a left subtree should be constructed,
     * we check the boundaries in the corresponding axis and
     * mark all ray segments that have to be shortened (are cut by the splitting plane).
     * The original ray boundaries are stored in the stack of
     * signed distances associated with each ray segment data structure.
     * For the marked rays we recompute the boundaries in a corresponding axis.
     * All ray boundaries are then resorted by quicksort before we start to evaluate the RDH.
     */
    // 1. select by axis
    RayBoundary *to_sort = rb_stack[bestAxis];

    // 2. mark rays depending on boundaries
    std::set<RaySegment *> marker;
    for (auto &i: *battery) {
        auto *segment = &i;
        // 4. mark segment to be shortened
        if (segment->tMin < tSplit || tSplit < segment->tMax) {
            marker.insert(segment);
        }
    }

    BoundaryStack bound_new{bound};
    for (int b = bound.leftMin; b < bound.rightMax; ++b) {
        auto &boundary = to_sort[b];

        // process boundary for marked segment
        if (marker.find(boundary.rs) != marker.end()) {
            double tPlane =
                    (tSplit - boundary.rs->ray->pos[bestAxis]) * (-1.0 * boundary.rs->ray->dir[bestAxis]); // inv_dir

            BoundaryFloat tVal = (boundary.flag == 0) ? boundary.rs->tMin : boundary.rs->tMax;
            // 4.1 store original boundaries
            {
                if (tPlane > tVal) {
                    boundary.rs->distances.emplace_back(Vector2d{boundary.rs->tMin, boundary.rs->tMax});
                    if (boundary.flag == 0) boundary.rs->tMin = tSplit;
                    else boundary.rs->tMax = tSplit;
                    boundary.flag = 1;
                } else if (tSplit < tVal) {
                    boundary.rs->distances.emplace_back(Vector2d{boundary.rs->tMin, boundary.rs->tMax});
                    if (boundary.flag == 0) boundary.rs->tMin = tSplit;
                    else boundary.rs->tMax = tSplit;
                    boundary.flag = 1;
                } else { // do nothing
                    //boundary.flag = 1;
                }
            }

            // 4.2 shorten
            /*if(boundary.flag == 0){
                boundary.rs->tMin = tSplit;
            }
            else if(boundary.flag == 2){
                boundary.rs->tMax = tSplit;
            }*/
        }
    }

    // All ray boundaries are then resorted by quicksort before we start to evaluate the RDH.
    std::sort(&to_sort[bound.leftMin], &to_sort[bound.rightMax],
              [](const RayBoundary &b0, const RayBoundary &b1) -> bool {
                  auto b0val = (b0.flag == 0) ? b0.rs->tMin : b0.rs->tMax;
                  auto b1val = (b1.flag == 0) ? b1.rs->tMin : b1.rs->tMax;
                  if (b0val == b1val)
                      return (int) b0.flag < (int) b1.flag;
                  else
                      return b0val < b1val;
              });

    // bound = [0,battery->size() * 2]
    for (int b = bound.leftMin; b < bound.rightMax; ++b) {
        auto &boundary = to_sort[b];

        if (boundary.flag == 0 && boundary.rs->tMin < tSplit) { // TODO: only flag check
            bound_new.leftMax = b;
        } else if (boundary.flag == 2 && boundary.rs->tMax > tSplit) { // TODO: only flag check
            bound_new.rightMin = b;
        } else if (boundary.flag == 1) {
            bound_new.leftMin = b;
            bound_new.rightMax = b;
        }
    }

    // 5. Evaluate RDH

    if (false/*battery_below.size() < HITCACHEMIN*/) {
        std::vector<double> dummy;
        buildTree(nodeNum + 1, bounds0, allPrimBounds, prims0, n0, depth - 1, edges,
                  prims0, prims1 + nPrimitives, badRefines, dummy, hybridWeight, bestAxis);
    } else {
        buildTreeRDH(nodeNum + 1, bounds0, allPrimBounds, prims0, n0, depth - 1, edges,
                     prims0, prims1 + nPrimitives, badRefines, hybridWeight, battery, rb_stack, bestAxis, tMin, tSplit,
                     primChance, bestCost, bound_new);
    }
    nodes[nodeNum + 1].parent = &nodes[nodeNum];

    /*
     * When returning to an interior node after a left subtree has been
     * constructed, we have to change the meaning of ray boundaries that
     * lie on the splitting plane. These originally right boundaries for the
     * left subtree now become the left boundaries for the right subtree to
     * be constructed. The correct restore of boundary values is achieved
     * by storing the depth of the node (together with the signed distance
     * to the splitting plane) where the splitting plane has changed the
     * boundary for each ray. The intervals of ray boundaries are also
     * restored during return from both left child and right child.
     */

    // find stacked distances
    for (int b = bound.leftMin; b < bound.rightMax; ++b) {
        auto &boundary = to_sort[b];
        // process boundary for marked segment
        if (marker.find(boundary.rs) != marker.end()) {
            // 4.1 store original boundaries
            auto dist = boundary.rs->distances.back();
            boundary.rs->distances.pop_back();

            if (boundary.flag == 0) { // left will be split
                boundary.rs->tMin = boundary.rs->tMax;
                boundary.flag = 1;
            } else if (boundary.flag == 1) { // split will be right
                boundary.rs->tMax = dist.v;
                boundary.flag = 2;
            } // right remains unchanged
        }
    }

    // revert distances


    int aboveChild = nextFreeNode;
    nodes[nodeNum].InitInterior(bestAxis, aboveChild, tSplit);
    if (false /*battery_above.size() < HITCACHEMIN*/) {
        std::vector<double> dummy;
        buildTree(aboveChild, bounds1, allPrimBounds, prims1, n1, depth - 1, edges,
                  prims0, prims1 + nPrimitives, badRefines, dummy, hybridWeight, bestAxis);
    } else {
        buildTreeRDH(aboveChild, bounds1, allPrimBounds, prims1, n1, depth - 1, edges,
                     prims0, prims1 + nPrimitives, badRefines, hybridWeight, battery, rb_stack, bestAxis, tSplit, tMax,
                     primChance, bestCost, bound_new);
    }
    nodes[aboveChild].parent = &nodes[nodeNum];
}

bool KdTreeAccel::Intersect(Ray &ray){
    if(restartFromNode && hasRopes){
        return IntersectRopeRestart(ray);
    }
    else if(hasRopes){
        return IntersectRope(ray);
    }
    else{
        return IntersectTrav(ray);
    }
}

bool KdTreeAccel::IntersectStat(RayStat &ray){
    if(restartFromNode && hasRopes){
        return IntersectRopeRestartStat(ray);
    }
    else if(hasRopes){
        return IntersectRopeStat(ray);
    }
    else{
        return IntersectTravStat(ray);
    }
}

bool KdTreeAccel::IntersectTrav(Ray &ray) {
    //ProfilePhase p(Prof::AccelIntersect);
    // Compute initial parametric range of ray inside kd-tree extent
    double tMin = 0.0, tMax = 1.0e99;

    // Neglect initial parametric test, since a ray starts always from inside
    /*if (!bounds.IntersectP(ray, &tMin, &tMax)) {
        return false;
    }*/

    // Prepare to traverse kd-tree for ray
    Vector3d invDir(1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z);
    constexpr int maxTodo = 64;
    KdToDo todo[maxTodo];
    int todoPos = 0;

    // Traverse kd-tree nodes in order for ray
    bool hit = false;
    const KdAccelNode *node = &nodes[0];
    while (node != nullptr) {
        // Bail out if we found a hit closer than the current node
        if (ray.tMax < tMin) break;
        if (!node->IsLeaf()) {
            // Process kd-tree interior node

            // Compute parametric distance along ray to split plane
            int axis = node->SplitAxis();
            double tPlane = (node->SplitPos() - ray.origin[axis]) * invDir[axis];

            // Get node children pointers for ray
            const KdAccelNode *firstChild, *secondChild;
            int belowFirst =
                    (ray.origin[axis] < node->SplitPos()) ||
                    (ray.origin[axis] == node->SplitPos() && ray.direction[axis] <= 0);
            if (belowFirst) {
                firstChild = node + 1;
                secondChild = &nodes[node->AboveChild()];
            } else {
                firstChild = &nodes[node->AboveChild()];
                secondChild = node + 1;
            }

            // Advance to next child node, possibly enqueue other child
            if (tPlane > tMax || tPlane <= 0)
                node = firstChild;
            else if (tPlane < tMin)
                node = secondChild;
            else {
                // Enqueue _secondChild_ in KDtodo list
                todo[todoPos].node = secondChild;
                todo[todoPos].tMin = tPlane;
                todo[todoPos].tMax = tMax;
                ++todoPos;
                node = firstChild;
                tMax = tPlane;
            }
        } else {
            // Check for intersections inside leaf node
            int nPrimitives = node->nPrimitives();
            if (nPrimitives == 1) {
                const std::shared_ptr<Primitive> &p =
                        primitives[node->onePrimitive];

                // Check one primitive inside leaf node
                if (p->globalId != ray.lastIntersected && p->Intersect(ray))
                    hit = true;
            } else {
                for (int i = 0; i < nPrimitives; ++i) {
                    int index =
                            primitiveIndices[node->primitiveIndicesOffset + i];
                    const std::shared_ptr<Primitive> &p = primitives[index];
                    // Check one primitive inside leaf node
                    if (p->globalId != ray.lastIntersected && p->Intersect(ray))
                        hit = true;
                }
            }

            // Grab next node to process from kdtodo list
            if (todoPos > 0) {
                --todoPos;
                node = todo[todoPos].node;
                tMin = todo[todoPos].tMin;
                tMax = todo[todoPos].tMax;
            } else
                break;
        }
    }

    return hit;
}

/*const KdAccelNode * KdAccelNode::getNeighboringNode(const Ray &ray) const
{
    // Check left face.
    double tmax[3];
    bbox.IntersectP(ray, nullptr, &tmax[0], 0);
    bbox.IntersectP(ray, nullptr, &tmax[1], 1);
    bbox.IntersectP(ray, nullptr, &tmax[2], 2);

    int sign[3];
    sign[0] = ray.direction.x < 0.0;
    sign[1] = ray.direction.y < 0.0;
    sign[2] = ray.direction.z < 0.0;

    if(tmax[1] > tmax[0])
    vals = mix(vec2(tmax[1], Y_AXIS + sign[1]), vec2(tmax[1], signs.x), tmax.y > tmax.x);
    vals = mix(vec2(tmax[2, Z_AXIS + sign[2]), vec2(vals.x, vals.y), tmax.z > vals.x);

    int smallestDiffPos = 0;
    double smallestDiff = tmax[0];
    for(int i = 1; i < 3; ++i){
        if(tmax[i] < smallestDiff){
            smallestDiff = tmax[i];
            smallestDiffPos = i;
        }
    }

    if(smallestDiffPos == X_AXIS){
        if(ray.direction.x < 0.0)
            smallestDiffPos = LEFT;
        else
            smallestDiffPos = RIGHT;
    }
    else if (smallestDiffPos == Y_AXIS){
        if(ray.direction.y < 0.0)
            smallestDiffPos = BOTTOM;
        else
            smallestDiffPos = TOP;
    }
    else if (smallestDiffPos == Z_AXIS){
        if(ray.direction.z < 0.0)
            smallestDiffPos = BACK;
        else
            smallestDiffPos = FRONT;
    }

    return ropes[smallestDiffPos];
}*/

const KdAccelNode * KdAccelNode::getNeighboringNode(const Ray &ray) const
{
    // Check left face.
    double tmax[3];
    for(int i = 0; i < 3; ++i)
        bbox.IntersectP(ray, nullptr, &tmax[i], i);

    int smallestDiffPos = 0;
    double smallestDiff = tmax[0];
    for(int i = 1; i < 3; ++i){
        if(tmax[i] < smallestDiff){
            smallestDiff = tmax[i];
            smallestDiffPos = i;
        }
    }

    if(smallestDiffPos == X_AXIS){
        if(ray.direction.x < 0.0)
            smallestDiffPos = LEFT;
        else
            smallestDiffPos = RIGHT;
        //std::cout << "Potential error: x axis very similar" << std::endl;
    }
    else if (smallestDiffPos == Y_AXIS){
        if(ray.direction.y < 0.0)
            smallestDiffPos = BOTTOM;
        else
            smallestDiffPos = TOP;
        //std::cout << "Potential error: z axis very similar" << std::endl;
    }
    else if (smallestDiffPos == Z_AXIS){
        if(ray.direction.z < 0.0)
            smallestDiffPos = BACK;
        else
            smallestDiffPos = FRONT;
        //std::cout << "Potential error: y axis very similar" << std::endl;
    }
    /*Vector3d p_exit = ray.origin + ( tmax[0] * ray.direction );
    Log::console_msg_master(4, "Intersection at: %lf, %lf, %lf -- %lf\n", p_exit.x, p_exit.y, p_exit.z, tmax[0]);
    p_exit = ray.origin + ( tmax[1] * ray.direction );
    Log::console_msg_master(4, "Intersection at: %lf, %lf, %lf -- %lf\n", p_exit.x, p_exit.y, p_exit.z, tmax[1]);
    p_exit = ray.origin + ( tmax[2] * ray.direction );
    Log::console_msg_master(4, "Intersection at: %lf, %lf, %lf -- %lf\n", p_exit.x, p_exit.y, p_exit.z, tmax[2]);

    if ( smallestDiff >= KD_TREE_EPSILON ) {
        std::cout << "ERROR: Node neighbor could not be returned." << std::endl;
        return NULL;
    }
    else{
    if(ropes[smallestDiffPos]!= nullptr){
        std::cout << "Next neighbor --> "<<smallestDiffPos<<"\n";
        ropes[smallestDiffPos]->prettyPrint();
    }
    else{
        std::cout << "No valid neighbor found --> "<<smallestDiffPos<<"\n";
        int i = 0;
        for(auto& rope : ropes){
            if(rope != nullptr) {
                std::cout << "Potential neighbor "<<i<<"\n";
                rope->prettyPrint();
            }
            ++i;
        }
    }*/
    return ropes[smallestDiffPos];
    //}
}

const KdAccelNode * KdAccelNode::getNeighboringNode(const Vector3d &p, const Vector3d &o) const
{
    // Check left face.
    double diff[6];
    diff[0] = fabs( p.x - bbox.min.x );
    diff[1] = fabs( p.z - bbox.max.z );
    diff[2] = fabs( p.x - bbox.max.x );
    diff[3] = fabs( p.z - bbox.min.z );
    diff[4] = fabs( p.y - bbox.max.y );
    diff[5] = fabs( p.y - bbox.min.y );

    int smallestDiffPos = -1;
    double smallestDiff = 1e60;
    for(int i = 0; i < 6; ++i){
        if(diff[i] < smallestDiff){
            smallestDiff = diff[i];
            smallestDiffPos = i;
        }
    }

    if(fabs(diff[0] - diff[2]) < KD_TREE_EPSILON){
        if(o.x < 0.0)
            smallestDiffPos = 0;
        else
            smallestDiffPos = 2;
        std::cout << "Potential error: x axis very similar" << std::endl;
    }
    else if (fabs(diff[1] - diff[3]) < KD_TREE_EPSILON){
        if(o.y < 0.0)
            smallestDiffPos = 0;
        else
            smallestDiffPos = 2;
        std::cout << "Potential error: y axis very similar" << std::endl;
    }
    else if (fabs(diff[4] - diff[5]) < KD_TREE_EPSILON){
        if(o.z < 0.0)
            smallestDiffPos = 0;
        else
            smallestDiffPos = 2;
        std::cout << "Potential error: z axis very similar" << std::endl;
    }

    if ( smallestDiff >= KD_TREE_EPSILON ) {
        std::cout << "ERROR: Node neighbor could not be returned." << std::endl;
        return NULL;
    }
    else{
        return ropes[smallestDiffPos];
    }

    /*if ( smallestDiffPos == 0 ) {
        return ropes[LEFT];
    }
        // Check front face.
    else if ( fabs( p.z - bbox.max.z ) < KD_TREE_EPSILON ) {
        return ropes[FRONT];
    }
        // Check right face.
    else if ( fabs( p.x - bbox.max.x ) < KD_TREE_EPSILON ) {
        return ropes[RIGHT];
    }
        // Check back face.
    else if ( fabs( p.z - bbox.min.z ) < KD_TREE_EPSILON ) {
        return ropes[BACK];
    }
        // Check top face.
    else if ( fabs( p.y - bbox.max.y ) < KD_TREE_EPSILON ) {
        return ropes[TOP];
    }
        // Check bottom face.
    else if ( fabs( p.y - bbox.min.y ) < KD_TREE_EPSILON ) {
        return ropes[BOTTOM];
    }
        // p should be a point on one of the faces of this node's bounding box, but in this case, it isn't.
    else {
        std::cout << "ERROR: Node neighbor could not be returned." << std::endl;
        return NULL;
    }*/
}

bool KdTreeAccel::IntersectRope(Ray &ray) {
    //ProfilePhase p(Prof::AccelIntersect);
    // Compute initial parametric range of ray inside kd-tree extent
    double tMin = 0.0, tMax = 1.0e99;
    double tMinHit = 1.0e99;
    int minId = -1;

    // (-0.796657,-0.591801,0.431659) --> (0.772294,0.467946,0.429638)
    /*ray.origin = {3.2301291053440027, -6.6826513460504922, 18.830288821438948};
    ray.direction = {0.034541292073762692, 0.59865232689214731, -0.80026388813221638};
    ray.lastIntersected = 88;
    if(!ray.pay)
        new RopePayload;
    ((RopePayload*)ray.pay)->lastNode = &nodes[14386];*/
    // Neglect initial parametric test, since a ray starts always from inside
    /*if (!bounds.IntersectP(ray, &tMin, &tMax)) {
        return false;
    }*/

    // Prepare to traverse kd-tree for ray
    Vector3d invDir(1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z);

    // Traverse kd-tree nodes in order for ray
    bool hit = false;
    const KdAccelNode *node = &nodes[0];
    // KdAccelNode *node = &nodes[0];

    node->bbox.IntersectP(ray, &tMin, nullptr);

    while (node != nullptr) {
        // Bail out if we found a hit closer than the current node
        if (ray.tMax < tMin) break;

        // Down traversal
        while(!node->IsLeaf()) {
            // Process kd-tree interior node
            //std::cout << "Moving down to find leaf:\n";
            //node->prettyPrint();

            //bool isLeft = false;
            Vector3d p_entry = ray.origin + ( tMin * ray.direction );

            // Choose left node or right node depending on splitting plane
            if (p_entry[node->SplitAxis()] < node->SplitPos())
                node = node + 1;
            else
                node = &nodes[node->AboveChild()];
        }

        // At a leaf
        // t_exit gets updated on hit with ray.tMax
        {
            // Check for intersections inside leaf node
            int nPrimitives = node->nPrimitives();
            if (nPrimitives == 1) {
                const std::shared_ptr<Primitive> &p =
                        primitives[node->onePrimitive];

                // Check one primitive inside leaf node
                if (p->globalId != ray.lastIntersected && p->Intersect(ray)) {
                    hit = true;
                }

                //Log::console_msg_master(4, " Testing against single prim %d -- %d: %d\n", node->onePrimitive, ray.lastIntersected, hit);

            } else {
                for (int i = 0; i < nPrimitives; ++i) {
                    int index =
                            primitiveIndices[node->primitiveIndicesOffset + i];
                    const std::shared_ptr<Primitive> &p = primitives[index];
                    // Check one primitive inside leaf node
                    if (p->globalId != ray.lastIntersected && p->Intersect(ray)) {
                        hit = true;
                    }
                    //Log::console_msg_master(4, " Testing against prim %d -- %d: %d\n", p->globalId, ray.lastIntersected, hit);
                }
            }
            // Exit leaf

            //std::cout << "Found leaf:\n";
            //node->prettyPrint();
            // TODO: Validate -> Get exit distance from current node
            //int axis = node->SplitAxis();
            //double tPlane = (node->SplitPos() - ray.origin[axis]) * invDir[axis];
            double tMax_l = 0.0; // placeholder, can be null
            // Get exit rope
            if ( node->bbox.IntersectP(ray, &tMax_l, &tMin) ) {
                // Set t_entry to be the entrance point of the next (neighboring) node.
                //t_entry = tmp_t_far;
            }
            else {
                // This should never happen.
                // If it does, then that means we're checking triangles in a node that the ray never intersects.

                break;
            }

            // Get neighboring node using ropes attached to current node.
            //Vector3d p_exit = ray.origin + ( tMin * ray.direction );
            /*auto* */node = node->getNeighboringNode(ray);
            /*if(nodetax!=NULL){
                node = nodetax;
                // TODO: Check if can be unwinded like this
                tMax = (node->SplitPos() - ray.origin[node->SplitAxis()]) * invDir[node->SplitAxis()];

            }
            else {
                // DEBUG
                *//*if (hit == false) {
                    node->bbox.IntersectP(ray, &tMax_l, &tMin);
                    node->getNeighboringNode(ray);
                }*//*
                // END DEBUG
                node = NULL;
            }*/


            // Break if neighboring node not found, meaning we've exited the kd-tree.
            if ( node == NULL ) {
                break;
            }
        }
    }

    /*if(!hit)
        Log::console_msg_master(4, "Couldn't find roped intersection for ray (%lf,%lf,%lf) --> (%lf,%lf,%lf) [%lf,%lf]\n",
               ray.origin.x, ray.origin.y, ray.origin.z, ray.direction.x, ray.direction.y, ray.direction.z, tMin, ray.tMax);
*/
    return hit;
}

bool KdTreeAccel::IntersectRopeStat(RayStat &ray) {
    //ProfilePhase p(Prof::AccelIntersect);
    time_rope.ReStart();
    // Compute initial parametric range of ray inside kd-tree extent
    double tMin = 0.0, tMax = 1.0e99;
    double tMinHit = 1.0e99;
    int minId = -1;
    ray.traversedNodes.clear();

    // Neglect initial parametric test, since a ray starts always from inside
    /*if (!bounds.IntersectP(ray, &tMin, &tMax)) {
        return false;
    }
*/
    // Prepare to traverse kd-tree for ray
    Vector3d invDir(1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z);

    // Traverse kd-tree nodes in order for ray
    bool hit = false;
    const KdAccelNode *node = &nodes[0];

    //node->bbox.IntersectP(ray, &tMin, nullptr);

    ray.stats.traversalSteps++;
    while (node != nullptr) {
        ray.traversedNodes.push_back(node->nodeId);

        // Bail out if we found a hit closer than the current node

        if (ray.tMax < tMin) break;

        // Down traversal
        while(!node->IsLeaf()) {
            ints[node->nodeId].nbChecks += ray.stats.traversalSteps;
            // Process kd-tree interior node
            //std::cout << "Moving down to find leaf:\n";
            //node->prettyPrint();

            ray.stats.nbDownTest++; // interior test
            ray.stats.traversalSteps++;
            ints[node->nodeId].nbIntersects++;

            //bool isLeft = false;
            Vector3d p_entry = ray.origin + ( tMin * ray.direction );

            // Advance to next child node, possibly enqueue other child
            if (p_entry[node->SplitAxis()] < node->SplitPos())
                node = node + 1;
            else
                node = &nodes[node->AboveChild()];
            ray.traversedNodes.push_back(node->nodeId);
        }

        ints[node->nodeId].nbChecks += ray.stats.traversalSteps;

        // At a leaf
        // t_exit gets updated on hit with ray.tMax
        {
            // Check for intersections inside leaf node
            int nPrimitives = node->nPrimitives();
            if (nPrimitives == 1) {
                ray.stats.traversalSteps++;

                const std::shared_ptr<Primitive> &p =
                        primitives[node->onePrimitive];

                // Check one primitive inside leaf node
                if (p->globalId != ray.lastIntersected && p->IntersectStat(ray)) {
                    hit = true;
                }

            } else {
                for (int i = 0; i < nPrimitives; ++i) {
                    ray.stats.traversalSteps++;
                    int index =
                            primitiveIndices[node->primitiveIndicesOffset + i];
                    const std::shared_ptr<Primitive> &p = primitives[index];
                    // Check one primitive inside leaf node
                    if (p->globalId != ray.lastIntersected && p->IntersectStat(ray)) {
                        hit = true;
                    }
                }
            }
            // Exit leaf

            //std::cout << "Found leaf:\n";
            //node->prettyPrint();
            // TODO: Validate -> Get exit distance from current node
            //int axis = node->SplitAxis();
            //double tPlane = (node->SplitPos() - ray.origin[axis]) * invDir[axis];
            double tMax_l = 0.0; // placeholder, can be null
            // Get exit rope
            if ( node->bbox.IntersectP(ray, &tMax_l, &tMin) ) {
                // Set t_entry to be the entrance point of the next (neighboring) node.
                //t_entry = tmp_t_far;
            }
            else {
                // This should never happen.
                // If it does, then that means we're checking triangles in a node that the ray never intersects.
                break;
            }

            // Get neighboring node using ropes attached to current node.
            //Vector3d p_exit = ray.origin + ( tMin * ray.direction );
            node = node->getNeighboringNode(ray);

            // Break if neighboring node not found, meaning we've exited the kd-tree.
            if ( node == NULL ) {
                break;
            }

            ray.stats.traversalSteps++; // inc at the end, as we do it once in the beginning out of the loop
        }
    }

    /*if(!hit)
        Log::console_msg_master(4, "Couldn't find roped intersection for ray (%lf,%lf,%lf) --> (%lf,%lf,%lf) [%lf,%lf]\n",
               ray.origin.x, ray.origin.y, ray.origin.z, ray.direction.x, ray.direction.y, ray.direction.z, tMin, ray.tMax);
*/
    perRayCount += ray.stats;
    ray.ResetStats();
    Profiling::intersectStatsKDRope.AddTime(time_rope.ElapsedMs());

    return hit;
}

bool KdTreeAccel::IntersectRopeRestart(Ray &ray) {
    //ProfilePhase p(Prof::AccelIntersect);
    // Compute initial parametric range of ray inside kd-tree extent
    double tMin = 0.0, tMax = 1.0e99;
    double tMinHit = 1.0e99;
    int minId = -1;

    // (-0.796657,-0.591801,0.431659) --> (0.772294,0.467946,0.429638)
    /*ray.origin = {3.2301291053440027, -6.6826513460504922, 18.830288821438948};
    ray.direction = {0.034541292073762692, 0.59865232689214731, -0.80026388813221638};
    ray.lastIntersected = 88;
    if(!ray.pay)
        new RopePayload;
    ((RopePayload*)ray.pay)->lastNode = &nodes[14386];*/
    // Neglect initial parametric test, since a ray starts always from inside
    /*if (!bounds.IntersectP(ray, &tMin, &tMax)) {
        return false;
    }*/

    RopePayload* ropePay = ray.pay ? (RopePayload*)ray.pay : nullptr;
    // Prepare to traverse kd-tree for ray
    Vector3d invDir(1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z);

    // Traverse kd-tree nodes in order for ray
    bool hit = false;
    const KdAccelNode *node;
    if(ropePay && ropePay->lastNode) {
        node = ropePay->lastNode;
#if defined(DEBUG)
        ((RopePayload *) ray.pay)->lastRay = ray;
#endif
        size_t nHops = 0;
        while(node != nullptr && !node->bbox.IsInside(ray.origin)){// not in node
            node = node->parent;
            nHops++;
        }
        if(node == nullptr)
            node = &nodes[0];
        /*else
            std::cout << "Sane starting point after " << nHops << " links\n";*/
    }
    else
        node = &nodes[0];
    // KdAccelNode *node = &nodes[0];

    node->bbox.IntersectP(ray, &tMin, nullptr);

    while (node != nullptr) {
        // Bail out if we found a hit closer than the current node
        if (ray.tMax < tMin) break;

        // Down traversal
        while(!node->IsLeaf()) {
            // Process kd-tree interior node
            //std::cout << "Moving down to find leaf:\n";
            //node->prettyPrint();

            //bool isLeft = false;
            Vector3d p_entry = ray.origin + ( tMin * ray.direction );

            // Choose left node or right node depending on splitting plane
            if (p_entry[node->SplitAxis()] < node->SplitPos())
                node = node + 1;
            else
                node = &nodes[node->AboveChild()];
        }

        // At a leaf
        // t_exit gets updated on hit with ray.tMax
        {
            // Check for intersections inside leaf node
            int nPrimitives = node->nPrimitives();
            if (nPrimitives == 1) {
                const std::shared_ptr<Primitive> &p =
                        primitives[node->onePrimitive];

                // Check one primitive inside leaf node
                if (p->globalId != ray.lastIntersected && p->Intersect(ray)) {
                    hit = true;
                    ropePay->lastNode = node;
                    /*if(ray.pay && ray.hardHit.hitId == p->globalId) {
                        double minLocal = 1e99;
                        node->bbox.IntersectP(ray, &minLocal, nullptr);
                        if (minId == p->globalId && tMinHit > minLocal) {
                            tMinHit = minLocal;
                            ((RopePayload *) ray.pay)->lastNode = node;
                        } else if(minId != p->globalId){
                            tMinHit = minLocal;
                            minId = p->globalId;
                            ((RopePayload *) ray.pay)->lastNode = node;
                        }
                    }*/
                }

                //Log::console_msg_master(4, " Testing against single prim %d -- %d: %d\n", node->onePrimitive, ray.lastIntersected, hit);

            } else {
                for (int i = 0; i < nPrimitives; ++i) {
                    int index =
                            primitiveIndices[node->primitiveIndicesOffset + i];
                    const std::shared_ptr<Primitive> &p = primitives[index];
                    // Check one primitive inside leaf node
                    if (p->globalId != ray.lastIntersected && p->Intersect(ray)) {
                        hit = true;
                        ropePay->lastNode = node;
                        /*if(ray.pay && ray.hardHit.hitId == p->globalId) {
                            double minLocal = 1e99;
                            node->bbox.IntersectP(ray, &minLocal, nullptr);
                            if (minId == p->globalId && tMinHit > minLocal) {
                                tMinHit = minLocal;
                                ((RopePayload *) ray.pay)->lastNode = node;
                            } else if(minId != p->globalId){
                                tMinHit = minLocal;
                                minId = p->globalId;
                                ((RopePayload *) ray.pay)->lastNode = node;
                            }
                        }*/
                    }
                    //Log::console_msg_master(4, " Testing against prim %d -- %d: %d\n", p->globalId, ray.lastIntersected, hit);
                }
            }
            // Exit leaf

            //std::cout << "Found leaf:\n";
            //node->prettyPrint();
            // TODO: Validate -> Get exit distance from current node
            //int axis = node->SplitAxis();
            //double tPlane = (node->SplitPos() - ray.origin[axis]) * invDir[axis];
            double tMax_l = 0.0; // placeholder, can be null
            // Get exit rope
            if ( node->bbox.IntersectP(ray, &tMax_l, &tMin) ) {
                // Set t_entry to be the entrance point of the next (neighboring) node.
                //t_entry = tmp_t_far;
            }
            else {
                // This should never happen.
                // If it does, then that means we're checking triangles in a node that the ray never intersects.
                // TODO: DEBUG
                /*node->prettyPrint();
                node->bbox.IntersectP(ray, &tMax_l, &tMin);
                double tmax_loc = -1.0;
                node->bbox.IntersectP(ray, nullptr, &tmax_loc, 0);
                Vector3d p_exit = ray.origin + ( tmax_loc * ray.direction );
                Log::console_msg_master(4, "[0] nIntersection at: %lf, %lf, %lf -- %lf\n", p_exit.x, p_exit.y, p_exit.z, tmax_loc);
                tmax_loc = -1.0;
                node->bbox.IntersectP(ray, nullptr, &tmax_loc, 1);
                p_exit = ray.origin + ( tmax_loc * ray.direction );
                Log::console_msg_master(4, "[1] nIntersection at: %lf, %lf, %lf -- %lf\n", p_exit.x, p_exit.y, p_exit.z, tmax_loc);
                tmax_loc = -1.0;
                node->bbox.IntersectP(ray, nullptr, &tmax_loc, 2);
                p_exit = ray.origin + ( tmax_loc * ray.direction );
                Log::console_msg_master(4, "[2] nIntersection at: %lf, %lf, %lf -- %lf\n", p_exit.x, p_exit.y, p_exit.z, tmax_loc);*/

                break;
            }

            // Get neighboring node using ropes attached to current node.
            //Vector3d p_exit = ray.origin + ( tMin * ray.direction );
            /*auto* */node = node->getNeighboringNode(ray);
            /*if(nodetax!=NULL){
                node = nodetax;
                // TODO: Check if can be unwinded like this
                tMax = (node->SplitPos() - ray.origin[node->SplitAxis()]) * invDir[node->SplitAxis()];

            }
            else {
                // DEBUG
                *//*if (hit == false) {
                    node->bbox.IntersectP(ray, &tMax_l, &tMin);
                    node->getNeighboringNode(ray);
                }*//*
                // END DEBUG
                node = NULL;
            }*/


            // Break if neighboring node not found, meaning we've exited the kd-tree.
            if ( node == NULL ) {
                break;
            }
        }
    }

    /*if(!hit)
        Log::console_msg_master(4, "Couldn't find roped intersection for ray (%lf,%lf,%lf) --> (%lf,%lf,%lf) [%lf,%lf]\n",
               ray.origin.x, ray.origin.y, ray.origin.z, ray.direction.x, ray.direction.y, ray.direction.z, tMin, ray.tMax);
*/
    return hit;
}

bool KdTreeAccel::IntersectRopeRestartStat(RayStat &ray) {
    //ProfilePhase p(Prof::AccelIntersect);
    time_roperestart.ReStart();

    // Compute initial parametric range of ray inside kd-tree extent
    double tMin = 0.0, tMax = 1.0e99;
    double tMinHit = 1.0e99;
    int minId = -1;
    ray.traversedNodes.clear();

    // Neglect initial parametric test, since a ray starts always from inside
    /*if (!bounds.IntersectP(ray, &tMin, &tMax)) {
        return false;
    }
*/
    // Prepare to traverse kd-tree for ray
    Vector3d invDir(1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z);
    RopePayload* ropePay = ray.pay ? (RopePayload*)ray.pay : nullptr;

    // Traverse kd-tree nodes in order for ray
    bool hit = false;
    const KdAccelNode *node;
    if(ray.pay && ropePay->lastNode) {
        node = ropePay->lastNode;
        size_t nHops = 0;
        //auto tmpNode = node;
        while(node != nullptr && !node->bbox.IsInside(ray.origin)){// not in node
            ray.traversedNodes.push_back(node->nodeId);
            //node = node->parent;
            node = node->getNeighboringNode(ray);
            nHops++;
            ray.stats.traversalSteps++;
            ray.stats.nbUpTest++;
        }
        /*size_t nHopsNeigh = 0;
        while(tmpNode != nullptr && !tmpNode->bbox.IsInside(ray.origin)){// not in node
            //ray.traversedNodes.push_back(node->nodeId);
            tmpNode = tmpNode->getNeighboringNode(ray);
            nHopsNeigh++;
        }*/
        if(node == nullptr)
            node = &nodes[0];
        /*else
            std::cout << "Sane starting point after " << nHops << " links\n";*/
    }
    else
        node = &nodes[0];

    node->bbox.IntersectP(ray, &tMin, nullptr);

    while (node != nullptr) {
        ray.stats.traversalSteps++;
        ray.traversedNodes.push_back(node->nodeId);

        // Bail out if we found a hit closer than the current node

        if (ray.tMax < tMin) break;

        // Down traversal
        while(!node->IsLeaf()) {
            ints[node->nodeId].nbChecks += ray.stats.traversalSteps;
            // Process kd-tree interior node
            //std::cout << "Moving down to find leaf:\n";
            //node->prettyPrint();
            ray.stats.nbDownTest++; // interior test
            ray.stats.traversalSteps++;
            ints[node->nodeId].nbIntersects++;

            //bool isLeft = false;
            Vector3d p_entry = ray.origin + ( tMin * ray.direction );

            // Advance to next child node, possibly enqueue other child
            if (p_entry[node->SplitAxis()] < node->SplitPos())
                node = node + 1;
            else
                node = &nodes[node->AboveChild()];
            ray.traversedNodes.push_back(node->nodeId);
        }

        ints[node->nodeId].nbChecks += ray.stats.traversalSteps;

        // At a leaf
        // t_exit gets updated on hit with ray.tMax
        {
            // Check for intersections inside leaf node
            int nPrimitives = node->nPrimitives();
            if (nPrimitives == 1) {
                ray.stats.traversalSteps++;

                const std::shared_ptr<Primitive> &p =
                        primitives[node->onePrimitive];

                // Check one primitive inside leaf node
                if (p->globalId != ray.lastIntersected && p->IntersectStat(ray)) {
                    hit = true;
                    ropePay->lastNode = node;
                }

            } else {
                for (int i = 0; i < nPrimitives; ++i) {
                    ray.stats.traversalSteps++;
                    int index =
                            primitiveIndices[node->primitiveIndicesOffset + i];
                    const std::shared_ptr<Primitive> &p = primitives[index];
                    // Check one primitive inside leaf node
                    if (p->globalId != ray.lastIntersected && p->IntersectStat(ray)) {
                        hit = true;
                        if(ropePay && ropePay->lastNode)
                            ropePay->lastNode = node;
                    }
                }
            }
            // Exit leaf

            //std::cout << "Found leaf:\n";
            //node->prettyPrint();
            // TODO: Validate -> Get exit distance from current node
            //int axis = node->SplitAxis();
            //double tPlane = (node->SplitPos() - ray.origin[axis]) * invDir[axis];
            double tMax_l = 0.0; // placeholder, can be null
            // Get exit rope
            if ( node->bbox.IntersectP(ray, &tMax_l, &tMin) ) {
                // Set t_entry to be the entrance point of the next (neighboring) node.
                //t_entry = tmp_t_far;
            }
            else {
                // This should never happen.
                // If it does, then that means we're checking triangles in a node that the ray never intersects.
                break;
            }

            // Get neighboring node using ropes attached to current node.
            //Vector3d p_exit = ray.origin + ( tMin * ray.direction );
            node = node->getNeighboringNode(ray);

            // Break if neighboring node not found, meaning we've exited the kd-tree.
            if ( node == NULL ) {
                break;
            }
        }
    }

    /*if(!hit)
        Log::console_msg_master(4, "Couldn't find roped intersection for ray (%lf,%lf,%lf) --> (%lf,%lf,%lf) [%lf,%lf]\n",
               ray.origin.x, ray.origin.y, ray.origin.z, ray.direction.x, ray.direction.y, ray.direction.z, tMin, ray.tMax);
*/
    perRayCount += ray.stats;
    ray.ResetStats();

    Profiling::intersectStatsKDRopeRestart.AddTime(time_roperestart.ElapsedMs());

    return hit;
}

bool KdTreeAccel::IntersectTravStat(RayStat &ray) {
    //ProfilePhase p(Prof::AccelIntersect);
    time_trav.ReStart();
    // Compute initial parametric range of ray inside kd-tree extent
    double tMin = 0.0, tMax = 1.0e99;
    ray.traversedNodes.clear();

    // Neglect initial parametric test, since a ray starts always from inside
    /*if (!bounds.IntersectP(ray, &tMin, &tMax)) {
        return false;
    }*/

    // Prepare to traverse kd-tree for ray
    Vector3d invDir(1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z);
    constexpr int maxTodo = 64;
    KdToDo todo[maxTodo];
    int todoPos = 0;

    // Traverse kd-tree nodes in order for ray
    bool hit = false;
    const KdAccelNode *node = &nodes[0];
    while (node != nullptr) {
        ray.traversedNodes.push_back(node->nodeId);
        ints[node->nodeId].nbChecks += ray.stats.traversalSteps;

        // Bail out if we found a hit closer than the current node
        if (ray.tMax < tMin) break;
        if (!node->IsLeaf()) {
            // Process kd-tree interior node
            ray.stats.traversalSteps++;
            ray.stats.nbDownTest++; // interior test
            ints[node->nodeId].nbIntersects++;

            // Compute parametric distance along ray to split plane
            int axis = node->SplitAxis();
            double tPlane = (node->SplitPos() - ray.origin[axis]) * invDir[axis];

            // Get node children pointers for ray
            const KdAccelNode *firstChild, *secondChild;
            int belowFirst =
                    (ray.origin[axis] < node->SplitPos()) ||
                    (ray.origin[axis] == node->SplitPos() && ray.direction[axis] <= 0);
            if (belowFirst) {
                firstChild = node + 1;
                secondChild = &nodes[node->AboveChild()];
            } else {
                firstChild = &nodes[node->AboveChild()];
                secondChild = node + 1;
            }

            // Advance to next child node, possibly enqueue other child
            if (tPlane > tMax || tPlane <= 0)
                node = firstChild;
            else if (tPlane < tMin)
                node = secondChild;
            else {
                // Enqueue _secondChild_ in KDtodo list
                todo[todoPos].node = secondChild;
                todo[todoPos].tMin = tPlane;
                todo[todoPos].tMax = tMax;
                ++todoPos;
                node = firstChild;
                tMax = tPlane;
            }
        } else {
            // Check for intersections inside leaf node
            int nPrimitives = node->nPrimitives();
            if (nPrimitives == 1) {
                const std::shared_ptr<Primitive> &p =
                        primitives[node->onePrimitive];

                // Check one primitive inside leaf node
                ray.stats.traversalSteps++;
                if (p->globalId != ray.lastIntersected && p->IntersectStat(ray))
                    hit = true;
            } else {
                for (int i = 0; i < nPrimitives; ++i) {
                    int index =
                            primitiveIndices[node->primitiveIndicesOffset + i];
                    const std::shared_ptr<Primitive> &p = primitives[index];
                    // Check one primitive inside leaf node
                    ray.stats.traversalSteps++;
                    if (p->globalId != ray.lastIntersected && p->IntersectStat(ray))
                        hit = true;
                }
            }

            // Grab next node to process from kdtodo list
            if (todoPos > 0) {
                --todoPos;
                node = todo[todoPos].node;
                tMin = todo[todoPos].tMin;
                tMax = todo[todoPos].tMax;
            } else
                break;
        }
    }

    perRayCount += ray.stats;
    ray.ResetStats();

    Profiling::intersectStatsKD.AddTime(time_trav.ElapsedMs());

    return hit;
}

RTStats KdTreeAccel::IntersectT(Ray &ray) {
    //ProfilePhase p(Prof::AccelIntersect);
    // Compute initial parametric range of ray inside kd-tree extent
    double tMin = 0.0, tMax = 1.0e99;
    /*if (!bounds.IntersectP(ray, &tMin, &tMax)) {
        return false;
    }*/

    // Prepare to traverse kd-tree for ray
    Vector3d invDir(1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z);
    constexpr int maxTodo = 64;
    KdToDo todo[maxTodo];
    int todoPos = 0;

    // Traverse kd-tree nodes in order for ray
    bool hit = false;
    const KdAccelNode *node = &nodes[0];

    RTStats stat;
    Chronometer time_trav(false);
    Chronometer time_int(false);
    while (node != nullptr) {
        // Bail out if we found a hit closer than the current node
        if (ray.tMax < tMin) break;
        if (!node->IsLeaf()) {

            // Process kd-tree interior node
            time_trav.Start();
            stat.nTraversedInner++;

            // Compute parametric distance along ray to split plane
            int axis = node->SplitAxis();
            double tPlane = (node->SplitPos() - ray.origin[axis]) * invDir[axis];

            // Get node children pointers for ray
            const KdAccelNode *firstChild, *secondChild;
            int belowFirst =
                    (ray.origin[axis] < node->SplitPos()) ||
                    (ray.origin[axis] == node->SplitPos() && ray.direction[axis] <= 0);
            if (belowFirst) {
                firstChild = node + 1;
                secondChild = &nodes[node->AboveChild()];
            } else {
                firstChild = &nodes[node->AboveChild()];
                secondChild = node + 1;
            }

            // Advance to next child node, possibly enqueue other child
            if (tPlane > tMax || tPlane <= 0)
                node = firstChild;
            else if (tPlane < tMin)
                node = secondChild;
            else {
                // Enqueue _secondChild_ in KDtodo list
                todo[todoPos].node = secondChild;
                todo[todoPos].tMin = tPlane;
                todo[todoPos].tMax = tMax;
                ++todoPos;
                node = firstChild;
                tMax = tPlane;
            }
            time_trav.Stop();
        } else {
            time_int.Start();

            // Check for intersections inside leaf node
            int nPrimitives = node->nPrimitives();

            stat.nTraversedLeaves++;
            stat.nIntersections += nPrimitives;
            if (nPrimitives == 1) {
                const std::shared_ptr<Primitive> &p =
                        primitives[node->onePrimitive];

                // Check one primitive inside leaf node
                if (p->globalId != ray.lastIntersected && p->IntersectT(ray))
                    hit = true;
            } else {
                for (int i = 0; i < nPrimitives; ++i) {
                    int index =
                            primitiveIndices[node->primitiveIndicesOffset + i];
                    const std::shared_ptr<Primitive> &p = primitives[index];
                    // Check one primitive inside leaf node
                    if (p->globalId != ray.lastIntersected && p->IntersectT(ray))
                        hit = true;
                }
            }

            // Grab next node to process from kdtodo list
            if (todoPos > 0) {
                --todoPos;
                node = todo[todoPos].node;
                tMin = todo[todoPos].tMin;
                tMax = todo[todoPos].tMax;
            } else
                break;
            time_int.Stop();
        }
    }

    stat.timeInt = time_int.ElapsedMs();
    stat.timeTrav = time_trav.ElapsedMs();
    return stat;
}

// Public-facing wrapper method.
/*bool KdTreeAccel::singleRayStacklessIntersect( Ray &ray ) const
{
    // Compute initial parametric range of ray inside kd-tree extent
    double tMin = 0.0, tMax = 1.0e99;

    // Neglect initial parametric test, since a ray starts always from inside
    *//*if (!bounds.IntersectP(ray, &tMin, &tMax)) {
        return false;
    }*//*

    //if ( intersects_root_node_bounding_box ) {
        bool hit = singleRayStacklessIntersect( root, ray_o, ray_dir, t_near, t_far, normal );
        if ( hit ) {
            t = t_far;
            hit_point = ray_o + ( t * ray_dir );
        }
        return hit;
    *//*}
    else {
        return false;
    }*//*
}*/

// Private method that does the heavy lifting.
/*bool KdTreeAccel::singleRayStacklessIntersect( KdAccelNode *curr_node, Ray &ray ) const
{
    bool intersection_detected = false;

    float t_entry_prev = -INFINITY;
    while ( t_entry < t_exit && t_entry > t_entry_prev ) {
        t_entry_prev = t_entry;

        // Down traversal - Working our way down to a leaf node.
        glm::vec3 p_entry = ray_o + ( t_entry * ray_dir );
        while ( !curr_node->is_leaf_node ) {
            curr_node = curr_node->isPointToLeftOfSplittingPlane( p_entry ) ? curr_node->left : curr_node->right;
        }

        // We've reached a leaf node.
        // Check intersection with triangles contained in current leaf node.
        for ( int i = 0; i < curr_node->num_tris; ++i ) {
            glm::vec3 tri = tris[curr_node->tri_indices[i]];
            glm::vec3 v0 = verts[( int )tri[0]];
            glm::vec3 v1 = verts[( int )tri[1]];
            glm::vec3 v2 = verts[( int )tri[2]];

            // Perform ray/triangle intersection test.
            float tmp_t = INFINITY;
            glm::vec3 tmp_normal( 0.0f, 0.0f, 0.0f );
            bool intersects_tri = Intersections::triIntersect( ray_o, ray_dir, v0, v1, v2, tmp_t, tmp_normal );

            if ( intersects_tri ) {
                if ( tmp_t < t_exit ) {
                    intersection_detected = true;
                    t_exit = tmp_t;
                    normal = tmp_normal;
                }
            }
        }

        // Compute distance along ray to exit current node.
        float tmp_t_near, tmp_t_far;
        bool intersects_curr_node_bounding_box = Intersections::aabbIntersect( curr_node->bbox, ray_o, ray_dir, tmp_t_near, tmp_t_far );
        if ( intersects_curr_node_bounding_box ) {
            // Set t_entry to be the entrance point of the next (neighboring) node.
            t_entry = tmp_t_far;
        }
        else {
            // This should never happen.
            // If it does, then that means we're checking triangles in a node that the ray never intersects.
            break;
        }

        // Get neighboring node using ropes attached to current node.
        glm::vec3 p_exit = ray_o + ( t_entry * ray_dir );
        curr_node = curr_node->getNeighboringNode( p_exit );

        // Break if neighboring node not found, meaning we've exited the kd-tree.
        if ( curr_node == NULL ) {
            break;
        }
    }

    return intersection_detected;
}*/

KdTreeAccel::KdTreeAccel(KdTreeAccel &&src)
noexcept: isectCost(src.isectCost),
          traversalCost(src.traversalCost),
          maxPrims(src.maxPrims),
          emptyBonus(src.emptyBonus),
          primitives(std::move(src.primitives)),
          primitiveIndices(std::move(src.primitiveIndices)) {
    nodes = src.nodes;
    nAllocedNodes = src.nAllocedNodes;
    nextFreeNode = src.nextFreeNode;
    src.nodes = nullptr;
    bb = src.bb;
}

KdTreeAccel::KdTreeAccel(
        const KdTreeAccel &src) noexcept: isectCost(src.isectCost),
                                          traversalCost(src.traversalCost),
                                          maxPrims(src.maxPrims),
                                          emptyBonus(src.emptyBonus),
                                          primitives(src.primitives),
                                          primitiveIndices(src.primitiveIndices) {
    if (nodes)
        exit(44);
}

KdTreeAccel &KdTreeAccel::operator=(const KdTreeAccel &src) noexcept {
    if (nodes)
        exit(44);

    return *this;
}

std::ostream &operator<<(std::ostream &os, KdTreeAccel::SplitMethod split_type) {
    switch (split_type) {
        case KdTreeAccel::SplitMethod::SAH : return os << "SAH" ;
        case KdTreeAccel::SplitMethod::ProbSplit: return os << "ProbSplit";
        case KdTreeAccel::SplitMethod::ProbHybrid: return os << "ProbHybrid";
        case KdTreeAccel::SplitMethod::TestSplit: return os << "TestSplit";
        case KdTreeAccel::SplitMethod::HybridSplit: return os << "HybridSplit";
        case KdTreeAccel::SplitMethod::HybridBin: return os << "HybridBin";
        default: return os << "Unknown";
            // omit default case to trigger compiler warning for missing cases
    };
    return os << static_cast<std::uint16_t>(split_type);
}

void KdTreeAccel::optimizeRopes( KdAccelNode *ropes[], AxisAlignedBoundingBox bbox )
{
    // Loop through ropes of all faces of node bounding box.
    for ( int i = 0; i < 6; ++i ) {
        KdAccelNode *rope_node = ropes[i];

        if ( rope_node == NULL ) {
            continue;
        }

        // Process until leaf node is reached.
        // The optimization - We try to push the ropes down into the tree as far as possible
        // instead of just having the ropes point to the roots of neighboring subtrees.
        while ( !rope_node->IsLeaf() ) {

            // Get node children pointers for recursive traversal
            KdAccelNode *leftChild, *rightChild;
            leftChild = rope_node + 1;
            rightChild = &nodes[rope_node->AboveChild()];

            // Case I.

            if ( i == LEFT || i == RIGHT ) {

                // Case I-A.

                // Handle parallel split plane case.
                if ( rope_node->SplitAxis() == X_AXIS ) {
                    rope_node = ( i == LEFT ) ? rightChild : leftChild;
                }

                    // Case I-B.

                else if ( rope_node->SplitAxis() == Y_AXIS ) {
                    if ( rope_node->SplitPos() < ( bbox.min.y - KD_TREE_EPSILON ) ) {
                        rope_node = rightChild;
                    }
                    else if ( rope_node->SplitPos() > ( bbox.max.y + KD_TREE_EPSILON ) ) {
                        rope_node = leftChild;
                    }
                    else {
                        break;
                    }
                }

                    // Case I-C.

                    // Split plane is Z_AXIS.
                else {
                    if ( rope_node->SplitPos() < ( bbox.min.z - KD_TREE_EPSILON ) ) {
                        rope_node = rightChild;
                    }
                    else if ( rope_node->SplitPos() > ( bbox.max.z + KD_TREE_EPSILON ) ) {
                        rope_node = leftChild;
                    }
                    else {
                        break;
                    }
                }
            }

                // Case II.

            else if ( i == FRONT || i == BACK ) {

                // Case II-A.

                // Handle parallel split plane case.
                if ( rope_node->SplitAxis() == Z_AXIS ) {
                    rope_node = ( i == BACK ) ? rightChild : leftChild;
                }

                    // Case II-B.

                else if ( rope_node->SplitAxis() == X_AXIS ) {
                    if ( rope_node->SplitPos() < ( bbox.min.x - KD_TREE_EPSILON ) ) {
                        rope_node = rightChild;
                    }
                    else if ( rope_node->SplitPos() > ( bbox.max.x + KD_TREE_EPSILON ) ) {
                        rope_node = leftChild;
                    }
                    else {
                        break;
                    }
                }

                    // Case II-C.

                    // Split plane is Y_AXIS.
                else {
                    if ( rope_node->SplitPos() < ( bbox.min.y - KD_TREE_EPSILON ) ) {
                        rope_node = rightChild;
                    }
                    else if ( rope_node->SplitPos() > ( bbox.max.y + KD_TREE_EPSILON ) ) {
                        rope_node = leftChild;
                    }
                    else {
                        break;
                    }
                }
            }

                // Case III.

                // TOP and BOTTOM.
            else {

                // Case III-A.

                // Handle parallel split plane case.
                if ( rope_node->SplitAxis() == Y_AXIS ) {
                    rope_node = ( i == BOTTOM ) ? rightChild : leftChild;
                }

                    // Case III-B.

                else if ( rope_node->SplitAxis() == Z_AXIS ) {
                    if ( rope_node->SplitPos() < ( bbox.min.z - KD_TREE_EPSILON ) ) {
                        rope_node = rightChild;
                    }
                    else if ( rope_node->SplitPos() > ( bbox.max.z + KD_TREE_EPSILON ) ) {
                        rope_node = leftChild;
                    }
                    else {
                        break;
                    }
                }

                    // Case III-C.

                    // Split plane is X_AXIS.
                else {
                    if ( rope_node->SplitPos() < ( bbox.min.x - KD_TREE_EPSILON ) ) {
                        rope_node = rightChild;
                    }
                    else if ( rope_node->SplitPos() > ( bbox.max.x + KD_TREE_EPSILON ) ) {
                        rope_node = leftChild;
                    }
                    else {
                        break;
                    }
                }
            }
        }
    }
}

// recursively attach ropes from root node in DF order
void KdTreeAccel::attachRopes(KdAccelNode *current, KdAccelNode *ropes[], bool optimize) {
// Base case.
    if ( current->IsLeaf() ) {
        //std::cout<<curr_node->id<<": "<<std::endl;
        for ( int i = 0; i < 6; ++i ) {
            current->ropes[i] = ropes[i];
        }
    }
    else {

        if(optimize)
            optimizeRopes( ropes, current->bbox );

        AABBFace SL, SR;
        if ( current->SplitAxis() == X_AXIS ) {
            SL = LEFT;
            SR = RIGHT;
        }
        else if ( current->SplitAxis() == Y_AXIS ) {
            SL = BOTTOM;
            SR = TOP;
        }
            // Split plane is Z_AXIS
        else {
            SL = BACK;
            SR = FRONT;
        }

        KdAccelNode* RS_left[6];
        KdAccelNode* RS_right[6];
        for ( int i = 0; i < 6; ++i ) {
            RS_left[i] = ropes[i];
            RS_right[i] = ropes[i];
        }

        // Get node children pointers for recursive traversal
        KdAccelNode *leftChild, *rightChild;
        leftChild = current + 1;
        rightChild = &nodes[current->AboveChild()];

        // Recurse
        RS_left[SR] = rightChild;
        attachRopes(leftChild, RS_left, true);

        RS_right[SL] = leftChild;
        attachRopes(rightChild, RS_right, true);
    }
}

// recursively attach ropes from root node in DF order
void KdTreeAccel::addParent(KdAccelNode *current, KdAccelNode *parent) {

    current->parent = parent;

    if ( !current->IsLeaf() ) {
        // Recurse
        addParent( current+1, current );
        addParent( &nodes[current->AboveChild()], current );
    }
}

void KdAccelNode::prettyPrint() const {
    std::cout << "bounding box min: ( " << bbox.min.x << ", " << bbox.min.y << ", " << bbox.min.z << " )" << std::endl;
    std::cout << "bounding box max: ( " << bbox.max.x << ", " << bbox.max.y << ", " << bbox.max.z << " )" << std::endl;
    if(IsLeaf()) {
        std::cout << "num_prim: " << nPrims << std::endl;
        if(nPrims > 1)
            std::cout << "first_prim_index: " << primitiveIndicesOffset << std::endl;
    }

    if(!IsLeaf()) {
        // Print split plane axis.
        if (SplitAxis() == X_AXIS) {
            std::cout << "split plane axis: X_AXIS" << std::endl;
        } else if (SplitAxis() == Y_AXIS) {
            std::cout << "split plane axis: Y_AXIS" << std::endl;
        } else if (SplitAxis() == Z_AXIS) {
            std::cout << "split plane axis: Z_AXIS" << std::endl;
        } else {
            std::cout << "split plane axis: invalid" << std::endl;
        }

        std::cout << "split plane value: " << SplitPos() << std::endl;
    }
    // Print whether or not node is a leaf node.
    if ( IsLeaf() ) {
        std::cout << "is leaf node: YES" << std::endl;
    }
    else {
        std::cout << "is leaf node: NO" << std::endl;
    }

    std::cout << "node index: " << nodeId << std::endl;

    // Print children indices.
    std::cout << "left child index: " << nodeId + 1 << std::endl;
    std::cout << "right child index: " << AboveChild() << std::endl;

    // Print neighboring nodes.
    for ( int i = 0; i < 6; ++i ) {
        if(ropes[i] != nullptr)
        std::cout << "neighbor node index " << i << ": " << ropes[i]->nodeId << std::endl;
    }

    // Print empty line.
    std::cout << std::endl;
}