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
 * BVH code based on pbrt-v3
 * Copyright (c) 1998-2015, Matt Pharr, Greg Humphreys, and Wenzel Jakob.
 */


#include <BoundingBox.h>
#include "BVH.h"
#include "Ray.h"
#include <cassert>
#include <Helper/ConsoleLogger.h>
#include <utility>
#include <random>
#include <Helper/Chronometer.h>

namespace Profiling {
    CummulativeBenchmark intersectStatsBVH{};
}

static Chronometer time_bvh;

namespace STATS {
    //STAT_MEMORY_COUNTER("Memory/BVH tree", treeBytes);
    //STAT_RATIO("BVH/Primitives per leaf node", totalPrimitives, totalLeafNodes);

    static int totalPrimitives = 0;
    static int totalLeafNodes = 0;
    static int interiorNodes = 0;
    static int leafNodes = 0;

    void _reset() {
        totalPrimitives = 0;
        totalLeafNodes = 0;
        interiorNodes = 0;
        leafNodes = 0;
    }
}

struct BVHBuildNode {
    // BVHBuildNode Public Methods
    void InitLeaf(int first, int n, const AxisAlignedBoundingBox &b) {
        firstPrimOffset = first;
        nPrimitives = n;
        bounds = b;
        children[0] = children[1] = nullptr;
        ++STATS::leafNodes;
        ++STATS::totalLeafNodes;
        STATS::totalPrimitives += n;
    }

    void InitInterior(int axis, BVHBuildNode *c0, BVHBuildNode *c1) {
        children[0] = c0;
        children[1] = c1;
        bounds = AxisAlignedBoundingBox::Union(c0->bounds, c1->bounds);
        splitAxis = axis;
        nPrimitives = 0;
        ++STATS::interiorNodes;
    }

    ~BVHBuildNode(){
        if(children[1]) delete children[1];
        if(children[0]) delete children[0];
    }
    AxisAlignedBoundingBox bounds;
    BVHBuildNode *children[2];
    int splitAxis, firstPrimOffset, nPrimitives;
};

/*struct LinearBVHNode {
    AxisAlignedBoundingBox bounds;
    union {
        int primitivesOffset;   // leaf
        int secondChildOffset;  // interior
    };
    uint16_t nPrimitives;  // 0 -> interior node
    uint8_t axis;          // interior node: xyz
    uint8_t pad[1];        // ensure 32 byte total size
};*/

int BVHAccel::SplitEqualCounts(std::vector<BVHPrimitiveInfo> &primitiveInfo, int start,
                               int end, int dim) {
    //<<Partition primitives into equally sized subsets>>
    int mid = (start + end) / 2;
    std::nth_element(&primitiveInfo[start], &primitiveInfo[mid],
                     &primitiveInfo[end - 1] + 1,
                     [dim](const BVHPrimitiveInfo &a, const BVHPrimitiveInfo &b) {
                         return a.centroid[dim] < b.centroid[dim];
                     });

    return mid;
}

int BVHAccel::SplitMiddle(std::vector<BVHPrimitiveInfo> &primitiveInfo, int start,
                          int end, int dim, AxisAlignedBoundingBox &centroidBounds) {
    //<<Partition primitives through node’s midpoint>>
    double pmid = (centroidBounds.min[dim] + centroidBounds.max[dim]) / 2;
    BVHPrimitiveInfo *midPtr =
            std::partition(&primitiveInfo[start], &primitiveInfo[end - 1] + 1,
                           [dim, pmid](const BVHPrimitiveInfo &pi) {
                               return pi.centroid[dim] < pmid;
                           });
    int mid = midPtr - &primitiveInfo[0];
    if (mid != start && mid != end)
        return mid;

    return SplitEqualCounts(primitiveInfo, start, end, dim);
}

int BVHAccel::SplitMiddleProb(std::vector<BVHPrimitiveInfo> &primitiveInfo, int start,
                              int end, int dim) {
    //<<Partition primitives through node’s midpoint>>
    std::sort(&primitiveInfo[start], &primitiveInfo[end],
              [](const BVHPrimitiveInfo &a, const BVHPrimitiveInfo &b) {
                  return a.probability < b.probability;
              });

    double sumProb = 0.0;
    for (int i = start; i < end; ++i) {
        sumProb += primitiveInfo[i].probability;
    }
    sumProb *= 0.5; // take half as split point

    double midProb = 0.0;
    int mid = start;
    for (; mid < end; mid += 2) {
        midProb += primitiveInfo[mid].probability;
        if (midProb > sumProb) break;
    }

    // If mid is valid, make swaps and return point
    if (mid != start && mid < end) {
        for (int i = start + 1, j = mid; i < mid && j < end; i += 2, j += 2) {
            std::swap(primitiveInfo[i], primitiveInfo[j]);
        }
        return mid;
    }

    return SplitEqualCounts(primitiveInfo, start, end, dim);
}

int BVHAccel::SplitProb(std::vector<BVHPrimitiveInfo> &primitiveInfo, int start,
                        int end, int dim, AxisAlignedBoundingBox &centroidBounds, AxisAlignedBoundingBox &bounds) {
    //<<Partition primitives using approximate SAH>>
    int nPrimitives = end - start;
    int mid = -1;
    if (nPrimitives <= 2) {
        //<<Partition primitives into equally sized subsets>>
        mid = SplitEqualCounts(primitiveInfo, start, end, dim);
    } else if (nPrimitives > maxPrimsInNode) { // else return -1 to just create a leaf node
        //<<Allocate BucketInfo for PROB partition buckets>>
        constexpr int nBuckets = 12;
        struct BucketInfo {
            int count = 0;
            double sumProb = 0.0;
            AxisAlignedBoundingBox bounds;
        };
        BucketInfo buckets[nBuckets];
        double maxProb = 0.0;
        double minProb = 1.0;
        //<<Initialize BucketInfo for PROB partition buckets>>
        for (int i = start; i < end; ++i) {
            int b = nBuckets *
                    centroidBounds.Offset(primitiveInfo[i].centroid)[dim];
            if (b == nBuckets) b = nBuckets - 1;
            assert(b >= 0);
            assert(b < nBuckets);
            buckets[b].count++;
            buckets[b].sumProb += primitiveInfo[i].probability;
            maxProb = std::max(maxProb, primitiveInfo[i].probability);
            minProb = std::min(minProb, primitiveInfo[i].probability);
            buckets[b].bounds = AxisAlignedBoundingBox::Union(buckets[b].bounds, primitiveInfo[i].bounds);
        }

        //<<Compute costs for splitting after each bucket>>
        double cost[nBuckets - 1];
        double cost_min[nBuckets - 1];
        double cost_mid[nBuckets - 1];
        double cost_sah[nBuckets - 1];

        double sumProb = 0.0;
        for (int i = 0; i < nBuckets - 1; ++i) {
            sumProb += buckets[i].sumProb;
        }

        int bestProbBucket = 0;
        int bestCount = 0;
        double prob = 0.0;
        for (int i = 0; i < nBuckets - 1; ++i) {
            prob += buckets[i].sumProb;
            bestCount += buckets[i].count;
            if (prob > sumProb * 0.5) {
                bestProbBucket = i;
                break;
            }
        }

        for (int i = 0; i < nBuckets - 1; ++i) {
            double prob0 = 0.0, prob1 = 0.0;
            AxisAlignedBoundingBox b0, b1;
            int count0 = 0, count1 = 0;
            for (int j = 0; j <= i; ++j) {
                count0 += buckets[j].count;
                prob0 += buckets[j].sumProb;
                b0 = AxisAlignedBoundingBox::Union(b0, buckets[j].bounds);
            }
            for (int j = i + 1; j < nBuckets; ++j) {
                count1 += buckets[j].count;
                prob1 += buckets[j].sumProb;
                b1 = AxisAlignedBoundingBox::Union(b1, buckets[j].bounds);
            }
            cost[i] = 0.5f + std::abs(prob0 - prob1) / maxProb;
            cost_min[i] = 0.5f + std::abs(prob0 - prob1) / minProb;
            cost_mid[i] = 0.5f + std::abs(prob0 - prob1) / (maxProb - minProb);
            cost_sah[i] = 0.5f + (count0 * b0.SurfaceArea() +
                                  count1 * b1.SurfaceArea()) / bounds.SurfaceArea();
        }

        //<<Find bucket to split at that minimizes PROB metric>>
        double minCost = cost_sah[0];
        int minCostSplitBucket = 0;
        for (int i = 1; i < nBuckets - 1; ++i) {
            if (cost_sah[i] < minCost) {
                minCost = cost_sah[i];
                minCostSplitBucket = i;
            }
        }
        //<<Either create leaf or split primitives at selected PROB bucket>>
        double leafCost = nPrimitives;
        if(sumProb <= 0.0) {
            if (minCost < leafCost) {
                BVHPrimitiveInfo *pmid = std::partition(&primitiveInfo[start],
                                                        &primitiveInfo[end - 1] + 1,
                                                        [=](const BVHPrimitiveInfo &pi) {
                                                            int b = nBuckets * centroidBounds.Offset(pi.centroid)[dim];
                                                            if (b == nBuckets) b = nBuckets - 1;
                                                            assert(b >= 0);
                                                            assert(b < nBuckets);
                                                            return b <= minCostSplitBucket;
                                                        });
                mid = pmid - &primitiveInfo[0];
            }
        }
        else{
            if (cost_sah[bestProbBucket] < leafCost) {
                BVHPrimitiveInfo *pmid = std::partition(&primitiveInfo[start],
                                                        &primitiveInfo[end - 1] + 1,
                                                        [=](const BVHPrimitiveInfo &pi) {
                                                            int b = nBuckets * centroidBounds.Offset(pi.centroid)[dim];
                                                            if (b == nBuckets) b = nBuckets - 1;
                                                            assert(b >= 0);
                                                            assert(b < nBuckets);
                                                            return b <= bestProbBucket;
                                                        });
                mid = pmid - &primitiveInfo[0];
            }
        }
    }
    return mid;
}

int BVHAccel::SplitSAH(std::vector<BVHPrimitiveInfo> &primitiveInfo, int start,
                       int end, int dim, AxisAlignedBoundingBox &centroidBounds, AxisAlignedBoundingBox &bounds) {
    //<<Partition primitives using approximate SAH>>
    int nPrimitives = end - start;
    int mid = -1;
    if (nPrimitives <= 2) {
        //<<Partition primitives into equally sized subsets>>
        mid = SplitEqualCounts(primitiveInfo, start, end, dim);
    } else if (nPrimitives > maxPrimsInNode) { // else return -1 to just create a leaf node
        //<<Allocate BucketInfo for SAH partition buckets>>
        constexpr int nBuckets = 12;
        struct BucketInfo {
            int count = 0;
            AxisAlignedBoundingBox bounds;
        };
        BucketInfo buckets[nBuckets];
        //<<Initialize BucketInfo for SAH partition buckets>>
        for (int i = start; i < end; ++i) {
            int b = nBuckets *
                    centroidBounds.Offset(primitiveInfo[i].centroid)[dim];
            if (b == nBuckets) b = nBuckets - 1;
            assert(b >= 0);
            assert(b < nBuckets);
            buckets[b].count++;
            buckets[b].bounds = AxisAlignedBoundingBox::Union(buckets[b].bounds, primitiveInfo[i].bounds);
        }

        //<<Compute costs for splitting after each bucket>>
        double cost[nBuckets - 1];
        for (int i = 0; i < nBuckets - 1; ++i) {
            AxisAlignedBoundingBox b0, b1;
            int count0 = 0, count1 = 0;
            for (int j = 0; j <= i; ++j) {
                b0 = AxisAlignedBoundingBox::Union(b0, buckets[j].bounds);
                count0 += buckets[j].count;
            }
            for (int j = i + 1; j < nBuckets; ++j) {
                b1 = AxisAlignedBoundingBox::Union(b1, buckets[j].bounds);
                count1 += buckets[j].count;
            }
            cost[i] = 0.5f + (count0 * b0.SurfaceArea() +
                              count1 * b1.SurfaceArea()) / bounds.SurfaceArea();
        }
        //<<Find bucket to split at that minimizes SAH metric>>
        double minCost = cost[0];
        int minCostSplitBucket = 0;
        for (int i = 1; i < nBuckets - 1; ++i) {
            if (cost[i] < minCost) {
                minCost = cost[i];
                minCostSplitBucket = i;
            }
        }
        //<<Either create leaf or split primitives at selected SAH bucket>>
        double leafCost = nPrimitives;
        if (minCost < leafCost) {
            BVHPrimitiveInfo *pmid = std::partition(&primitiveInfo[start],
                                                    &primitiveInfo[end - 1] + 1,
                                                    [=](const BVHPrimitiveInfo &pi) {
                                                        int b = nBuckets * centroidBounds.Offset(pi.centroid)[dim];
                                                        if (b == nBuckets) b = nBuckets - 1;
                                                        assert(b >= 0);
                                                        assert(b < nBuckets);
                                                        return b <= minCostSplitBucket;
                                                    });
            mid = pmid - &primitiveInfo[0];
        }
    }
    return mid;
}

int BVHAccel::SplitRDH(std::vector<BVHPrimitiveInfo> &primitiveInfo, int start,
                       int end, int dim, AxisAlignedBoundingBox &centroidBounds, std::vector<TestRay> &local_battery, AxisAlignedBoundingBox &bounds) {
    //<<Partition primitives using approximate SAH>>
    int nPrimitives = end - start;
    int mid = -1;
    if (nPrimitives <= 2) {
        //<<Partition primitives into equally sized subsets>>
        mid = SplitEqualCounts(primitiveInfo, start, end, dim);
    } else if (nPrimitives > maxPrimsInNode) { // else return -1 to just create a leaf node
        //<<Allocate BucketInfo for SAH partition buckets>>
        constexpr int nBuckets = 12;
        struct BucketInfo {
            int count = 0;
            AxisAlignedBoundingBox bounds;
        };
        BucketInfo buckets[nBuckets];
        //<<Initialize BucketInfo for SAH partition buckets>>
        for (int i = start; i < end; ++i) {
            int b = nBuckets *
                    centroidBounds.Offset(primitiveInfo[i].centroid)[dim];
            if (b == nBuckets) b = nBuckets - 1;
            assert(b >= 0);
            assert(b < nBuckets);
            buckets[b].count++;
            buckets[b].bounds = AxisAlignedBoundingBox::Union(buckets[b].bounds, primitiveInfo[i].bounds);
        }

        //<<Compute costs for splitting after each bucket>>
        double cost[nBuckets - 1];
        double hitcost[nBuckets - 1];

        std::vector<bool> rayHasHit(local_battery.size(), false);
        for (int i = 0; i < nBuckets - 1; ++i) {
            AxisAlignedBoundingBox b0, b1;
            int count0 = 0, count1 = 0;
            for (int j = 0; j <= i; ++j) {
                b0 = AxisAlignedBoundingBox::Union(b0, buckets[j].bounds);
                count0 += buckets[j].count;
            }
            for (int j = i + 1; j < nBuckets; ++j) {
                b1 = AxisAlignedBoundingBox::Union(b1, buckets[j].bounds);
                count1 += buckets[j].count;
            }

            int hitCount1 = 0;
            int hitCount0 = 0;
            int hitCountBoth = 0;
            auto ray = Ray();
#pragma omp parallel default(none) firstprivate(ray) shared(rayHasHit, hitCount0, hitCount1, b0, b1, hitCountBoth, local_battery)
            {

                int hitCount0_loc = 0;
                int hitCount1_loc = 0;
                int hitCountBoth_loc = 0;

#pragma omp for
                for (int sample_id = 0; sample_id < local_battery.size(); sample_id++) {
                    auto &test = local_battery[sample_id];
                    ray.origin = test.pos;
                    ray.direction = test.dir;
                    Vector3d invDir(1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z);
                    int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};

                    bool hit0 = b0.IntersectBox(ray, invDir, dirIsNeg) ? true : false;
                    bool hit1 = b1.IntersectBox(ray, invDir, dirIsNeg) ? true : false;
                    if (hit0 && hit1)
                        hitCountBoth_loc++;
                    else {
                        if (hit0) {
                            hitCount0_loc++;
                        }
                        if (hit1) {
                            hitCount1_loc++;
                        }
                    }
                    rayHasHit[sample_id] = rayHasHit[sample_id] | (hit0 || hit1);
                }

                // parallel reduce
#pragma omp critical
                {
                    hitCount0 += hitCount0_loc;
                    hitCount1 += hitCount1_loc;
                    hitCountBoth += hitCountBoth_loc;
                }
            }
            cost[i] = 0.5f + (count0 * b0.SurfaceArea() +
                    count1 * b1.SurfaceArea()) / bounds.SurfaceArea();
            //hitcost[i] = 0.5f + (double)(hitcount0 + hitcount1) / (double)local_battery.size();
            double pAbove = (double) (hitCount0) / (double) local_battery.size();
            double pBelow = (double) (hitCount1) / (double) local_battery.size();
            double pBoth = (double) (hitCountBoth) / (double) local_battery.size(); // as a penalty

            /*hitcost[i] = 1.0f * std::exp(-2.30f * (std::pow(pAbove,2) + std::pow(pBelow,2)));
            hitcost[i] *= nPrimitives;
            hitcost[i] -= nPrimitives;
            hitcost[i] = std::abs(hitcost[i]);*/
            //hitcost[i] = 0.5f + 1.0f - std::abs(pAbove - pBelow) + pBoth;
            hitcost[i] = 0.5f + 1.0f *
                            (pBelow * count0 + pAbove * count1);

        }

        //<<Find bucket to split at that minimizes SAH metric>>
        double minCost = hitcost[0];
        int minCostSplitBucket = 0;
        for (int i = 1; i < nBuckets - 1; ++i) {
            if (hitcost[i] < minCost) {
                minCost = hitcost[i];
                minCostSplitBucket = i;
            }
        }

        /*size_t ray_n = 0;
        for(auto iter = local_battery.begin(); iter != local_battery.end(); ){
            if(!rayHasHit[ray_n]) iter = local_battery.erase(iter);
            else iter++;
            ray_n++;
        }*/

        //<<Either create leaf or split primitives at selected SAH bucket>>
        double leafCost = nPrimitives;
        if (minCost < leafCost) {
            BVHPrimitiveInfo *pmid = std::partition(&primitiveInfo[start],
                                                    &primitiveInfo[end - 1] + 1,
                                                    [=](const BVHPrimitiveInfo &pi) {
                int b = nBuckets * centroidBounds.Offset(pi.centroid)[dim];
                if (b == nBuckets) b = nBuckets - 1;
                assert(b >= 0);
                assert(b < nBuckets);
                return b <= minCostSplitBucket;
            });
            mid = pmid - &primitiveInfo[0];
        }
    }
    return mid;
}

BVHAccel::BVHAccel(const std::vector<TestRay> &battery, std::vector<std::shared_ptr<Primitive>> p, int maxPrimsInNode,
                   SplitMethod splitMethod)
                   : maxPrimsInNode(std::min(255, maxPrimsInNode)),
                   splitMethod(splitMethod),
                   primitives(std::move(p)),
                   battery(battery) {

    /*std::unique_ptr<double[]> primChance(new double[primitives.size()]);
    for (size_t i = 0; i < primitives.size(); ++i) {
        primChance[i] = 0.0;
    }

    auto ray = Ray();
    ray.rng = new MersenneTwister();
    for (auto &test : this->battery) {
        ray.origin = test.pos;
        ray.direction = test.dir;
        Vector3d invDir(1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z);
        int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};
        for (size_t i = 0; i < primitives.size(); ++i) {
            primChance[i] += primitives[i]->Intersect(ray);
        }
    }
    delete ray.rng;*/

    // if battery is too large, only use a random sample
    /*if(this->battery.size() > HITCACHELIMIT / 128) {
        std::vector<TestRay> sampleRays;
        std::sample(this->battery.begin(), this->battery.end(), std::back_inserter(sampleRays),
                    HITCACHELIMIT / 128, std::mt19937{std::random_device{}()});
        this->battery.swap(sampleRays);
    }*/

    /*for (size_t i = 0; i < primitives.size(); ++i) {
        primChance[i] /= this->battery.size();
    }*/

    std::vector<BVHPrimitiveInfo> primitiveInfo(0);
    if (splitMethod == SplitMethod::RDH) {
        if (this->battery.empty())
            this->splitMethod = SplitMethod::SAH;
    }

    Log::console_msg_master(4, "--- BVH with HitBattery ---\n");
    Log::console_msg_master(4, " Battery size   : {}\n", battery.size());
    Log::console_msg_master(4, " Sampled battery: {}\n", this->battery.size());

    construct(primitiveInfo);
}

BVHAccel::BVHAccel(const std::vector<double> &probabilities, std::vector<std::shared_ptr<Primitive>> p,
                   int maxPrimsInNode,
                   SplitMethod splitMethod)
        : maxPrimsInNode(std::min(255, maxPrimsInNode)),
          splitMethod(splitMethod),
          primitives(std::move(p)) {

    std::vector<BVHPrimitiveInfo> primitiveInfo;

    if (splitMethod == SplitMethod::ProbSplit) {
        if (primitives.size() != probabilities.size())
            this->splitMethod = SplitMethod::SAH;
    }

    if (splitMethod == SplitMethod::ProbSplit) {
        primitiveInfo.resize(primitives.size());
        for (size_t i = 0; i < primitives.size(); ++i)
            primitiveInfo[i] = {i, primitives[i]->sh.bb, probabilities[primitives[i]->globalId]};
        construct(primitiveInfo);
    }
    else{
        construct(primitiveInfo);
    }

}

BVHAccel::~BVHAccel() {
    if (nodes) {
        delete[] nodes;
        nodes = nullptr;
    }
};

BVHAccel::BVHAccel(std::vector<std::shared_ptr<Primitive>> p,
                   int maxPrimsInNode, SplitMethod splitMethod)
                   : maxPrimsInNode(std::min(255, maxPrimsInNode)),
                   splitMethod(splitMethod),
                   primitives(std::move(p)) {

    construct(std::vector<BVHPrimitiveInfo>());
}

void BVHAccel::construct(std::vector<BVHPrimitiveInfo> primitiveInfo){
    nodes = nullptr;
    if (primitives.empty())
        return;
    STATS::_reset();

    //<<Build BVH from primitives>>
    //<<Initialize primitiveInfo array for primitives>>

    if(primitiveInfo.empty()) {
        std::vector<BVHPrimitiveInfo>(primitives.size()).swap(primitiveInfo);
        if (splitMethod != SplitMethod::ProbSplit) {
            for (size_t i = 0; i < primitives.size(); ++i)
                primitiveInfo[i] = {i, primitives[i]->sh.bb};
        } else {
            return;
        }
    }

    //<<Build BVH tree for primitives using primitiveInfo>>
    //MemoryArena arena(1024 * 1024);
    int totalNodes = 0;
    std::vector<std::shared_ptr<Primitive>> orderedPrims;
    BVHBuildNode *root;
    /*if (splitMethod == SplitMethod::HLBVH)
        root = HLBVHBuild(arena, primitiveInfo, &totalNodes, orderedPrims);
    else
        */

    std::vector<TestRayLoc> indices;
    indices.reserve(battery.size());
    for (int i = 0; i < battery.size(); ++i)
        indices.emplace_back(i, 0, 1.0e99);

    root = recursiveBuild(/*arena,*/ primitiveInfo, 0, primitives.size(),
                                     &totalNodes, orderedPrims, battery, indices);
    primitives.swap(orderedPrims);

    Log::console_msg_master(4, "BVH created with {} nodes for {} "
           "primitives ({:.2f} MB)\n",
           totalNodes, (int) primitives.size(),
           float(totalNodes * sizeof(LinearBVHNode)) /
           (1024.f * 1024.f));
    //<<Compute representation of depth-first traversal of BVH tree>>
    nodes = new LinearBVHNode[totalNodes]; //AllocAligned<LinearBVHNode>(totalNodes);
    int offset = 0;

    std::vector<IntersectCount>(totalNodes).swap(ints);
    flattenBVHTreeStats(root, &offset, 0);
    assert(totalNodes == offset);
    delete root;

    Log::console_msg_master(4,"--- BVH STATS ---\n");
    Log::console_msg_master(4," Total Primitives: {}\n", STATS::totalPrimitives);
    Log::console_msg_master(4," Total Leaf Nodes: {}\n", STATS::totalLeafNodes);
    Log::console_msg_master(4," Interior Nodes:   {}\n", STATS::interiorNodes);
    Log::console_msg_master(4," Leaf Nodes:       {}\n", STATS::leafNodes);
};

BVHBuildNode *BVHAccel::recursiveBuild(std::vector<BVHPrimitiveInfo> &primitiveInfo, int start, int end,
                                       int *totalNodes,
                                       std::vector<std::shared_ptr<Primitive>> &orderedPrims,
                                       std::vector<TestRay> &testBattery,
                                       const std::vector<TestRayLoc> &local_battery) {
    assert(start != end);

    BVHBuildNode *node = new BVHBuildNode();//arena.Alloc<BVHBuildNode>();
    (*totalNodes)++;
    //<<Compute bounds of all primitives in BVH node>>
    AxisAlignedBoundingBox bounds;
    for (int i = start; i < end; ++i)
        bounds = AxisAlignedBoundingBox::Union(bounds, primitiveInfo[i].bounds);

    int nPrimitives = end - start;
    if (nPrimitives == 1) {
        //<<Create leaf BVHBuildNode>>
        int firstPrimOffset = orderedPrims.size();
        for (int i = start; i < end; ++i) {
            int primNum = primitiveInfo[i].primitiveNumber;
            orderedPrims.push_back(primitives[primNum]);
        }
        node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
        return node;
    } else {
        //<<Compute bound of primitive centroids, choose split dimension dim>>
        AxisAlignedBoundingBox centroidBounds;
        for (int i = start; i < end; ++i)
            centroidBounds = AxisAlignedBoundingBox::Union(centroidBounds, primitiveInfo[i].centroid);
        int dim = centroidBounds.MaximumExtent();

        //<<Partition primitives into two sets and build children>>
        int mid = (start + end) / 2;
        if (centroidBounds.max[dim] == centroidBounds.min[dim]) {
            //<<Create leaf BVHBuildNode>>
            int firstPrimOffset = orderedPrims.size();
            for (int i = start; i < end; ++i) {
                int primNum = primitiveInfo[i].primitiveNumber;
                orderedPrims.push_back(primitives[primNum]);
            }
            node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
            return node;
        } else {
            //<<Partition primitives based on splitMethod>>
            switch (splitMethod) {
                case SplitMethod::MolflowSplit: {
                    // Mix of Middle and EqualCounts
                    // Find best middle cut (most equal counts) in all dimensions
                    int tmpMid = 0;
                    double bestDeviation = 1.0e99;
                    for (auto itDim : {0, 1, 2}) {
                        tmpMid = SplitMiddle(primitiveInfo, start, end, itDim, centroidBounds);
                        double dev = std::abs((double) tmpMid - (double) ((end - start) / 2.0 + start));
                        if (dev < bestDeviation) {
                            bestDeviation = dev;
                            mid = tmpMid;
                            dim = itDim;
                        }
                    }
                    break;
                }
                case SplitMethod::Middle: {
                    mid = SplitMiddle(primitiveInfo, start, end, dim, centroidBounds);
                    break;
                }
                case SplitMethod::EqualCounts: {
                    //<<Partition primitives into equally sized subsets>>
                    mid = SplitEqualCounts(primitiveInfo, start, end, dim);
                    break;
                }
                case SplitMethod::ProbSplit: {
                    //<<Partition primitives into equally sized subsets>>
                    mid = SplitProb(primitiveInfo, start, end, dim, centroidBounds, bounds);
                    //mid = SplitMiddleProb(primitiveInfo, start, end, dim);
                    break;
                }
                case SplitMethod::RDH: {
                    if(local_battery.size() > HITCACHEMIN) {
                        //local_battery = testBattery;
                        //<<Partition primitives into equally sized subsets>>
                        mid = SplitRDH(primitiveInfo, start, end, dim, centroidBounds,
                                       const_cast<std::vector<TestRay> &>(testBattery), bounds);
                        //mid = SplitMiddleProb(primitiveInfo, start, end, dim);
                        break;
                    }
                }
                case SplitMethod::SAH:
                default: {
                    mid = SplitSAH(primitiveInfo, start, end, dim, centroidBounds, bounds);
                    break;
                }
            }

            if (mid < 0 || nPrimitives <= maxPrimsInNode) {
                // Create leaf node, when max prims are reached
                int firstPrimOffset = orderedPrims.size();
                for (int i = start; i < end; ++i) {
                    int primNum = primitiveInfo[i].primitiveNumber;
                    orderedPrims.push_back(primitives[primNum]);
                }
                node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
                return node;
            }

            bool skipL = false;
            bool skipR = false;
            std::vector<TestRayLoc> battery_below, battery_above;
            if(!skipL || !skipR)
            {
                AxisAlignedBoundingBox b0{bounds};
                AxisAlignedBoundingBox b1{bounds};
                b0.max[dim] = b1.min[dim] = mid;
                auto ray = Ray();
#pragma omp parallel default(none) firstprivate(ray) shared(local_battery, battery_above, battery_below, b0, b1)
                {
                    std::vector<TestRayLoc> battery_below_local, battery_above_local;
#pragma omp for
                    for (int sample_id = 0; sample_id < local_battery.size(); sample_id++) {
                        auto& ind = local_battery[sample_id];
                        ray.origin = battery[ind.index].pos;
                        ray.direction = battery[ind.index].dir;
                        Vector3d invDir(1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z);
                        int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};
                        bool hit0 = b0.IntersectBox(ray, invDir, dirIsNeg) ? true : false;
                        bool hit1 = b1.IntersectBox(ray, invDir, dirIsNeg) ? true : false;

                        {

                            // Advance to next child node, possibly enqueue other child
                            {
                                if (hit0) {
                                    battery_above_local.emplace_back(ind.index, ind.tMin, ind.tMax);
                                }
                                if (hit1) {
                                    battery_below_local.emplace_back(ind.index, ind.tMin, ind.tMax);
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

            node->InitInterior(dim,
                               recursiveBuild(primitiveInfo, start, mid,
                                              totalNodes, orderedPrims, testBattery, battery_above),
                               recursiveBuild(primitiveInfo, mid, end,
                                              totalNodes, orderedPrims, testBattery, battery_below));
        }
    }
    return node;
}

int BVHAccel::flattenBVHTree(BVHBuildNode *node, int *offset) {
    LinearBVHNode *linearNode = &nodes[*offset];
    linearNode->bounds = node->bounds;
    int myOffset = (*offset)++;
    if (node->nPrimitives > 0) {
        assert(!node->children[0] && !node->children[1]);
        assert(node->nPrimitives < 65536);
        linearNode->primitivesOffset = node->firstPrimOffset;
        linearNode->nPrimitives = node->nPrimitives;
    } else {
        // Create interior flattened BVH node
        linearNode->axis = node->splitAxis;
        linearNode->nPrimitives = 0;
        flattenBVHTree(node->children[0], offset);
        linearNode->secondChildOffset =
                flattenBVHTree(node->children[1], offset);
    }
    return myOffset;
}

int BVHAccel::flattenBVHTreeStats(BVHBuildNode *node, int *offset, int level) {
    LinearBVHNode *linearNode = &nodes[*offset];
    linearNode->bounds = node->bounds;

    ints[*offset].nbPrim = node->nPrimitives;
    ints[*offset].level = level;

    int myOffset = (*offset)++;
    if (node->nPrimitives > 0) {
        assert(!node->children[0] && !node->children[1]);
        assert(node->nPrimitives < 65536);
        linearNode->primitivesOffset = node->firstPrimOffset;
        linearNode->nPrimitives = node->nPrimitives;
    } else {
        // Create interior flattened BVH node
        linearNode->axis = node->splitAxis;
        linearNode->nPrimitives = 0;
        flattenBVHTreeStats(node->children[0], offset, level+1);
        linearNode->secondChildOffset =
                flattenBVHTreeStats(node->children[1], offset, level+1);
    }
    return myOffset;
}

void BVHAccel::ComputeBB() {
    bb = nodes ? nodes[0].bounds : AxisAlignedBoundingBox();
}

bool BVHAccel::Intersect(Ray &ray) {
    if (!nodes) return false;

    bool hit = false;
    Vector3d invDir(1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z);
    int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};
    // Follow ray through BVH nodes to find primitive intersections
    int toVisitOffset = 0, currentNodeIndex = 0;
    int nodesToVisit[64];
    while (true) {
        const LinearBVHNode *node = &nodes[currentNodeIndex];
        // Check ray against BVH node
        if (node->bounds.IntersectBox(ray, invDir, dirIsNeg)) {
            if (node->nPrimitives > 0) {
                // Intersect ray with primitives in leaf BVH node
                for (int i = 0; i < node->nPrimitives; ++i) {

                    const std::shared_ptr<Primitive> &p = primitives[node->primitivesOffset + i];
                    // Do not check last collided facet to prevent self intersections
                    if (p->globalId != ray.lastIntersected && p->Intersect(ray)) {
                        hit = true;
                    }
                }
                if (toVisitOffset == 0) break;
                currentNodeIndex = nodesToVisit[--toVisitOffset];
            } else {
                // Put far BVH node on _nodesToVisit_ stack, advance to near
                // node
                if (dirIsNeg[node->axis]) {
                    nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
                    currentNodeIndex = node->secondChildOffset;
                } else {
                    nodesToVisit[toVisitOffset++] = node->secondChildOffset;
                    currentNodeIndex = currentNodeIndex + 1;
                }
            }
        } else {
            if (toVisitOffset == 0) break;
            currentNodeIndex = nodesToVisit[--toVisitOffset];
        }
    }

    return hit;
}

bool BVHAccel::IntersectStat(RayStat &ray) {
    if (!nodes) return false;

    /*Chronometer time_outerloop; // 0
Chronometer time_intersectbox; // 1
Chronometer time_intersect; // 2
if(algorithm_times.size() != 3)
    std::vector<CummulativeBenchmark>(3,CummulativeBenchmark()).swap(algorithm_times);
*/
    time_bvh.ReStart();

    bool hit = false;
    Vector3d invDir(1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z);
    int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};
    // Follow ray through BVH nodes to find primitive intersections
    int toVisitOffset = 0, currentNodeIndex = 0;
    int nodesToVisit[64];
    while (true) {
        //time_outerloop.ReStart();
        const LinearBVHNode *node = &nodes[currentNodeIndex];
        // Check ray against BVH node
        ray.stats.traversalSteps++;

        // nbChecks count box and primitive intersection tests
        ints[currentNodeIndex].nbChecks += ray.stats.traversalSteps;
        if (node->bounds.IntersectBox(ray, invDir, dirIsNeg)) {
            //time_intersectbox.ReStart();
            ints[currentNodeIndex].nbIntersects++;
            if (node->nPrimitives > 0) {
                // Intersect ray with primitives in leaf BVH node
                for (int i = 0; i < node->nPrimitives; ++i) {

                    const std::shared_ptr<Primitive> &p = primitives[node->primitivesOffset + i];
                    // Do not check last collided facet to prevent self intersections
                    ray.stats.traversalSteps++;
                    if (p->globalId != ray.lastIntersected && p->IntersectStat(ray)) {
                        hit = true;
                    }
                }
                if (toVisitOffset == 0) break;
                currentNodeIndex = nodesToVisit[--toVisitOffset];
            } else {
                // Put far BVH node on _nodesToVisit_ stack, advance to near
                // node
                if (dirIsNeg[node->axis]) {
                    nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
                    currentNodeIndex = node->secondChildOffset;
                } else {
                    nodesToVisit[toVisitOffset++] = node->secondChildOffset;
                    currentNodeIndex = currentNodeIndex + 1;
                }
            }
            //this->algorithm_times[1].AddTime(time_intersectbox.ElapsedMs());
        } else {
            if (toVisitOffset == 0) break;
            currentNodeIndex = nodesToVisit[--toVisitOffset];
        }

        //this->algorithm_times[0].AddTime(time_outerloop.ElapsedMs());
    }

    perRayCount += ray.stats;
    ray.ResetStats();
    Profiling::intersectStatsBVH.AddTime(time_bvh.ElapsedMs());

    return hit;
}

BVHAccel::BVHAccel(BVHAccel &&src) noexcept: maxPrimsInNode(src.maxPrimsInNode), splitMethod(src.splitMethod) {
    primitives = std::move(src.primitives);
    nodes = src.nodes;
    src.nodes = nullptr;
    bb = src.bb;
}

BVHAccel::BVHAccel(const BVHAccel &src) noexcept: maxPrimsInNode(src.maxPrimsInNode), splitMethod(src.splitMethod) {
    if (nodes)
        exit(44);
}

BVHAccel &BVHAccel::operator=(const BVHAccel &src) noexcept {
    if (nodes)
        exit(44);

    return *this;
}

std::ostream &operator<<(std::ostream &os, BVHAccel::SplitMethod split_type) {
    switch (split_type) {
        case BVHAccel::SplitMethod::SAH : return os << "SAH" ;
        case BVHAccel::SplitMethod::HLBVH: return os << "HLBVH";
        case BVHAccel::SplitMethod::Middle: return os << "Middle";
        case BVHAccel::SplitMethod::EqualCounts: return os << "EqualCounts";
        case BVHAccel::SplitMethod::MolflowSplit: return os << "MolflowSplit";
        case BVHAccel::SplitMethod::ProbSplit: return os << "ProbSplit";
        case BVHAccel::SplitMethod::RDH: return os << "RDH";
        default: return os << "Unknown";
            // omit default case to trigger compiler warning for missing cases
    };
    return os << static_cast<std::uint16_t>(split_type);
}