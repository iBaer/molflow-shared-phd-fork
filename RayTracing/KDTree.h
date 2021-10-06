//
// Created by pascal on 4/21/21.
//

#ifndef MOLFLOW_PROJ_KDTREE_H
#define MOLFLOW_PROJ_KDTREE_H

#include <vector>
#include <memory>
#include <FacetData.h>
#include "Primitive.h"

using Primitive = Facet;

// KdTreeAccel Forward Declarations
// KdTreeAccel Local Declarations
struct KdAccelNode {
    // KdAccelNode Methods
    void InitLeaf(int *primNums, int np, std::vector<int> *primitiveIndices);
    void InitInterior(int axis, int ac, double s);
    double SplitPos() const { return split; }
    int nPrimitives() const { return nPrims >> 2; }
    int SplitAxis() const { return flags & 3; }
    bool IsLeaf() const { return (flags & 3) == 3; }
    int AboveChild() const { return aboveChild >> 2; }
    union {
        double split;                 // Interior
        int onePrimitive;            // Leaf
        int primitiveIndicesOffset;  // Leaf
    };

    // stats
    int nodeId{0};

private:
    union {
        int flags;       // Both
        int nPrims;      // Leaf
        int aboveChild;  // Interior
    };

    friend class Geometry;
};

struct BoundEdge;

struct RTStats {
    size_t nti{0}; // traversed interior nodes
    size_t ntl{0}; // traversed leaves
    size_t nit{0}; // ray-bject intersection
    double timeTrav{0};
    double timeInt{0};
};

struct TestRayLoc {
    TestRayLoc(size_t ind, double min, double max) {
        index = ind;
        tMin = min;
        tMax = max;
    }
    size_t index{0};
    double tMin{0.0};
    double tMax{0.0};
};

struct SplitCandidate {
    SplitCandidate() = default;
    SplitCandidate(double c, int a, int o){
        cost=c;
        axis=a;
        offset=o;
    }
    double cost{1.0e99};
    int axis{-1};
    int offset{-1};
};

class KdTreeAccel : public RTPrimitive {
public:
    using BoundaryFloat = double;
    struct RaySegment{
        BoundaryFloat tMin{-1.0e99};
        BoundaryFloat tMax{1.0e99};
        const TestRay* ray{nullptr};
        std::vector<Vector2d> distances;
    };
    struct RayBoundary{
        RaySegment* rs{nullptr};
        int flag{-1}; // 0: left, 1: on splitting plane, 2: right
    };
    // local stack shows which _RayBoundary/ies from the whole set of boundaries belong to
    //the left, to the right, and to both of them
    struct BoundaryStack{
        BoundaryStack(int ll, int lh, int rl, int rh) :
            leftMin(ll), leftMax(lh), rightMin(rl), rightMax(rh) {};
        int leftMin{-1};
        int leftMax{-1};
        int rightMin{-1};
        int rightMax{-1};
    };

    // KdTreeAccel Public Types
    enum class SplitMethod {
        SAH, ProbSplit, TestSplit, HybridSplit, HybridBin
    };

public:
    // KdTreeAccel Public Methods
    /*KdTreeAccel(SplitMethod splitMethod, std::vector<std::shared_ptr<Primitive>> p,
                const std::vector<double> &probabilities = std::vector<double>{},
                int isectCost = 80, int traversalCost = 1,
                double emptyBonus = 0.5, int maxPrims = 1, int maxDepth = -1);
    */
    KdTreeAccel(SplitMethod splitMethod, std::vector<std::shared_ptr<Primitive>> p,
                const std::vector<double> &probabilities = std::vector<double>{}, const std::vector<TestRay> &battery = std::vector<TestRay>{},
                int isectCost = 80, int traversalCost = 1,
                double emptyBonus = 0.5, int maxPrims = 1, int maxDepth = -1);
    KdTreeAccel(KdTreeAccel && src) noexcept;
    KdTreeAccel(const KdTreeAccel & src) noexcept;

    KdTreeAccel& operator=(const KdTreeAccel & src) noexcept;
    ~KdTreeAccel() override;

    bool Intersect(Ray &ray) override;
    RTStats IntersectT(Ray &ray);

    std::vector<IntersectCount> ints;
private:
    void ComputeBB() override;
    // KdTreeAccel Private Methods
    void buildTree(int nodeNum, const AxisAlignedBoundingBox &nodeBounds,
                   const std::vector<AxisAlignedBoundingBox> &allPrimBounds, int *primNums, int nPrimitives,
                   int depth, const std::unique_ptr<BoundEdge[]> edges[3], int *prims0, int *prims1,
                   int badRefines, const std::vector<double> &probabilities, int prevSplitAxis);
    void buildTree(int nodeNum, const AxisAlignedBoundingBox &nodeBounds,
                   const std::vector<AxisAlignedBoundingBox> &allPrimBounds, int *primNums, int nPrimitives,
                   int depth, const std::unique_ptr<BoundEdge[]> edges[3], int *prims0, int *prims1,
                   int badRefines, const std::vector<TestRay> &battery,
                   const std::vector<TestRayLoc> &local_battery, int prevSplitAxis, double tMin, double tMax,
                   const std::vector<double> &primChance, double oldCost);
    void buildTreeRDH(int nodeNum, const AxisAlignedBoundingBox &nodeBounds,
                      const std::vector<AxisAlignedBoundingBox> &allPrimBounds, int *primNums, int nPrimitives,
                      int depth, const std::unique_ptr<BoundEdge[]> edges[3], int *prims0, int *prims1,
                      int badRefines, const std::unique_ptr<std::vector<RaySegment>> &battery,
                      RayBoundary ** rb_stack, int prevSplitAxis, double tMin,
                      double tMax, const std::vector<double> &primChance, double oldCost, BoundaryStack bound);
private:
    // KdTreeAccel Private Data
    SplitMethod splitMethod;
    const int isectCost, traversalCost, maxPrims;
    const double emptyBonus;
    std::vector<std::shared_ptr<Primitive>> primitives;
    std::vector<int> primitiveIndices;
    KdAccelNode *nodes;
    int nAllocedNodes, nextFreeNode;
    AxisAlignedBoundingBox bounds;

    std::tuple<double, int, int> SplitSAH(int axis, const AxisAlignedBoundingBox &nodeBounds,
                            const std::vector<AxisAlignedBoundingBox> &allPrimBounds, int *primNums,
                            int nPrimitives, const std::unique_ptr<BoundEdge[]> edges[3]);
    std::tuple<double, int, int> SplitProb(int axis, const AxisAlignedBoundingBox &nodeBounds,
                            const std::vector<AxisAlignedBoundingBox> &allPrimBounds, int *primNums, int nPrimitives,
                            const std::unique_ptr<BoundEdge[]> edges[3],
                            const std::vector<double> &probabilities);

    std::tuple<double, int, int> SplitHybrid(int axis, const AxisAlignedBoundingBox &nodeBounds,
                                             const std::vector<AxisAlignedBoundingBox> &allPrimBounds,
                                             int *primNums, int nPrimitives,
                                             const std::unique_ptr<BoundEdge[]> edges[3],
                                             const std::vector<TestRay> &battery,
                                             const std::vector<TestRayLoc> &local_battery,
                                             const std::vector<double> &primChance, double tMax) const;
    std::tuple<double, int, int, bool, bool> SplitHybridBin(int axis, const AxisAlignedBoundingBox &nodeBounds,
                                                            const std::vector<AxisAlignedBoundingBox> &allPrimBounds,
                                                            int *primNums, int nPrimitives,
                                                            const std::unique_ptr<BoundEdge[]> edges[3],
                                                            const std::vector<TestRay> &battery,
                                                            const std::vector<TestRayLoc> &local_battery,
                                                            const std::vector<double> &primChance, double tMax) const;
    std::tuple<double, int, int>
    SplitTest(int axis, const AxisAlignedBoundingBox &nodeBounds,
              const std::vector<AxisAlignedBoundingBox> &allPrimBounds,
              int *primNums, int nPrimitives,
              const std::unique_ptr<BoundEdge[]> edges[3],
              const std::vector<TestRay> &battery, const std::vector<TestRayLoc> &local_battery,
              const std::vector<double> &primChance, double tMax);

              std::tuple<double, int, int>
    SplitTest(int axis, const AxisAlignedBoundingBox &nodeBounds,
              const std::vector<AxisAlignedBoundingBox> &allPrimBounds,
              int *primNums, int nPrimitives,
              const std::unique_ptr<BoundEdge[]> edges[3],
              const std::unique_ptr<std::vector<RaySegment>> &battery,
              RayBoundary **rb_stack,
              const std::vector<double> &primChance, double tMax,
              const BoundaryStack &boundary_stack);

    friend class Geometry;

    void PrintTreeInfo();

    bool IntersectStat(RayStat &ray);
};

struct KdToDo {
    const KdAccelNode *node;
    double tMin, tMax;
};

#endif //MOLFLOW_PROJ_KDTREE_H