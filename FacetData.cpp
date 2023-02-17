//
// Created by pbahr on 08/02/2021.
//

#include "FacetData.h"
#include "Polygon.h"
#include "Helper/MathTools.h"
#include "RayTracing/RTHelper.h" // SubProcessFacetTempVar
#include "RayTracing/Ray.h" // hitlink

bool Facet::Intersect(Ray &ray) {
    //++iSCount;
    Vector3d rayDirOpposite(-1.0 * ray.direction);
    double det = Dot(this->sh.Nuv, rayDirOpposite);

    // Eliminate "back facet"
    if ((this->sh.is2sided) || (det > 0.0)) { //If 2-sided or if ray going opposite facet normal

        double u, v, d;
        // Ray/rectangle instersection. Find (u,v,dist) and check 0<=u<=1, 0<=v<=1, dist>=0

        if (det != 0.0) {

            double iDet = 1.0 / det;
            Vector3d intZ = ray(0) - this->sh.O;

            u = iDet * DET33(intZ.x, this->sh.V.x, rayDirOpposite.x,
                             intZ.y, this->sh.V.y, rayDirOpposite.y,
                             intZ.z, this->sh.V.z, rayDirOpposite.z);

            if (u >= 0.0 && u <= 1.0) {

                v = iDet * DET33(this->sh.U.x, intZ.x, rayDirOpposite.x,
                                 this->sh.U.y, intZ.y, rayDirOpposite.y,
                                 this->sh.U.z, intZ.z, rayDirOpposite.z);

                if (v >= 0.0 && v <= 1.0) {

                    d = iDet * Dot(this->sh.Nuv, intZ);

                    if (d>0.0) {

                        // Now check intersection with the facet polygon (in the u,v space)
                        // This check could be avoided on rectangular facet.
                        if (IsInPoly(u, v, vertices2)) {
                            bool hardHit = this->surf->IsHardHit(ray);
/*#if defined(SYNRAD)

                            if(typeid(*surf) == typeid(MaterialSurface))
                            hardHit &= (this->sh.reflectType > 10 //Material reflection
                                        && currentParticle.model->materials[this->sh.reflectType - 10].hasBackscattering //Has complex scattering
                                        && currentParticle.model->materials[this->sh.reflectType - 10].GetReflectionType(currentParticle.energy,
                                                                                                                         acos(Dot(currentParticle.direction, this->sh.N)) - PI / 2, currentParticle.randomGenerator.rnd()) == REFL_TRANS));

#endif*/
                            /*#if defined(MOLFLOW)
                            double time = ray.time + d / 100.0 / currentParticle.velocity;
                            double currentOpacity = currentParticle.model->GetOpacityAt(f, time);
                            hardHit = ((currentOpacity == 1.0) || (currentParticle.randomGenerator.rnd()<currentOpacity));
#endif

#if defined(SYNRAD)
                            hardHit = !((this->sh.opacity < 0.999999 //Partially transparent facet
										&& currentParticle.randomGenerator.rnd()>this->sh.opacity)
										|| (this->sh.reflectType > 10 //Material reflection
										&& currentParticle.model->materials[this->sh.reflectType - 10].hasBackscattering //Has complex scattering
										&& currentParticle.model->materials[this->sh.reflectType - 10].GetReflectionType(currentParticle.energy,
										acos(Dot(currentParticle.direction, this->sh.N)) - PI / 2, currentParticle.randomGenerator.rnd()) == REFL_TRANS));
#endif*/
                            if (hardHit) {

                                // Hard hit
                                if (d < ray.tMax) {
                                    ray.tMax = d;
                                    //ray.lastIntersected = this->globalId;
                                    /*if(ray.hitChain->hit){
                                        ray.hitChain->next = new HitChain();
                                        ray.hitChain = ray.hitChain->next;
                                    }
                                    ray.hitChain->hit = new SubProcessFacetTempVar();
                                    ray.hitChain->hit->isHit = true;
                                    ray.hitChain->hit->colU = u;
                                    ray.hitChain->hit->colV = v;
                                    ray.hitChain->hit->colDistTranspPass = d;
                                    ray.hitChain->hitId = globalId;*/

                                    ray.hits.emplace_back(globalId, SubProcessFacetTempVar());
                                    auto& hit = ray.hits.back().hit;
                                    hit.isHit = true;
                                    hit.colU = u;
                                    hit.colV = v;
                                    hit.colDistTranspPass = d;

                                }
                            }
                            else {
                                /*if(ray.hitChain->hit){
                                    ray.hitChain->next = new HitChain();
                                    ray.hitChain = ray.hitChain->next;
                                }
                                ray.hitChain->hit = new SubProcessFacetTempVar();
                                ray.hitChain->hit->isHit = false;
                                ray.hitChain->hit->colU = u;
                                ray.hitChain->hit->colV = v;
                                ray.hitChain->hit->colDistTranspPass = d;
                                ray.hitChain->hitId = globalId;*/

                                ray.hits.emplace_back(globalId, SubProcessFacetTempVar());
                                auto& hit = ray.hits.back().hit;
                                hit.isHit = false;
                                hit.colU = u;
                                hit.colV = v;
                                hit.colDistTranspPass = d;
                            }
                            return hardHit;
                        } // IsInFacet
                    } // d range
                } // u range
            } // v range
        } // det==0
    } // dot<0

    return false;
}

inline double FMA(double a, double b, double c) {
    return a * b + c;
}

inline auto DifferenceOfProducts(double a, double b, double c, double d) {
    auto cd = c * d;
    auto differenceOfProducts = FMA(a, b, -cd);
    auto error = FMA(-c, d, cd);
    return differenceOfProducts + error;
}

constexpr double machEps =
        std::numeric_limits<double>::epsilon() * 0.5;

constexpr double gamma(int n)
{
    return (n * machEps) / (1 - n * machEps);
}

// Transform barycentrics to UV coordinates via texture coordinates
static inline
Vector2d getHitLocation(const Vector2d barycentrics, Vector2d* texCoord, unsigned int texIndex)
{
    const double u = barycentrics.u; // beta and gamma
    const double v = barycentrics.v; // and gamma
    const double w = 1.0f - u - v; // and alpha

    Vector2d tex = w * texCoord[texIndex]
            +         u * texCoord[texIndex+1]
            +         v * texCoord[texIndex+2];

    return tex; //hitLocationU , hitLocationV
}

// Transform barycentrics to UV coordinates via texture coordinates
static inline
Vector2d getHitLocation(double u, double v, double w, Vector2d* texCoord, unsigned int texIndex)
{
    /*const double u = barycentrics.u; // beta and gamma
    const double v = barycentrics.v; // and gamma
    const double w = 1.0f - u - v; // and alpha*/

    Vector2d tex = w * texCoord[texIndex]
            +         u * texCoord[texIndex+1]
            +         v * texCoord[texIndex+2];

    return tex; //hitLocationU , hitLocationV
}

// Transformation to UV coordinates via Cramer's rule
// without direction info
static inline
Vector2d getHitLocation_old(const FacetProperties& poly, const Vector3d& rayOrigin)
{
    const Vector3d b = rayOrigin - poly.O;

    double det = poly.U.x * poly.V.y - poly.U.y * poly.V.x; // Could be precalculated
    double detU = b.x * poly.V.y - b.y * poly.V.x;
    double detV = poly.U.x * b.y - poly.U.y * b.x;

    if(abs(det)<=1.0e-12f){
        det = poly.U.y * poly.V.z - poly.U.z * poly.V.y;
        detU = b.y * poly.V.z - b.z * poly.V.y;
        detV = poly.U.y * b.z - poly.U.z * b.y;
        if(abs(det)<=1.0e-12f){
            det = poly.U.z * poly.V.x - poly.U.x * poly.V.z;
            detU = b.z * poly.V.x - b.x * poly.V.z;
            detV = poly.U.z * b.x - poly.U.x * b.z;
#ifdef DEBUG
            if(fabsf(det)<=std::numeric_limits<float>::epsilon()){
                printf("[HitLoc] Dangerous determinant calculated: %e : %e : %e -> %e : %e\n",det,detU,detV,detU/det,detV/det);
            }
#endif
        }
    }

    return Vector2d(detU/det, detV/det); //hitLocationU , hitLocationV
}

bool TriangleFacet::Intersect(Ray &ray) {
    auto& vert = *this->vertices3;
    auto p0 = vert[indices[0]], p1 = vert[indices[1]], p2 = vert[indices[2]];

    // Return no intersection if triangle is degenerate
    {
        auto cross = CrossProduct(p2 - p0, p1 - p0);
        if (Dot(cross, cross) == 0)
            return {};
    }

    // Transform triangle vertices to ray coordinate space
    // Translate vertices based on ray origin
    auto p0t = p0 - ray.origin;
    auto p1t = p1 - ray.origin;
    auto p2t = p2 - ray.origin;

    // Permute components of triangle vertices and ray direction
    int kz = MaxComponentIndex(Abs(ray.direction));
    int kx = kz + 1;
    if (kx == 3)
        kx = 0;
    int ky = kx + 1;
    if (ky == 3)
        ky = 0;
    auto d = Permute(ray.direction, {kx, ky, kz});
    p0t = Permute(p0t, {kx, ky, kz});
    p1t = Permute(p1t, {kx, ky, kz});
    p2t = Permute(p2t, {kx, ky, kz});

    // Apply shear transformation to translated vertex positions
    double Sx = -d.x / d.z;
    double Sy = -d.y / d.z;
    double Sz = 1.0 / d.z;
    p0t.x += Sx * p0t.z;
    p0t.y += Sy * p0t.z;
    p1t.x += Sx * p1t.z;
    p1t.y += Sy * p1t.z;
    p2t.x += Sx * p2t.z;
    p2t.y += Sy * p2t.z;

    // Compute edge function coefficients _e0_, _e1_, and _e2_
    double e0 = DifferenceOfProducts(p1t.x, p2t.y, p1t.y, p2t.x);
    double e1 = DifferenceOfProducts(p2t.x, p0t.y, p2t.y, p0t.x);
    double e2 = DifferenceOfProducts(p0t.x, p1t.y, p0t.y, p1t.x);

    // Perform triangle edge and determinant tests
    if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
        return {};
    double det = e0 + e1 + e2;
    if (det == 0)
        return {};
    
    // Compute scaled hit distance to triangle and test against ray $t$ range
    p0t.z *= Sz;
    p1t.z *= Sz;
    p2t.z *= Sz;
    double tScaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
    if (det < 0 && (tScaled >= 0 || tScaled < ray.tMax * det))
        return {};
    else if (det > 0 && (tScaled <= 0 || tScaled > ray.tMax * det))
        return {};
    
    // Compute barycentric coordinates and $t$ value for triangle intersection
    double invDet = 1.0 / det;
    double b0 = e0 * invDet, b1 = e1 * invDet, b2 = e2 * invDet;
    double t = tScaled * invDet;
    //DCHECK(!IsNaN(t));

    // Ensure that computed triangle $t$ is conservatively greater than zero
    // Compute $\delta_z$ term for triangle $t$ error bounds
    double maxZt = MaxComponentValue(Abs(Vector3d(p0t.z, p1t.z, p2t.z)));
    double deltaZ = gamma(3) * maxZt;

    // Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
    double maxXt = MaxComponentValue(Abs(Vector3d(p0t.x, p1t.x, p2t.x)));
    double maxYt = MaxComponentValue(Abs(Vector3d(p0t.y, p1t.y, p2t.y)));
    double deltaX = gamma(5) * (maxXt + maxZt);
    double deltaY = gamma(5) * (maxYt + maxZt);

    // Compute $\delta_e$ term for triangle $t$ error bounds
    double deltaE = 2 * (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);

    // Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
    double maxE = MaxComponentValue(Abs(Vector3d(e0, e1, e2)));
    double deltaT =
            3 * (gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) * std::abs(invDet);
    if (t <= deltaT)
        return {};

    // Return _TriangleIntersection_ for intersection
    //return TriangleIntersection{b0, b1, b2, t};
    
    {
        bool hardHit = this->surf->IsHardHit(ray);
        /*#if defined(SYNRAD)
        
                                    if(typeid(*surf) == typeid(MaterialSurface))
                                    hardHit &= (this->sh.reflectType > 10 //Material reflection
                                                && currentParticle.model->materials[this->sh.reflectType - 10].hasBackscattering //Has complex scattering
                                                && currentParticle.model->materials[this->sh.reflectType - 10].GetReflectionType(currentParticle.energy,
                                                                                                                                 acos(Dot(currentParticle.direction, this->sh.N)) - PI / 2, currentParticle.randomGeneratoray.rnd()) == REFL_TRANS));
        
        #endif*/
        /*#if defined(MOLFLOW)
        double time = ray.time + d / 100.0 / currentParticle.velocity;
        double currentOpacity = currentParticle.model->GetOpacityAt(f, time);
        hardHit = ((currentOpacity == 1.0) || (currentParticle.randomGeneratoray.rnd()<currentOpacity));
#endif

#if defined(SYNRAD)
        hardHit = !((this->sh.opacity < 0.999999 //Partially transparent facet
                    && currentParticle.randomGeneratoray.rnd()>this->sh.opacity)
                    || (this->sh.reflectType > 10 //Material reflection
                    && currentParticle.model->materials[this->sh.reflectType - 10].hasBackscattering //Has complex scattering
                    && currentParticle.model->materials[this->sh.reflectType - 10].GetReflectionType(currentParticle.energy,
                    acos(Dot(currentParticle.direction, this->sh.N)) - PI / 2, currentParticle.randomGeneratoray.rnd()) == REFL_TRANS));
#endif*/
        if (hardHit) {

            // Hard hit
            if (t > deltaT) {
                ray.tMax = t;

                ray.hits.emplace_back(globalId, SubProcessFacetTempVar());
                auto& hit = ray.hits.back().hit;
                hit.isHit = true;
                auto coord = getHitLocation(b0,b1,b2,this->texCoord, 0);
                //auto coord = getHitLocation_old(sh, ray.origin);
                hit.colU = coord.u;
                hit.colV = coord.v;
                hit.colDistTranspPass = t;

            }
        }
        else {

            ray.hits.emplace_back(globalId, SubProcessFacetTempVar());
            auto& hit = ray.hits.back().hit;
            hit.isHit = false;
            auto coord = getHitLocation(b0,b1,b2,this->texCoord, 0);
            //auto coord = getHitLocation_old(sh, ray.origin);
            hit.colU = coord.u;
            hit.colV = coord.v;
            hit.colDistTranspPass = t;
        }
        return hardHit;
    }
    return true;
}
