﻿/*
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
#if defined(MOLFLOW)
#include "../src/MolFlow.h"
#endif

#if defined(SYNRAD)
#include "../src/SynRad.h"
#endif
#include "Facet_shared.h"
#include "Polygon.h"
//#include <malloc.h>

#include <string.h>
#include <math.h>
#include "GLApp/GLToolkit.h"
#include "Helper/MathTools.h"
#include <sstream>

using namespace pugi;

#if defined(MOLFLOW)
extern MolFlow *mApp;
#endif

#if defined(SYNRAD)
extern SynRad*mApp;
#endif

/**
* \brief Constructor with initialisation based on the number of indices/facets
* \param nbIndex number of indices/facets
*/
InterfaceFacet::InterfaceFacet(size_t nbIndex) : sh(0) {
	indices.resize(nbIndex);                    // Ref to Geometry Vector3d
	vertices2.resize(nbIndex);
	visible.resize(nbIndex);

	memset(&facetHitCache, 0, sizeof(FacetHitBuffer));

	sh.nbIndex = nbIndex;

	sh.sticking = 0.0;
	sh.opacity = 1.0;

	sh.profileType = PROFILE_NONE;
	
	sh.texWidth = 0;
	sh.texHeight = 0;
	sh.texWidth_precise = 0.0;
	sh.texHeight_precise = 0.0;
	sh.center.x = 0.0;
	sh.center.y = 0.0;
	sh.center.z = 0.0;

	sh.is2sided = false;
	sh.isProfile = false;
	//wp.isOpaque = true;
	sh.isTextured = false;
	sh.countAbs = false;
	sh.countRefl = false;
	sh.countTrans = false;
	sh.countDirection = false;

	sh.superIdx = 0;
	sh.superDest = 0;
	sh.teleportDest = 0;
	sh.isVolatile = false;

	textureVisible = true;
	volumeVisible = true;

	texDimW = 0;
	texDimH = 0;
    tRatioU = 0.0;
    tRatioV = 0.0;

	//mesh = NULL;
	//meshPts = NULL;
	cellPropertiesIds.clear();
	meshvector.clear();
	meshvectorsize = 0;
	hasMesh = false;
	//nbElem = 0;
	selectedElem.u = 0;
	selectedElem.v = 0;
	selectedElem.width = 0;
	selectedElem.height = 0;
	dirCache = nullptr;
	textureError = false;

	glTex = 0;
	glList = 0;
	glElem = 0;
	glSelElem = 0;
	selected = false;

#if defined(MOLFLOW)
	

	sh.temperature = 293.15; // 20degC
	sh.outgassing = 0.0;           // 1 unit*l/s //will be outgasssing
	sh.desorbType = DES_NONE;
	sh.desorbTypeN = 0.0;

	sh.reflection.diffusePart = 1.0; //totally diffuse reflection
	sh.reflection.specularPart = 0.0;
	sh.reflection.cosineExponent = 0.0; //Cos^0 = uniform

	sh.countDes = false;
	sh.countACD = false;
	sh.useOutgassingFile = false;
	sh.accomodationFactor = 1.0;

	sh.enableSojournTime = false;
	sh.sojournFreq = 1E13;
	sh.sojournE = 100;

	sh.outgassing_paramId = -1;
	sh.opacity_paramId = -1;
	sh.sticking_paramId = -1;

	sh.isMoving = false;

	hasOutgassingFile = false;
	//outgassingMapWindow = NULL;

	sh.anglemapParams.record = false;

	sh.anglemapParams.phiWidth = sh.anglemapParams.thetaLowerRes = sh.anglemapParams.thetaHigherRes = 0;
	sh.anglemapParams.thetaLimit = 1.570796326; //slightly lower than PI/2

	angleMapCache = std::vector<size_t>(); //SAFE_DELETE called on it, must initialize

	//sh.facetHistogramParams.record = false;

    ogMap.totalFlux = sh.totalOutgassing = ogMap.totalDose = 0.0;

	userOutgassing = "";
	userOpacity = "";
	userSticking = "";
#endif

#if defined(SYNRAD)
	sh.doScattering = false;
	sh.rmsRoughness = 100.0E-9; //100nm
	sh.autoCorrLength = 100 * 100E-9; //tau=autoCorr/RMS=100

	sh.reflectType = REFLECTION_SPECULAR;
	sh.recordSpectrum = false;
#endif
}

/**
* \brief Destructor for safe deletion
*/
InterfaceFacet::~InterfaceFacet() {
	  //SAFE_DELETE(cellPropertiesIds);
	  SAFE_FREE(dirCache);
	  DELETE_TEX(glTex);
	  DELETE_LIST(glList);
	  DELETE_LIST(glElem);
	  DELETE_LIST(glSelElem);
	  /*for (size_t i = 0; i < meshvectorsize; i++)
          SAFE_DELETE(meshvector[i].points);*/
	  //SAFE_DELETE(meshvector);
#if defined(MOLFLOW)
	  //SAFE_FREE(outgassingMapWindow);
	  //SAFE_FREE(angleMapCache);
#endif
}

/*
void Facet::DetectOrientation() {

	// Detect polygon orientation (clockwise or counter clockwise)
	// p= 1.0 => The second vertex is convex and vertex are counter clockwise.
	// p=-1.0 => The second vertex is concave and vertex are clockwise.
	// p= 0.0 => The polygon is not a simple one and orientation cannot be detected.

	GLAppPolygon p;
	p.pts=vertices2;
	p.sign = 1.0;

	bool convexFound = false;
	size_t i = 0;
	while (i < p.pts.size() && !convexFound) {
		
		auto [empty,center] = EmptyTriangle(p, (int)i - 1, (int)i, (int)i + 1);
		if (empty || sh.nbIndex == 3) {
			size_t _i1 = Previous(i, p.pts.size());
			size_t _i2 = IDX(i, p.pts.size());
			size_t _i3 = Next(i, p.pts.size());
			if (IsInPoly(center, p.pts)) {
				convexFound = true;
				// Orientation
				if (IsConvex(p, i)) p.sign = 1.0;
				else                p.sign = -1.0;
			}
		}
		i++;
	}

	if (!convexFound) {
		// Not a simple polygon
		sh.sign = 0.0;
	}
	else {
		sh.sign = p.sign;
	}
}
*/

/**
* \brief Restore texture and geometry information
* \return if restoration was okay
*/
int InterfaceFacet::RestoreDeviceObjects() {

	// Initialize scene objects (OpenGL)
	if (sh.isTextured) {
		glGenTextures(1, &glTex);
		glList = glGenLists(1);
	}

	//BuildMeshGLList();
	BuildSelElemList();

	return GL_OK;

}

/**
* \brief Invalidate texture and geometry information
* \return if invalidation was okay
*/
int InterfaceFacet::InvalidateDeviceObjects() {

	// Free all alocated resource (OpenGL)
	DELETE_TEX(glTex);
	DELETE_LIST(glList);
	DELETE_LIST(glElem);
	DELETE_LIST(glSelElem);
	return GL_OK;

}

/**
* \brief Set texture on facet
* \param width width of the texture
* \param height height of the texture
* \param useMesh true if a new mesh needs to be created (if none exists f->hasMesh)
* \return true if texture was set
*/
bool InterfaceFacet::SetTextureProperties(double width, double height, bool useMesh) {

	bool dimOK = (width*height > 0.0000001);

	if (dimOK) {
        const double ceilCutoff = 0.9999999;
        sh.texWidth_precise = width;
		sh.texHeight_precise = height;
		sh.texWidth = (int)ceil(width * ceilCutoff); //0.9999999: cut the last few digits (convert rounding error 1.00000001 to 1, not 2)
		sh.texHeight = (int)ceil(height * ceilCutoff);
		dimOK = (sh.texWidth > 0 && sh.texHeight > 0);
	}
	else {
		sh.texWidth = 0;
		sh.texHeight = 0;
		sh.texWidth_precise = 0.0;
		sh.texHeight_precise = 0.0;
	}

	UpdateFlags(); //set hasMesh to true if everything was OK
	return true;

}

/**
* \brief Set texture on facet
* \param width width of the texture
* \param height height of the texture
* \param useMesh true if a new mesh needs to be created (if none exists f->hasMesh)
* \return true if texture was set
*/
bool InterfaceFacet::SetTexture(double width, double height, bool useMesh) {

	bool dimOK = (width*height > 0.0000001);

	texDimW = 0;
	texDimH = 0;
	hasMesh = false;
	//SAFE_FREE(mesh);
	/*for (size_t i = 0; i < meshvectorsize; i++)
        SAFE_DELETE(meshvector[i].points);*/
	//SAFE_DELETE(meshvector);
	meshvectorsize = 0;
	SAFE_FREE(dirCache);
	DELETE_TEX(glTex);
	DELETE_LIST(glList);
	DELETE_LIST(glElem);
	/*if (meshPts) {
		for (size_t i = 0; i < nbElem; i++)
			SAFE_FREE(meshPts[i].pts);
	}*/

	//SAFE_FREE(meshPts);
	//SAFE_FREE(cellPropertiesIds);
	//nbElem = 0;
	UnselectElem();

	if (dimOK) {

		// Add a 1 texel border for bilinear filtering (rendering purpose)
		texDimW = GetPower2(sh.texWidth + 2);
		texDimH = GetPower2(sh.texHeight + 2);
		if (texDimW < 4) texDimW = 4;
		if (texDimH < 4) texDimH = 4;
		glGenTextures(1, &glTex);
		glList = glGenLists(1);
		if (useMesh)
			if (!BuildMesh()) return false;
		if (sh.countDirection) {
			dirCache = (DirectionCell *)calloc(sh.texWidth*sh.texHeight, sizeof(DirectionCell));
			if (!dirCache) return false;
			//memset(dirCache,0,dirSize); //already done by calloc
		}

	}

    UpdateFlags(); //set hasMesh to true if everything was OK
    return true;

}

/**
* \brief Converts (u,v) Vertex to 3D Vertex
* \param u local u coordinate of facet
* \param v local v coordinate of facet
*/
void InterfaceFacet::glVertex2u(double u, double v) {

	glVertex3d(sh.O.x + sh.U.x*u + sh.V.x*v,
		sh.O.y + sh.U.y*u + sh.V.y*v,
		sh.O.z + sh.U.z*u + sh.V.z*v);

}

/**
* \brief Allocate memory for mesh and initialise
* \return true if mesh properly build
*/
bool InterfaceFacet::BuildMesh() {
	try{
        cellPropertiesIds.resize(sh.texWidth * sh.texHeight, 0);
	    meshvector.resize(sh.texWidth * sh.texHeight, CellProperties()); //will shrink at the end
	}
    catch (const std::exception &e) {
		std::cerr << "Couldn't allocate memory" << std::endl;
		return false;
		//throw Error("malloc failed on Facet::BuildMesh()");
	}
	//memset(meshvector, 0, sh.texWidth * sh.texHeight * sizeof(CellProperties));
	
	meshvectorsize = 0;
	hasMesh = true;
	
	GLAppPolygon P1, P2;
	double sx, sy;
	double iw = 1.0 / (double)sh.texWidth_precise;
	double ih = 1.0 / (double)sh.texHeight_precise;
	double rw = sh.U.Norme() * iw;
	double rh = sh.V.Norme() * ih;
	double fullCellArea = iw*ih;

	std::vector<Vector2d>(4).swap(P1.pts);
	//P1.sign = 1;
	P2.pts = vertices2; 
	//P2.sign = -sign;

	for (size_t j = 0;j < sh.texHeight;j++) {
		sy = (double)j;
		for (size_t i = 0;i < sh.texWidth;i++) {
			sx = (double)i;

			bool allInside = false;
			double u0 = sx * iw;
			double v0 = sy * ih;
			double u1 = (sx + 1.0) * iw;
			double v1 = (sy + 1.0) * ih;
			//mesh[i + j*wp.texWidth].elemId = -1;

			if (sh.nbIndex <= 4) {

				// Optimization for quad and triangle
                allInside = IsInPoly(Vector2d(u0,v0), vertices2)
                            && IsInPoly(Vector2d(u0, v1), vertices2)
                            && IsInPoly(Vector2d(u1, v0), vertices2)
                            && IsInPoly(Vector2d(u1, v1), vertices2);

			}

			if (!allInside) {
				CellProperties cellprop{};

				// Intersect element with the facet (facet boundaries)
				P1.pts[0].u = u0;
				P1.pts[0].v = v0;
				P1.pts[1].u = u1;
				P1.pts[1].v = v0;
				P1.pts[2].u = u1;
				P1.pts[2].v = v1;
				P1.pts[3].u = u0;
				P1.pts[3].v = v1;
				auto [A,center,vList] = GetInterArea(P1, P2, visible);
				if (!IsZero(A)) {

					if (A > (fullCellArea + 1e-10)) {

						// Polyon intersection error !
						// Switch back to brute force
						auto [bfArea,center] = GetInterAreaBF(P2, Vector2d(u0, v0), Vector2d(u1, v1));
						bool fullElem = IsZero(fullCellArea - bfArea);
						if (!fullElem) {
							cellprop.area = (bfArea*(rw*rh) / (iw*ih));
							cellprop.uCenter = (float)center.u;
							cellprop.vCenter = (float)center.v;
							cellprop.nbPoints = 0;
							cellprop.points.clear();
							cellPropertiesIds[i + j*sh.texWidth] = (int)meshvectorsize;
							meshvector[meshvectorsize++] = cellprop;
						}
						else {
							cellPropertiesIds[i + j*sh.texWidth] = -1;
						}

						//cellprop.full = IsZero(fullCellArea - A);

					}
					else {

						bool fullElem = IsZero(fullCellArea - A);
						if (!fullElem) {
							// !! P1 and P2 are in u,v coordinates !!
							cellprop.area = (A*(rw*rh) / (iw*ih));
							cellprop.uCenter = (float)center.u;
							cellprop.vCenter = (float)center.v;
							//cellprop.full = IsZero(fullCellArea - A);
							//cellprop.elemId = nbElem;

							// Mesh coordinates
							cellprop.points = //(Vector2d*)malloc(sizeof(vList[0])*vList.size());
							    vList;
							//memcpy(cellprop.points, vList.data(), sizeof(vList[0])*vList.size());
							cellprop.nbPoints = vList.size();
							cellPropertiesIds[i + j*sh.texWidth] = (int)meshvectorsize;
							meshvector[meshvectorsize++] = cellprop;
							//nbElem++;

						}
						else {
							cellPropertiesIds[i + j*sh.texWidth] = -1;
						}

					}

				}
				else cellPropertiesIds[i + j*sh.texWidth] = -2; //zero element

			}
			else {  //All indide and triangle or quad
				cellPropertiesIds[i + j*sh.texWidth] = -1;

				/*mesh[i + j*wp.texWidth].area = (float)(rw*rh);
				mesh[i + j*wp.texWidth].uCenter = (float)(u0 + u1) / 2.0f;
				mesh[i + j*wp.texWidth].vCenter = (float)(v0 + v1) / 2.0f;
				mesh[i + j*wp.texWidth].full = true;
				mesh[i + j*wp.texWidth].elemId = nbElem;

				// Mesh coordinates
				meshPts[nbElem].nbPts = 4;
				meshPts[nbElem].pts = (Vector2d *)malloc(4 * sizeof(Vector2d));

				if (!meshPts[nbElem].pts) {
					throw Error("Couldn't allocate memory for texture mesh points.");
				}
				meshPts[nbElem].pts[0].u = u0;
				meshPts[nbElem].pts[0].v = v0;
				meshPts[nbElem].pts[1].u = u1;
				meshPts[nbElem].pts[1].v = v0;
				meshPts[nbElem].pts[2].u = u1;
				meshPts[nbElem].pts[2].v = v1;
				meshPts[nbElem].pts[3].u = u0;
				meshPts[nbElem].pts[3].v = v1;
				nbElem++;*/

			}

			//tA += mesh[i + j*wp.texWidth].area;

		}
	}
	//Shrink mesh vector
	//meshvector = (CellProperties*)realloc(meshvector, sizeof(CellProperties)*meshvectorsize);
    meshvector.resize(meshvectorsize);

	// Check meshing accuracy (TODO)
	/*
	int p = (int)(ceil(log10(wp.area)));
	double delta = pow(10.0,(double)(p-5));
	if( fabs(wp.area - tA)>delta ) {
	}
	*/

	if (mApp->needsMesh) BuildMeshGLList();
	return true;

}

/**
* \brief Build OpenGL geometry for meshing
*/
void InterfaceFacet::BuildMeshGLList() {

	if (cellPropertiesIds.empty())
		return;

	DELETE_LIST(glElem);

	// Build OpenGL geometry for meshing
	glElem = glGenLists(1);
	glNewList(glElem, GL_COMPILE);

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	size_t nb = sh.texWidth*sh.texHeight;
	for (size_t i = 0; i < nb; i++) {
		if (cellPropertiesIds[i] != -2) {
			glBegin(GL_POLYGON);
			size_t nbPts = GetMeshNbPoint(i);
			/*size_t nbDrawn = 0;
			size_t n;
			if (mApp->leftHandedView) {
				n = 0;
			}
			else {
				n = nbPts - 1;
			}
			for (; nbDrawn < nbPts; nbDrawn++) {*/
			for (size_t n=0;n<nbPts;n++) {
				glEdgeFlag(true);
				Vector2d pt = GetMeshPoint(i, n);
				glVertex2u(pt.u, pt.v);
				/*if (mApp->leftHandedView) {
					n++;
				}
				else {
					n--;
				}*/
			}
			glEnd();
		}

	}

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEndList();

}

/**
* \brief Build GL List for selected elements
*/
void InterfaceFacet::BuildSelElemList() {

	DELETE_LIST(glSelElem);
	int nbSel = 0;

	if (!cellPropertiesIds.empty() && selectedElem.width != 0 && selectedElem.height != 0) {

		glSelElem = glGenLists(1);
		glNewList(glSelElem, GL_COMPILE);
		glColor3f(1.0f, 1.0f, 1.0f);
		glLineWidth(1.0f);
		glEnable(GL_LINE_SMOOTH);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glEnable(GL_POLYGON_OFFSET_LINE);
		glPolygonOffset(-1.0f, -1.0f);
		for (size_t i = 0; i < selectedElem.width; i++) {
			for (size_t j = 0; j < selectedElem.height; j++) {
				size_t add = (selectedElem.u + i) + (selectedElem.v + j)*sh.texWidth;
				//int elId = mesh[add].elemId;

				//if (cellPropertiesIds[add]!=-1 && meshvector[cellPropertiesIds[add]].elemId>=0) {
				if (cellPropertiesIds[add] != -2) {

					glBegin(GL_POLYGON);
					/*for (int n = 0; n < meshPts[elId].nbPts; n++) {
						glEdgeFlag(true);
						glVertex2u(meshPts[elId].pts[n].u, meshPts[elId].pts[n].v);
					}*/
					for (size_t p = 0;p < GetMeshNbPoint(add);p++) {
						Vector2d point = GetMeshPoint(add, p);
						glEdgeFlag(true);
						glVertex2u(point.u, point.v);
					}
					glEnd();
					nbSel++;
				}

			}
		}

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDisable(GL_POLYGON_OFFSET_LINE);
		glDisable(GL_LINE_SMOOTH);
		glEndList();

		// Empty selection
		if (nbSel == 0) UnselectElem();
	}
}

/**
* \brief Clear elements from selected List
*/
void InterfaceFacet::UnselectElem() {

	DELETE_LIST(glSelElem);
	selectedElem.width = 0;
	selectedElem.height = 0;

}

/**
* \brief Clear selection and select new element/s (pixel from a texture) based on local coordinates and rectangle size
* \param u local u coordinate of element
* \param v local v coordinate of element
* \param width width of element
* \param height height of element
*/
void InterfaceFacet::SelectElem(size_t u, size_t v, size_t width, size_t height) {

	UnselectElem();

	if (!cellPropertiesIds.empty() && u >= 0 && u < sh.texWidth && v >= 0 && v < sh.texHeight) {

		size_t maxW = sh.texWidth - u;
		size_t maxH = sh.texHeight - v;
		selectedElem.u = u;
		selectedElem.v = v;
		selectedElem.width = Min(maxW, width);
		selectedElem.height = Min(maxH, height);
		BuildSelElemList();

	}

}

/**
* \brief For specific rendering a selected element
*/
void InterfaceFacet::RenderSelectedElem() {
	if (glSelElem) glCallList(glSelElem);
}

/**
* \brief Fill vertex array
* \param v Array of vertices
*/
void InterfaceFacet::FillVertexArray(InterfaceVertex *v) {

	int nb = 0;
	for (size_t i = 0;i < sh.texHeight*sh.texWidth;i++) {
		if (cellPropertiesIds[i] != -2) {
			for (size_t j = 0; j < GetMeshNbPoint(i); j++) {
				Vector2d p = GetMeshPoint(i, j);
				v[nb].x = sh.O.x + sh.U.x*p.u + sh.V.x*p.v;
				v[nb].y = sh.O.y + sh.U.y*p.u + sh.V.y*p.v;
				v[nb].z = sh.O.z + sh.U.z*p.u + sh.V.z*p.v;
				nb++;
			}
		}
	}

}

/**
* \brief Get Texture Swap size
* \param useColormap if a colormap is used or not
* \return texture swap size
*/
size_t InterfaceFacet::GetTexSwapSize(bool useColormap) {

	size_t tSize = texDimW*texDimH;
	if (useColormap) tSize = tSize * 4;
	return tSize;

}

/**
* \brief Get Texture Swap size calculated with a size ratio
* \param ratio ratio used for size conversion
* \param useColor if a colormap is used or not
* \return texture swap size
*/
size_t InterfaceFacet::GetTexSwapSizeForRatio(double ratio, bool useColor) {

	double nU = sh.U.Norme();
	double nV = sh.V.Norme();
	double width = nU*ratio;
	double height = nV*ratio;

	bool dimOK = (width*height > 0.0000001);

	if (dimOK) {

		int iwidth = (int)ceil(width);
		int iheight = (int)ceil(height);
		int m = Max(iwidth, iheight);
		size_t tDim = GetPower2(m);
		if (tDim < 16) tDim = 16;
		size_t tSize = tDim*tDim;
		if (useColor) tSize *= 4;
		return tSize;

	}
	else {
		return 0;
	}
}

/**
* \brief Get number of texture cells
* \return pair of number of texture cells in both directions
*/
std::pair<size_t, size_t> InterfaceFacet::GetNbCell() {
	return std::make_pair(sh.texHeight, sh.texWidth);
}

/**
* \brief Get number of cells calculated with a size ratio
* \param ratio ratio used for size conversion
* \return number of texture cells
*/
size_t InterfaceFacet::GetNbCellForRatio(double ratio) {

	double nU = sh.U.Norme();
	double nV = sh.V.Norme();
	double width = nU*ratio;
	double height = nV*ratio;

	bool dimOK = (width*height > 0.0000001);

	if (dimOK) {
        int iWidth = (int)ceil(width);
        int iHeight = (int)ceil(height);
        return iWidth * iHeight;
	}
    return 0;
}

/**
* \brief Get number of cells calculated with a size ratio
* \param ratioU ratio in U direction used for size conversion
* \param ratioV ratio in V direction used for size conversion
* \return number of texture cells
*/
std::pair<size_t, size_t> InterfaceFacet::GetNbCellForRatio(double ratioU, double ratioV) {

    double nU = sh.U.Norme();
    double nV = sh.V.Norme();
    double width = nU*ratioU;
    double height = nV*ratioV;

    bool dimOK = (width*height > 0.0000001);
    if (dimOK) {
        const double ceilCutoff = 0.9999999;
        int iWidth = (int)ceil(width * ceilCutoff); //0.9999999: cut the last few digits (convert rounding error 1.00000001 to 1, not 2)
        int iHeight = (int)ceil(height * ceilCutoff);
        return std::make_pair(iWidth, iHeight);
    }
    else {
        return std::make_pair(0, 0);
    }
}

/**
* \brief Revert vertex order
*/
void InterfaceFacet::SwapNormal() {

	// Revert vertex order (around the second point)

	/*
	size_t* tmp = (size_t *)malloc(sh.nbIndex * sizeof(size_t));
	for (size_t i = sh.nbIndex, j = 0; i > 0; i--, j++) //Underrun-safe
		tmp[(i + 1) % sh.nbIndex] = GetIndex((int)j + 1);
	free(indices);
	indices = tmp;
	*/

	std::reverse(indices.begin(), indices.end());

	/* normal recalculated at reinitialize
	// Invert normal
	wp.N.x = -wp.N.x;
	wp.N.y = -wp.N.y;
	wp.N.z = -wp.N.z;*/

}

/**
* \brief Shift vertex order by offset to the left
* \param offset offset for left shift
*/
void InterfaceFacet::ShiftVertex(const int& offset) {
	// Shift vertex
	/*
	size_t *tmp = (size_t *)malloc(sh.nbIndex * sizeof(size_t));
	for (size_t i = 0; i < sh.nbIndex; i++)
		tmp[i] = GetIndex((int)i + offset);
	free(indices);
	indices = tmp;
	*/

	std::rotate(indices.begin(), indices.begin() + offset, indices.end());
}

/**
* \brief Detect non visible edge (for polygon which contains holes)
*/
void InterfaceFacet::InitVisibleEdge() {

	// Detect non visible edge (for polygon which contains holes)
	std::fill(visible.begin(), visible.end(), true);

	for (int i = 0;i < sh.nbIndex;i++) {

		size_t p11 = GetIndex(i);
		size_t p12 = GetIndex(i + 1);

		for (size_t j = i + 1;j < sh.nbIndex;j++) {

			size_t p21 = GetIndex((int)j);
			size_t p22 = GetIndex((int)j + 1);

			if ((p11 == p22 && p12 == p21) || (p11 == p21 && p12 == p22)) {
				// Invisible edge found
				visible[i] = false;
				visible[j] = false;
			}
		}
	}
}

/**
* \brief Get vertex index from buffer for an idx
* \param idx index
* \return vertex index
*/
size_t InterfaceFacet::GetIndex(int idx) {
	if (idx < 0) {
		return indices[(sh.nbIndex + idx) % sh.nbIndex];
	}
	else {
		return indices[idx % sh.nbIndex];
	}
}

/**
* \brief Get vertex index from buffer for an idx
* \param idx index
* \return vertex index
*/
size_t InterfaceFacet::GetIndex(size_t idx) {
		return indices[idx % sh.nbIndex];
}

/**
* \brief Calculate mesh area and consider the usage of 2 sided meshes
* \param index cell index
* \param correct2sides if correction for 2 sided meshes should be applied (use factor 2)
* \return mesh area
*/
double InterfaceFacet::GetMeshArea(size_t index, bool correct2sides) {
	if (cellPropertiesIds.empty()) return -1.0f;
	if (cellPropertiesIds[index] == -1) {
		return ((correct2sides && sh.is2sided) ? 2.0 : 1.0) / (tRatioU*tRatioV);
	}
	else if (cellPropertiesIds[index] == -2) {
		return 0.0;
	}
	else {
		return ((correct2sides && sh.is2sided) ? 2.0 : 1.0) * meshvector[cellPropertiesIds[index]].area;
	}
}

/**
* \brief Get number of mesh points depending on if its describing a full element, outside of polygon or not
* \param index of mesh
* \return number of mesh points
*/
size_t InterfaceFacet::GetMeshNbPoint(size_t index) {
	size_t nbPts;
	if (cellPropertiesIds[index] == -1) nbPts = 4;
	else if (cellPropertiesIds[index] == -2) nbPts = 0;
	else nbPts = meshvector[cellPropertiesIds[index]].nbPoints;
	return nbPts;
}

/**
* \brief Get the uv coordinate of a point in a mesh
* \param index of mesh
* \param pointId id of the point in the mesh
* \return Vector for mesh point
*/
Vector2d InterfaceFacet::GetMeshPoint(size_t index, size_t pointId) {
	Vector2d result;
	if (cellPropertiesIds.empty()) {
		result.u = 0.0;
		result.v = 0.0;
		return result;
	}
	else {
		int id = cellPropertiesIds[index];
		if (id == -2) {
			result.u = 0.0;
			result.v = 0.0;
			return result;
		}
		else if (id != -1) {
			if (pointId < meshvector[id].nbPoints)
				return meshvector[id].points[pointId];
			else {
				result.u = 0.0;
				result.v = 0.0;
				return result;
			}

		}

		else { //full elem
			double iw = 1.0 / (double)sh.texWidth_precise;
			double ih = 1.0 / (double)sh.texHeight_precise;
			double sx = (double)(index%sh.texWidth);
			double sy = (double)(index / sh.texWidth);
			if (pointId == 0) {
				double u0 = sx * iw;
				double v0 = sy * ih;
				result.u = u0;
				result.v = v0;
				return result;
			}
			else if (pointId == 1) {
				double u1 = (sx + 1.0) * iw;
				double v0 = sy * ih;
				result.u = u1;
				result.v = v0;
				return result;
			}
			else if (pointId == 2) {
				double u1 = (sx + 1.0) * iw;
				double v1 = (sy + 1.0) * ih;
				result.u = u1;
				result.v = v1;
				return result;
			}
			else if (pointId == 3) {
				double u0 = sx * iw;
				double v1 = (sy + 1.0) * ih;
				result.u = u0;
				result.v = v1;
				return result;
			}
			else {
				result.u = 0.0;
				result.v = 0.0;
				return result;
			}
		}
	}
}

/**
* \brief Get the uv coordinate of the central point in a mesh
* \param index of mesh
* \return Vector for point in the center of a mesh
*/
Vector2d InterfaceFacet::GetMeshCenter(size_t index) {
	Vector2d result;
	if (cellPropertiesIds.empty()) {
		result.u = 0.0;
		result.v = 0.0;
		return result;
	}
	if (cellPropertiesIds[index] != -1) {
		if (cellPropertiesIds[index] == -2) {
			result.u = 0.0;
			result.v = 0.0;
			return result;
		}
		else {
			result.u = meshvector[cellPropertiesIds[index]].uCenter;
			result.v = meshvector[cellPropertiesIds[index]].vCenter;
			return result;
		}
	}
	else {
		double iw = 1.0 / (double)sh.texWidth_precise;
		double ih = 1.0 / (double)sh.texHeight_precise;
		double sx = (double)(index%sh.texWidth);
		double sy = (double)(index / sh.texWidth);
		double u0 = sx * iw;
		double v0 = sy * ih;
		double u1 = (sx + 1.0) * iw;
		double v1 = (sy + 1.0) * ih;
		result.u = (float)(u0 + u1) / 2.0f;
		result.v = (float)(v0 + v1) / 2.0f;
		return result;
	}
}

/**
* \brief Get calculated area of a facet (depends on one or double sided)
* \return facet area
*/
double InterfaceFacet::GetArea() {
	return sh.area*(sh.is2sided ? 2.0 : 1.0);
}

/**
* \brief Check if facet is a link facet
* \return true if link facet
*/
bool InterfaceFacet::IsTXTLinkFacet() {
	return ((sh.opacity == 0.0) && (sh.sticking >= 1.0));
}

/**
* \brief Real center coordinate in global space
* \return 3d vector for real center coordinate
*/
Vector3d InterfaceFacet::GetRealCenter() {
	return Project(sh.center, sh.O, sh.N);
}

/**
* \brief Update if profile and texture flag
*/
void InterfaceFacet::UpdateFlags() {

	sh.isProfile = (sh.profileType != PROFILE_NONE);
	//wp.isOpaque = (wp.opacity != 0.0);
	//sh.isTextured = ((texDimW*texDimH) > 0);
    sh.isTextured = ((sh.texWidth_precise*sh.texHeight_precise) > 0);
}

/**
* \brief Detect if 2 facets are in the same plane (orientation preserving) and have same parameters (used by collapse) within a certain threshold
* \param f second facet
* \param threshold threshold for comparison
* \return true if coplanar and equal
*/
bool InterfaceFacet::IsCoplanarAndEqual(InterfaceFacet *f, double threshold) {

	// Detect if 2 facets are in the same plane (orientation preserving)
	// and have same parameters (used by collapse)
	bool equal =
		(fabs(a - f->a) < threshold) &&
		(fabs(b - f->b) < threshold) &&
		(fabs(c - f->c) < threshold) &&
		(fabs(d - f->d) < threshold) &&
		IsEqual(sh.sticking, f->sh.sticking) &&
		IsEqual(sh.opacity, f->sh.opacity) &&
		(sh.is2sided == f->sh.is2sided);

#if defined(MOLFLOW)
	equal = equal &&
		(sh.desorbType == f->sh.desorbType) &&
		IsEqual(sh.reflection.diffusePart, f->sh.reflection.diffusePart) &&
		IsEqual(sh.reflection.specularPart, f->sh.reflection.specularPart) &&
		IsEqual(sh.reflection.cosineExponent, f->sh.reflection.cosineExponent) &&
		(sh.temperature == f->sh.temperature);
	if (sh.area > 0.0 && f->sh.area > 0.0) { //Compare per-area outgassing
		equal = equal &&  IsEqual(sh.outgassing/sh.area,f->sh.outgassing/f->sh.area);
	}
#endif
#if defined(SYNRAD)
	equal = equal &&
		(sh.reflectType == f->sh.reflectType);
#endif
		
	return equal;
}

/**
* \brief Copy properties from another facet
* \param f second facet
* \param copyMesh if mesh values (counters) should also be copied
*/
void InterfaceFacet::CopyFacetProperties(InterfaceFacet *f, bool copyMesh) {
	sh.sticking = f->sh.sticking;
	sh.opacity = f->sh.opacity;

	if (copyMesh) {
		sh.profileType = f->sh.profileType;
	}
	else {
		sh.profileType = PROFILE_NONE;
	}
	sh.is2sided = f->sh.is2sided;
#if defined(MOLFLOW)
	sh.outgassing = f->sh.outgassing;
	sh.desorbType = f->sh.desorbType;
	sh.desorbTypeN = f->sh.desorbTypeN;
	sh.reflection = f->sh.reflection;
	sh.temperature = f->sh.temperature;
#endif
#if defined(SYNRAD)
	sh.reflectType = f->sh.reflectType;
	sh.doScattering = f->sh.doScattering;
	sh.rmsRoughness = f->sh.rmsRoughness;
	sh.autoCorrLength = f->sh.autoCorrLength;
#endif

	sh.superIdx = f->sh.superIdx;
	sh.superDest = f->sh.superDest;
	sh.teleportDest = f->sh.teleportDest;

	if (copyMesh) {
		sh.countAbs = f->sh.countAbs;
		sh.countRefl = f->sh.countRefl;
		sh.countTrans = f->sh.countTrans;
#if defined(MOLFLOW)
		sh.countDes = f->sh.countDes;
		sh.countACD = f->sh.countACD;
#endif
		sh.countDirection = f->sh.countDirection;
		hasMesh = f->hasMesh;
        tRatioU = f->tRatioU;
        tRatioV = f->tRatioV;
    }
	this->UpdateFlags();
	textureVisible = f->textureVisible;
	volumeVisible = f->volumeVisible;
	selected = f->selected;
	
	
	//These are required for the collapse routine
	
	a = f->a;
	b = f->b;
	c = f->c;
	d = f->d;
	
	//wp.area = f->wp.area;
	//planarityError = f->planarityError;
	sh.N = f->sh.N;
	
}

/**
* \brief Divide a facet to new facets, each inheriting their parent's parameters. Create one new facet of each mesh cell.
* \return Facetgroup 
*/
FacetGroup InterfaceFacet::Explode() {
	FacetGroup result;
	result.nbV = 0;
#ifdef MOLFLOW
result.originalPerAreaOutgassing = (sh.area > 0.0) ? sh.outgassing / sh.area : 0.0;
#endif // MOLFLOW

	size_t nonZeroElems = 0, nb = 0;
	for (size_t i = 0; i < sh.texHeight*sh.texWidth; i++) {
		if (cellPropertiesIds[i] != -2) {
			try {
				size_t nbPoints = GetMeshNbPoint(i);
				result.nbV += nbPoints;
				InterfaceFacet *f = new InterfaceFacet(nbPoints);
				f->CopyFacetProperties(this); //Copies absolute outgassing
				result.facets.push_back(f);
			}
			catch (...) {
				for (size_t d = 0; d < i; d++)
					SAFE_DELETE(result.facets[d]);
				throw Error("Cannot reserve memory for new facet(s)");
			}
		}
	}
	return result;
}
