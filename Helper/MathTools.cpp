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
#include "MathTools.h"
#include "Random.h"

#include <cmath>
#include <cstdio>
#include <algorithm> //std::Lower_bound
#include <iterator>
#include <chrono>
#include <cstring> //strdup

bool
IsEqual(const double &a, const double &b, double toleranceRatio) {
	return std::abs(a - b) < Max(1E-99, Max(std::abs(a),std::abs(b))*toleranceRatio);
}

size_t IDX(const int& i, const size_t& nb) {
	//Return circular index restrained within [0..nb[, allows negative index (Python logics: -1=last)
    int ret = i%(int)nb;
    return (ret>=0)?(ret):(ret+nb);
}

size_t IDX(const size_t& i, const size_t& nb) {
	//Return circular index restrained within [0..nb[
	return i%nb;
}


size_t Next(const int& i, const size_t& nb, const bool& inverseDir) {
	//Returns the next element of a circular index (next of last is first)
	//inverseDir is a helper: when true, returns the previous
	return Next((size_t)i,nb,inverseDir);
}

size_t Next(const size_t& i, const size_t& nb, const bool& inverseDir) {
	//Returns the next element of a circular index (next of last is first)
	//inverseDir is a helper: when true, returns the previous
	if (!inverseDir) {
		size_t next=i+1;
		if (next==nb) next = 0;
		return next;
	} else return Previous(i,nb,false);
}

size_t Previous(const size_t& i, const size_t& nb, const bool& inverseDir) {
	//Returns the previous element of a circular index (previous of first is last)
	//inverseDir is a helper: when true, returns the next
	if (!inverseDir) {
		if (i==0) return nb-1;
		else return i-1;
	} else return Next(i,nb,false);
}

size_t Previous(const int& i, const size_t& nb, const bool& inverseDir) {
	return Previous((size_t)i,nb,inverseDir);
}



size_t GetPower2(size_t n) {
// Return a power of 2 which is greater or equal than n
  if((n & (n-1))==0) {
    // already a power of 2
    return n;
  } else {
    // Get the power of 2 above
    int p = 0;
    while(n!=0) { n = n >> 1; p++; }
    return (size_t)1 << p;
  }

}

double RoundAngle(double a) {
// Return a in [-PI,PI]
  double r=a;
  while(r<-PI) r+=2.0*PI;
  while(r> PI) r-=2.0*PI;
  return r;

}

// Converts an angle from degree to rad
double  DegToRad(double a) {
    return a / 180.0 * PI;
}

// Converts an angle from degree to rad
double RadToDeg(double a) {
    return a * 180.0 / PI;
}

char* FormatMemory(size_t size) {
	return FormatMemoryLL((long long)size);
}

char* FormatMemoryLL(long long size) {

	static char ret[256];
	const char *suffixStr[] = { "KB", "MB", "GB", "TB", "PB" };
	auto dSize = static_cast<double>(size);
	int suffix = 0;

	while (dSize >= 1024.0 && suffix < 4) {
		dSize /= 1024.0;
		suffix++;
	}

	if (suffix == 0) {
		sprintf(ret, "%u bytes", static_cast<unsigned int>(size));
	}
	else {
		if (fabs(dSize - floor(dSize)) < 1e-3)
			sprintf(ret, "%.0f%s", dSize, suffixStr[suffix - 1]);
		else
			sprintf(ret, "%.2f%s", dSize, suffixStr[suffix - 1]);
	}
	return ret;

}

[[maybe_unused]] double my_erf(double x)
{
	// constants
	double a1 = 0.254829592;
	double a2 = -0.284496736;
	double a3 = 1.421413741;
	double a4 = -1.453152027;
	double a5 = 1.061405429;
	double p = 0.3275911;

	// Save the sign of x
	int sign = 1;
	if (x < 0)
		sign = -1;
	x = fabs(x);

	// A&S formula 7.1.26
	double t = 1.0 / (1.0 + p*x);
	double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

	return sign*y;
}

double Weigh(const double & a, const double & b, const double & weigh)
{
	return a + (b - a)*weigh;
}

double InterpolateY(const double& x, const std::vector<std::pair<double, double>>& table, const bool& logX, const bool& logY, const bool& allowExtrapolate) {
	return InterpolateXY(x, table, true, logX, logY, allowExtrapolate);
}

double InterpolateX(const double& y, const std::vector<std::pair<double, double>>& table, const bool& logX, const bool& logY, const bool& allowExtrapolate) {
	return InterpolateXY(y, table, false, logX, logY, allowExtrapolate);
}

double InterpolateXY(const double & lookupValue, const std::vector<std::pair<double, double>>& table, const bool & first, const bool& logX, const bool& logY, const bool & allowExtrapolate) {
	//InterpolateX and InterpolateY, only differing by first param
	//Avoids repeated code with minor changes only
	//returns double
	//First: if true, interpolate Y for supplied X, otherwise vica versa

	//firstToSecond: either unused or selector between first and second element of pairs to search for lookupValue
	//param2: either size of pointer A or index of element to search lookupValue in the second element of the pairs

	if (table.size() == 1) return GetElement(table[0], !first);

	int lowerIndex = my_lower_bound(lookupValue, table, first);

	if (lowerIndex == -1) {
		lowerIndex = 0;
		if (!allowExtrapolate) return GetElement(table[lowerIndex], !first); //return first element
	}
	else if (lowerIndex == (static_cast<int>(table.size()) - 1)) {
		if (allowExtrapolate) {
			lowerIndex = static_cast<int>(table.size()) - 2;
		}
		else return GetElement(table[lowerIndex], !first); //return last element
	}

	double delta = (first ? logX : logY) ? log10(GetElement(table[lowerIndex + 1], first)) - log10(GetElement(table[lowerIndex], first)) : GetElement(table[lowerIndex + 1], first) - GetElement(table[lowerIndex], first);
	double overshoot = (first ? logX : logY) ? log10(lookupValue) - log10(GetElement(table[lowerIndex], first)) : lookupValue - GetElement(table[lowerIndex], first);

	if (first ? logY : logX) return Pow10(Weigh(log10(GetElement(table[lowerIndex], !first)),log10(GetElement(table[lowerIndex + 1], !first)),overshoot / delta));
	else return Weigh(GetElement(table[lowerIndex], !first),GetElement(table[lowerIndex + 1], !first),overshoot / delta);

}

double InterpolateVectorX(const double& y, const std::vector < std::pair<double, std::vector<double>>> & table, const size_t& elementIndex, const bool& logX, const bool& logY, const bool& allowExtrapolate) {
	//Interpolate X belonging to supplied Y
	//Avoids repeated code with minor changes only
	//returns double

	//firstToSecond: either unused or selector between first and second element of pairs to search for lookupValue
	//param2: either size of pointer A or index of element to search lookupValue in the second element of the pairs
	//elementIndex: which index of the vector to use for interpolation
	const bool first = false;
	double lookupValue = y;

	if (table.size() == 1) return GetElement(table[0], !first, elementIndex);

	int lowerIndex = my_lower_bound(lookupValue, table, first, elementIndex);

	if (lowerIndex == -1) {
		lowerIndex = 0;
		if (!allowExtrapolate) return GetElement(table[lowerIndex], !first, elementIndex); //return first element
	}
	else if (lowerIndex == (static_cast<int>(table.size()) - 1)) {
		if (allowExtrapolate) {
			lowerIndex = static_cast<int>(table.size()) - 2;
		}
		else return GetElement(table[lowerIndex], !first, elementIndex); //return last element
	}

	double delta = (first ? logX : logY) ? log10(GetElement(table[lowerIndex + 1], first, elementIndex)) - log10(GetElement(table[lowerIndex], first, elementIndex)) : GetElement(table[lowerIndex + 1], first, elementIndex) - GetElement(table[lowerIndex], first, elementIndex);
	double overshoot = (first ? logX : logY) ? log10(lookupValue) - log10(GetElement(table[lowerIndex], first, elementIndex)) : lookupValue - GetElement(table[lowerIndex], first, elementIndex);

	if (first? logY : logX) return Pow10(Weigh(log10(GetElement(table[lowerIndex], !first, elementIndex)),log10(GetElement(table[lowerIndex + 1], !first, elementIndex)),overshoot / delta));
	else return Weigh(GetElement(table[lowerIndex], !first, elementIndex),GetElement(table[lowerIndex + 1], !first, elementIndex),overshoot / delta);
}

std::vector<double> InterpolateVectorY(const double& x, const std::vector<std::pair<double, std::vector<double>>>& table, const bool& logX, const bool& logY, const bool& allowExtrapolate) {
	//Same as InterpolateY but returns a vector.
	//Must repeat most of code because C++ doesn't allow runtime evaluated return-type (and only 'bool first' decides what to return)
	if (table.size() == 1) return table[0].second;

	int lowerIndex = my_lower_bound(x, table, true, 0);

	if (lowerIndex == -1) {
		lowerIndex = 0;
		if (!allowExtrapolate) return table[lowerIndex].second; //return first element
	}
	else if (lowerIndex == (static_cast<int>(table.size()) - 1)) {
		if (allowExtrapolate) {
			lowerIndex = static_cast<int>(table.size()) - 2;
		}
		else return table[lowerIndex].second; //return last element
	}

	double delta = (logX) ? log10(table[lowerIndex + 1].first) - log10(table[lowerIndex].first) : table[lowerIndex + 1].first - table[lowerIndex].first;
	double overshoot = (logX) ? log10(x) - log10(table[lowerIndex].first) : x - table[lowerIndex].first;

	size_t distrYsize = table[0].second.size();
	std::vector<double> result; result.resize(distrYsize);
	for (size_t e = 0; e < distrYsize; e++)
	{
		if (logY) result[e]=Pow10(Weigh(log10(table[lowerIndex].second[e]),log10(table[lowerIndex + 1].second[e]), overshoot / delta)); //log-log interpolation
		else result[e]=Weigh(table[lowerIndex].second[e],table[lowerIndex + 1].second[e],overshoot / delta);
	}
	return result;
}




/*
double FastLookupY(const double& x, const std::vector<std::pair<double, double>>& table, const bool& limitToBounds) {
	//Function inspired by http://stackoverflow.com/questions/11396860/better-way-than-if-else-if-else-for-linear-interpolation
	_ASSERTE(table.size());
	if (table.size() == 1) return table[0].second; //constant value

												   // Assumes that table .first is SORTED AND EQUIDISTANT
												   // Check if x is out of bound
	std::vector<std::pair<double, double> >::const_iterator lower, upper;
	bool outOfLimits = false;

	if (x >= table.back().first) {
		if (limitToBounds) return table.back().second;
		else {
			outOfLimits = true;
			lower = upper = table.end() - 1;
			lower--;
		}
	}
	else if (x < table[0].first) {
		if (limitToBounds) return table[0].second;
		else {
			outOfLimits = true;
			lower = upper = table.begin();
			upper++;
		}
	}

	if (!outOfLimits) {
		double distanceX = table[1].first - table[0].first;
		size_t lowerIndex = (int)((x - table[0].first) / distanceX);
		lower = upper = table.begin() + (lowerIndex + 1);
		// Corner case
		if (upper == table.begin()) return upper->second;
		lower--;
	}
	double result = lower->second + (upper->second - lower->second)*(x - lower->first) / (upper->first - lower->first);
	return result;
}
*/


int my_lower_bound(const double & key, const double* A,const size_t& size)
//"iterative" version of algorithm, modified from https://en.wikipedia.org/wiki/Binary_search_algorithm
//key: searched value
//A: the lookup table
//first: either unused or selector between first and second element of pairs
//returns index of last lower value, or -1 if key not found

// GetElement: chooses first or second member of a pair (avoids writing this function twice)

{
	int l = 0;
	int h = (int)size; // Not n - 1
	while (l < h) {
		int mid = (l + h) / 2;
		if (key <= A[mid]) {
			h = mid;
		}
		else {
			l = mid + 1;
		}
	}
	return l - 1;
}

int my_lower_bound(const double & key, const std::vector<double>& A)
//"iterative" version of algorithm, modified from https://en.wikipedia.org/wiki/Binary_search_algorithm
//key: searched value
//A: the lookup table
//first: either unused or selector between first and second element of pairs
//returns index of last lower value, or -1 if key not found

// GetElement: chooses first or second member of a pair (avoids writing this function twice)

{
	int l = 0;
	int h = (int)A.size(); // Not n - 1
	while (l < h) {
		int mid = (l + h) / 2;
		if (key <= A[mid]) {
			h = mid;
		}
		else {
			l = mid + 1;
		}
	}
	return l - 1;
}

int my_lower_bound(const double & key, const std::vector<std::pair<double, double>>& A, const bool & first)
//"iterative" version of algorithm, modified from https://en.wikipedia.org/wiki/Binary_search_algorithm
//key: searched value
//A: the lookup table
//first: either unused or selector between first and second element of pairs
//returns index of last lower value, or -1 if key not found

// GetElement: chooses first or second member of a pair (avoids writing this function twice)

{
	int l = 0;
	int h = (int)A.size(); // Not n - 1
	while (l < h) {
		int mid = (l + h) / 2;
		if (key <= GetElement(A[mid], first)) {
			h = mid;
		}
		else {
			l = mid + 1;
		}
	}
	return l - 1;
}

int my_lower_bound(const double & key, const std::vector<std::pair<double, std::vector<double>>>& A, const bool & first, const size_t & elementIndex)
//"iterative" version of algorithm, modified from https://en.wikipedia.org/wiki/Binary_search_algorithm
//key: searched value
//A: the lookup table
//first: either unused or selector between first and second element of pairs
//returns index of last lower value, or -1 if key not found

// GetElement: chooses first or second member of a pair (avoids writing this function twice)

{
	int l = 0;
	int h = (int)A.size(); // Not n - 1
	while (l < h) {
		int mid = (l + h) / 2;
		if (key <= GetElement(A[mid],first,elementIndex)) {
			h = mid;
		}
		else {
			l = mid + 1;
		}
	}
	return l-1;
}

/*!
 * @brief Lookup the index of the interval related to a given key
 * @param key specific moment
 * @param moments vector of time intervals
 * @return -1 if moment doesnt relate to an interval, else index of moment (+1 to account for [0]== steady state)
 */
[[maybe_unused]] int LookupMomentIndex(const double & key, const std::vector<std::pair<double, double>>& moments){
    /*int lowerBound = my_lower_bound(key, moments, true);
    if(lowerBound != -1 && lowerBound < moments.size()){
        if(moments[lowerBound].first <= key && key < moments[lowerBound].second){
            return lowerBound + 1;
        }
    }
    return -1;*/
    auto lowerBound = std::lower_bound(moments.begin(), moments.end(), std::make_pair(key,key));
    --lowerBound; //even moments.end() can be a bound

    if(lowerBound->first <= key && key < lowerBound->second){
        return static_cast<int>(std::distance(moments.begin(), lowerBound) + 1);
    }
    return -1;
}

/*!
 * @brief Lookup the index of the interval related to a given key and a start position for accelerated lookup
 * @param key specific moment
 * @param moments vector of time intervals
 * @param startIndex offset to only look in a subset of moments
 * @return -1 if moment doesnt relate to an interval, else index of moment (+1 to account for [0]== steady state)
 */
int LookupMomentIndex(const double & key, const std::vector<std::pair<double, double>>& moments, const size_t &startIndex){

    if(!moments.empty()) {
        auto lowerBound = std::lower_bound(moments.begin() + startIndex, moments.end(), std::make_pair(key, key));
        if(lowerBound == moments.begin())
            return -1;
        --lowerBound; //even moments.end() can be a bound

        if (lowerBound->first <= key && key < lowerBound->second) {
            return static_cast<int>(std::distance(moments.begin() + startIndex, lowerBound) + startIndex + 1);
        }
    }
    return -1;
}


/*double QuadraticInterpolateX(const double& y,
	const double& a, const double& b, const double& c,
	const double& FA, const double& FB, const double& FC, MersenneTwister& randomGenerator) {
	double amb = a - b;
	double amc = a - c;
	double bmc = b - c;
	double amb_amc = amb*amc;
	double amc_bmc = amc*bmc;
	double divisor = (2 * (-(FA / (amb_amc)) + FB / (amb_amc)+FB / (amc_bmc)-FC / (amc_bmc)));
	                 
	if (fabs(divisor) < 1e-30) {
		//Divisor is 0 when the slope is 1st order (a straight line) i.e. (FC-FB)/(c-b) == (FB-FA)/(b-a)
		if ((FB - FA) < 1e-30) {
			//FA==FB, shouldn't happen
			return a + rnd()*(b - a);
		}
		else {
			//Inverse linear interpolation
			return a + (y - FA) * (b - a)/ (FB - FA);
		}
	}
	else {
		//(reverse interpolate y on a 2nd order polynomial fitted on {a,FA},{b,FB},{c,FC} where FA<y<FB):
		//Root of Lagrangian polynomials, solved by Mathematica
		return (FA / (amb)-(a*FA) / (amb_amc)-(b*FA) / (amb_amc)-FB / (amb)+(a*FB) / (amb_amc)+(b*FB) / (amb_amc)+(a*FB) / (amc_bmc)+(b*	FB) / (amc_bmc)-(a*FC) / (amc_bmc)
			-(b*FC) / (amc_bmc)-sqrt(Sqr(-(FA / (amb)) + (a*FA) / (amb_amc)+(b*FA) / (amb_amc)+FB / (amb)-(a*FB) / (amb_amc)-(b*FB) / (amb_amc)-(a*FB) / (amc_bmc)-(b*FB)
				/ (amc_bmc)+(a*FC) / (amc_bmc)+(b*FC) / (amc_bmc)) - 4 * (-(FA / (amb_amc)) + FB / (amb_amc)+FB / (amc_bmc)-FC / (amc_bmc))*(-FA + (a*FA) / (amb)-(a*b*FA)
					/ (amb_amc)-(a*FB) / (amb)+(a*b*FB) / (amb_amc)+(a*b*FB) / (amc_bmc)-(a*b*FC) / (amc_bmc)+y))) / divisor;
	}
}*/

int weighed_lower_bound_X(const double & key, const double & weigh, double * A, double * B, const size_t & size)
{
	//interpolates among two lines of a cumulative distribution
	//all elements of line 1 and line 2 must be monotonously increasing (except equal consecutive values)
	//key: lookup value
	//weigh: between 0 and 1 (0: only first distribution, 1: only second distribution)
	//A* and B* : pointers to arrays of 'size' number of CDF values. The first value (not included) is assumed to be 0, the last (not included) is assumed to be 1
	//return value: lower index. If -1, then key is smaller than first element, if 'size-1', then key is larger than last element
	
		if (size == 0) return -1;
		if (size == 1) {
			double weighed = Weigh(A[0], B[0], weigh);
			if (key < weighed) return -1;
			else return 0;
		}
		int L = 0;
		int R = (int)(size - 1);
		// continue searching while [imin,imax] is not empty
		int M; double weighed,nextWeighed;
		while (L<=R)
		{
			M = (L + R) / 2;
			weighed = Weigh(A[M], B[M], weigh);
			nextWeighed = Weigh(A[M + 1], B[M + 1], weigh);
			if (weighed <= key && key < nextWeighed) {
				// key found at index M
				return M;
			}
			else if (weighed < key) {
				L = M + 1;
			}
			else  {
				R = M - 1;
			}
		}
		//Not found
		if (M == 0)
			return -1; //key lower than first element
		else
			return (int)size - 1; //key larger than last element
	
}

inline double GetElement(const std::pair<double, double>& pair, const bool & first) {
	return first ? pair.first : pair.second;
}

inline double GetElement(const std::pair<double, std::vector<double>>& pair, const bool & first, const size_t & elementIndex) {
	return first ? pair.first : pair.second[elementIndex];
}

// Deprecated timer method, use Chronometer class now
[[maybe_unused]] size_t GetSysTimeMs()
{
	std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(
		std::chrono::system_clock::now().time_since_epoch()
		);
	return ms.count();
}

double Pow10(const double& a) {
	return pow(10,a);
}

std::tuple<double, double> CartesianToPolar(const Vector3d& incidentDir, const Vector3d& normU, const Vector3d& normV, const Vector3d& normN) {

    //input vectors need to be normalized

    // Get polar coordinates of the incoming particule direction in the (U,V,N) facet space.
    // Note: The facet is parallel to (U,V), we use its (nU,nV,N) orthonormal basis here.
    // (nU,nV,N) and (x,y,z) are both left handed

    // Cartesian(x,y,z) to polar in (nU,nV,N) transformation

    // Basis change (x,y,z) -> (nU,nV,N)
    // We use the fact that (nU,nV,N) belongs to SO(3)
    double u = Dot(incidentDir, normU);
    double v = Dot(incidentDir, normV);
    double n = Dot(incidentDir, normN);
    Saturate(n, -1.0, 1.0); //sometimes rounding errors do occur, 'acos' function would return no value for theta

    // (u,v,n) -> (theta,phi)

    double inTheta = acos(n);              // Angle to normal (PI/2 => PI
    //double rho = sqrt(v*v + u*u);
    //double inPhi = asin(v / rho);     //At this point, -PI/2 < inPhi < PI/2
    //if (u < 0.0) inPhi = PI - inPhi;  // Angle to U
    double inPhi = atan2(v, u); //-PI .. PI, and the angle is 0 when pointing towards u
    return { inTheta, inPhi };
}

Vector3d
PolarToCartesian(const Vector3d &nU, const Vector3d &nV, const Vector3d &nN, const double &theta, const double &phi,
                 const bool &reverse) {

    //returns sHandle->currentParticle.direction

    //Vector3d U, V, N;
    //double u, v, n;

    // Polar in (nU,nV,N) to Cartesian(x,y,z) transformation  ( nU = U/|U| , nV = V/|V| )
    // tetha is the angle to the normal of the facet N, phi to U
    // ! See Geometry::InitializeGeometry() for further informations on the (U,V,N) basis !
    // (nU,nV,N) and (x,y,z) are both left handed

    double u = sin(theta)*cos(phi);
    double v = sin(theta)*sin(phi);
    double n = cos(theta);
    //#endif

    // Get the (nU,nV,N) orthonormal basis of the facet
    Vector3d U = nU; // nU
    Vector3d V = nV; // nV
    Vector3d N = nN; // nN
    if (reverse) {
        N = -1.0 * N;
    }
    // Basis change (nU,nV,N) -> (x,y,z)
    return u*U + v*V + n * N;
}

/*
int isinf(double x)
{
#if defined _MSC_VER
	typedef unsigned __int64 uint64;
#else
	typedef uint64_t uint64;
#endif
	union { uint64 u; double f; } ieee754;
	ieee754.f = x;
	return ((unsigned)(ieee754.u >> 32) & 0x7fffffff) == 0x7ff00000 &&
		((unsigned)ieee754.u == 0);
}

int isnan(double x)
{
#if defined _MSC_VER
	typedef unsigned __int64 uint64;
#else
	typedef uint64_t uint64;
#endif
	union { uint64 u; double f; } ieee754;
	ieee754.f = x;
	return ((unsigned)(ieee754.u >> 32) & 0x7fffffff) +
		((unsigned)ieee754.u != 0) > 0x7ff00000;
}*/