#include "MathTools.h"
#include "GLTypes.h" //BOOL
#include <math.h>
#include <cstdio>
#include <algorithm> //std::Lower_bound

int IsEqual(const double &a, const double &b, double tolerance) {
	return fabs(a - b) < tolerance;
}

int GetPower2(int n) {
// Return a power of 2 which is greater or equal than n
  if((n & (n-1))==0) {
    // already a power of 2
    return n;
  } else {
    // Get the power of 2 above
    int p = 0;
    while(n!=0) { n = n >> 1; p++; }
    return 1 << p;
  }

}

double RoundAngle(double a) {
// Return a in [-PI,PI]
  double r=a;
  while(r<-PI) r+=2.0*PI;
  while(r> PI) r-=2.0*PI;
  return r;

}

char* FormatMemory(size_t size) {
	return FormatMemoryLL((long long)size);
}

char* FormatMemoryLL(long long size) {

	static char ret[256];
	const char *suffixStr[] = { "KB", "MB", "GB", "TB", "PB" };
	double dSize = (double)size;
	int suffix = 0;

	while (dSize >= 1024.0 && suffix < 4) {
		dSize /= 1024.0;
		suffix++;
	}

	if (suffix == 0) {
		sprintf(ret, "%u bytes", (unsigned int)size);
	}
	else {
		if (fabs(dSize - floor(dSize)) < 1e-3)
			sprintf(ret, "%.0f%s", dSize, suffixStr[suffix - 1]);
		else
			sprintf(ret, "%.2f%s", dSize, suffixStr[suffix - 1]);
	}
	return ret;

}

int my_binary_search(const double& key, double* A, const size_t& size)
//"iterative" version of algorithm, modified from https://en.wikipedia.org/wiki/Binary_search_algorithm
//key: searched value
//A: ordered arrray of lookup values
//size: array length
//returns index of last lower value, or -1 if key not found

{
	int imin = 0;
	int imax = size - 1;
	// continue searching while [imin,imax] is not empty
	while (imin <= imax)
	{
		// calculate the midpoint for roughly equal partition
		int imid = (imin + imax) / 2;
		if (imid == size - 1 || imid == 0 || (A[imid] <= key && key < A[imid + 1])) {
			// key found at index imid
			return imid;
		}
		// determine which subarray to search
		else if (A[imid] < key) {
			// change min index to search upper subarray
			imin = imid + 1;
		}
		else
		{
			// change max index to search lower subarray
			imax = imid - 1;
		}
	}
	// key was not found
	return -1;
}

int my_binary_search(const double& key, std::vector<double> A, const size_t& size) {
	return my_binary_search(key, &(A[0]), size);
}

double my_erf(double x)
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

BOOL compare_second(const std::pair<double, double>& lhs, const std::pair<double, double>& rhs) {
	return (lhs.second<rhs.second);
}


double InterpolateY(double x, const std::vector<std::pair<double, double>>& table, BOOL limitToBounds, BOOL logarithmic) {
	//Function inspired by http://stackoverflow.com/questions/11396860/better-way-than-if-else-if-else-for-linear-interpolation
	_ASSERTE(table.size());
	if (table.size() == 1) return table[0].second; //constant value
												   // Assumes that "table" is sorted by .first
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

	// INFINITY is defined in math.h in the glibc implementation
	if (!outOfLimits) {
		lower = upper = std::lower_bound(table.begin(), table.end(), std::make_pair(x, -MY_INFINITY));
		// Corner case
		if (upper == table.begin()) return upper->second;
		lower--;
	}
	if (logarithmic) return exp(log(lower->second) + (log(upper->second) - log(lower->second))
		*(log(x) - log(lower->first)) / (log(upper->first) - log(lower->first)));
	else return lower->second + (upper->second - lower->second)*(x - lower->first) / (upper->first - lower->first);

}

double InterpolateX(double y, const std::vector<std::pair<double, double>>& table, BOOL limitToBounds) {
	//Function inspired by http://stackoverflow.com/questions/11396860/better-way-than-if-else-if-else-for-linear-interpolation
	_ASSERTE(table.size());
	if (table.size() == 1) return table[0].second; //constant value

												   // Assumes that "table" is sorted by .second
												   // Check if y is out of bound
	std::vector<std::pair<double, double> >::const_iterator lower, upper;
	BOOL outOfLimits = FALSE;

	if (y >= table.back().second) {
		if (limitToBounds) return table.back().first;
		else {
			outOfLimits = TRUE;
			lower = upper = table.end() - 1;
			lower--;
		}
	}
	else if (y < table[0].second) {
		if (limitToBounds) return table[0].first;
		else {
			outOfLimits = TRUE;
			lower = upper = table.begin();
			upper++;
		}
	}

	// INFINITY is defined in math.h in the glibc implementation
	if (!outOfLimits) {
		lower = upper = std::lower_bound(table.begin(), table.end(), std::make_pair(MY_INFINITY, y), compare_second);
		// Corner case
		if (upper == table.begin()) return upper->first;
		lower--;
	}
	return lower->first + (upper->first - lower->first)*(y - lower->second) / (upper->second - lower->second);
}

double FastLookupY(double x, const std::vector<std::pair<double, double>>& table, BOOL limitToBounds) {
	//Function inspired by http://stackoverflow.com/questions/11396860/better-way-than-if-else-if-else-for-linear-interpolation
	_ASSERTE(table.size());
	if (table.size() == 1) return table[0].second; //constant value

												   // Assumes that table .first is SORTED AND EQUIDISTANT
												   // Check if x is out of bound
	std::vector<std::pair<double, double> >::const_iterator lower, upper;
	BOOL outOfLimits = FALSE;

	if (x >= table.back().first) {
		if (limitToBounds) return table.back().second;
		else {
			outOfLimits = TRUE;
			lower = upper = table.end() - 1;
			lower--;
		}
	}
	else if (x < table[0].first) {
		if (limitToBounds) return table[0].second;
		else {
			outOfLimits = TRUE;
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