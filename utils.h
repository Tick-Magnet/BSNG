
#ifndef _UTILS_
#define _UTILS_

#include <omp.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib> // Include the standard library for C functions
#include <getopt.h>

using namespace std;

// Define a type alias for the coordinate type
typedef float point_coord_type;

// Define a type alias for a 2D array of floating-point numbers
typedef vector<vector<point_coord_type>> array2dfloat;

// Function to find the Kth median in a vector of floats
// Parameters:
//   A - Input vector of floats
//   K - The desired Kth median
// Returns:
//   The Kth median value

float findKMedian(vector<float>& A, int K);

#endif
