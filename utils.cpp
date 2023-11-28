
#include "utils.h"

// Function to find the Kth element in a vector without recursion
// Parameters:
//   A - Input vector of floats
//   K - The index of the desired Kth element
// Returns:
//   The Kth element in the sorted vector

float findKMedian(vector<float>& A, int K) {
    int l = 0;
    int m = A.size() - 1;

    while (l < m) {
        float x = A[K];
        int i = l;
        int j = m;

        // Partition the vector into two segments
        do {
            while (A[i] < x) i++;
            while (x < A[j]) j--;
            
            if (i <= j) {
                swap(A[i], A[j]);
                i++;
                j--;
            }
        } while (i <= j);

        // Determine which segment contains the desired Kth element
        if (j < K) {
            l = i;
        } else if (K < i) {
            m = j;
        } else {
            // Kth element found
            return A[K];
        }
    }

    // Return the Kth element
    return A[K];
}