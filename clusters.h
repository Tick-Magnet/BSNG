#ifndef _CLUSTER_
#define _CLUSTER_

#include "utils.h"
#include "kdtree2.hpp"

namespace NWUClustering {

    // Struct to represent a set of points
    struct Points {
        array2dfloat m_points; // 2D array of points
        int m_i_dims;          // Number of dimensions
        int m_i_num_points;    // Number of points
    };

    // Class for handling clusters of points
    class Clusters {
    public:
        Clusters() : m_pts(nullptr), m_kdtree(nullptr) { }
        virtual ~Clusters();

        // Read data from a file
        int read_file(char* infilename, int isBinaryFile);

        // Build a k-d tree for the points
        int build_kdtree();

    public:
        Points* m_pts;                  // Pointer to points data
        kdtree2* m_kdtree;              // Pointer to k-d tree
        vector<int> m_pid_to_cid;       // Mapping from point id to cluster id
        vector<vector<int>> m_clusters; // Vector of clusters
        int m_parcent_of_data;          // Percentage of data
        char* csvOutputFilename;		// Filename of CSV output file
    };
}

#endif
