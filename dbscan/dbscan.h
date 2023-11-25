
#ifndef _DBSCAN_
#define _DBSCAN_

#include "utils.h"
#include "clusters.h"
#include <queue>
#include <random>

namespace NWUClustering
{
	// Class representing the SNG (Sow & Grow) Focused Clustering Algorithm
	class ClusteringAlgoDBS : public Clusters {
	
	public:
		// Constructor
		ClusteringAlgoDBS(){ }

		// Destructor
		virtual ~ClusteringAlgoDBS();

		//Set parameters for the S&G Algorithm
		void set_dbs_params(double eps, int minPts);
		
		// Write Clusters using the regular SNG Algorithm
		void writeClustersDBS(ostream& o); 

		// Write Clusters using the union-find SNG Algorithm
		void writeClusters_ufDBS(ostream& o); 

	public:
		
		// Parameters for running the SNG algorithm
		double 	m_epsSquare;
		int 	m_minPts;
  
		// Noise vector to mark noise points
        vector<bool> m_noise;
	       	
		// ??? (noise vector) or (Visited vector to keep track of visited points)
        vector<bool> m_visited;

        // Vector to store parent points
        vector<int> m_parents;

        // Vector to mark core points
        vector<int> m_corepoint;

        // Vector to store cluster members
        vector<int> m_member;
	};	

    // Run the union-find SNG algorithm
    void run_dbs_algo_uf(ClusteringAlgoDBS& dbs);

    // Run the regular SNG algorithm
    void run_dbs_algo(ClusteringAlgoDBS& dbs);
};

#endif