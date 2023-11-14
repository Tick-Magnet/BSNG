
#ifndef _SNG_
#define _SNG_

#include "utils.h"
#include "clusters.h"
#include <queue>
#include <random>

namespace NWUClustering
{
	// Class representing the SNG (Sow & Grow) Focused Clustering Algorithm
	class ClusteringAlgo : public Clusters {
	
	public:
		// Constructor
		ClusteringAlgo(){ }

		// Destructor
		virtual ~ClusteringAlgo();

		//Set parameters for the S&G Algorithm
		void set_sng_params(double eps, int minPts, int seeds);
		
		// Write Clusters using the regular SNG Algorithm
		void writeClusters(ostream& o); 

		// Write Clusters using the union-find SNG Algorithm
		void writeClusters_uf(ostream& o); 

	public:
		
		// Parameters for running the SNG algorithm
		double 	m_epsSquare;
		int 	m_minPts;
		int 	m_seeds;
  
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
    void run_sng_algo_uf(ClusteringAlgo& dbs);

    // Run the regular SNG algorithm
    void run_sng_algo(ClusteringAlgo& dbs);
};

#endif