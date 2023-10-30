#include "sng.h"
#include <random>

namespace NWUClustering
{
	// Set SNG algorithm parameters
    void ClusteringAlgo::set_sng_params(double eps, int minPts, int seeds) {
		m_epsSquare =  eps * eps;
		m_minPts =  minPts;
		m_seeds = seeds;
	}

	// Destructor to clean up resources
	ClusteringAlgo::~ClusteringAlgo() {
		m_noise.clear();
		m_visited.clear();
		m_parents.clear();
		m_corepoint.clear();
		m_member.clear();
	}
	
	// Write clusters to the output stream
	void ClusteringAlgo::writeClusters(ostream& o) {
		// Writing point id and cluster id pairs per line; noise has cluster id 0	
		int iMaxID = m_clusters.size(), id, i, j;
		for(i = 0; i < m_pts->m_i_num_points; i++) {

			id = m_pid_to_cid[i];
			o << i << " " << id << endl;
		}

		int sum_points = 0;
		int noise = 0;
		for(i = 0; i < m_clusters.size(); i++) {
			sum_points += m_clusters[i].size();
		}
	
		for (i = 0; i < m_pts->m_i_num_points; i++) {
			if(m_noise[i])
				noise++;
		}	
		
		// Output summary information
		cout << "Total points " << noise + sum_points << " pt_in_cls " << sum_points << " noise " << noise << endl;
		cout << "Number of clusters: " << m_clusters.size() << endl;
	}

	void ClusteringAlgo::writeClusters_uf(ostream& o)
	{
		// Writing point id and cluster id pairs per line; noise has cluster id 0	
		vector <int> clusters;
		clusters.resize(m_pts->m_i_num_points, 0);

		int i, j, sum_points = 0, noise = 0, root, rootcount = 0, tmp;

		// Calculate cluster information
		for(i = 0; i < m_pts->m_i_num_points; i++) {
			root = m_parents[i];

			// Get the number of trees
			if(i == m_parents[i])
				rootcount++;

			// Get the root of the tree containing i	
			while(root != m_parents[root])
				root = m_parents[root];

			// Compress the tree to reduce the height of the branch of the tree to 1
       		j = i;
			while(m_parents[j] != root) {
				tmp  = m_parents[j];
				m_parents[j] = root;
				j = tmp;
			}

			// Count the number of vertices in this tree
			clusters[root]++;
		}

		int count = 0;
		// Process cluster information
		for (i = 0; i < m_pts->m_i_num_points; i++) {
			if (clusters[i] == 1) {
				// Vertex i is noise
				clusters[i] = 0;
				noise++;
			} else if (clusters[i] > 1) {
				// Get the size of cluster count; this will happen only if i is a root
				count++;
				sum_points += clusters[i];
				clusters[i] = count;
			}
			// Skip if i is not a root
   		}

		// Write point id and cluster ids to the output stream
		for (i = 0; i < m_pts->m_i_num_points; i++) {
			o << i << " " << clusters[m_parents[i]] << endl;
		}

		// Output summary information
		cout << "Total points " << noise + sum_points << " pt_in_cls " << sum_points << " noise " << noise << endl;
		cout << "Number of clusters: " << count << endl;

		clusters.clear();
	}

	// Run the Union-Find version of the SNG clustering algorithm
	void run_sng_algo_uf(ClusteringAlgo& sng)
	{			
		int tid, i, pid, j, k, npid, root, root1, root2;
 
        // Initialize clustering parameters
		sng.m_clusters.clear();
		kdtree2_result_vector ne;
			
		// assign parent to itestf
		sng.m_parents.resize(sng.m_pts->m_i_num_points, -1);
		sng.m_member.resize(sng.m_pts->m_i_num_points, 0);
		sng.m_corepoint.resize(sng.m_pts->m_i_num_points, 0);

		int sch, maxthreads = omp_get_max_threads();
		
		// Calculate the thread distribution
		if(sng.m_pts->m_i_num_points % maxthreads == 0)
			sch = sng.m_pts->m_i_num_points/maxthreads;
		else
			sch = sng.m_pts->m_i_num_points/maxthreads + 1;
		
		vector < vector <int > > merge;
		vector <int> init;
		merge.resize(maxthreads, init);

		// Reserve space for merge vector
		for(tid = 0; tid < maxthreads; tid++)
			merge[tid].reserve(sng.m_pts->m_i_num_points);
		
		vector < int > prID;
		prID.resize(sng.m_pts->m_i_num_points, -1);

		vector<int>* ind = sng.m_kdtree->getIndex();		

		double start = omp_get_wtime();	

		#pragma omp parallel private(root, root1, root2, tid, ne, npid, i, j, pid) shared(sch, ind) //, prID)
		{
			int lower, upper;
			tid = omp_get_thread_num();

        	lower = sch * tid;
	        upper = sch * (tid + 1);

        	if(upper > sng.m_pts->m_i_num_points)
		        upper = sng.m_pts->m_i_num_points;

	    	for(i = lower; i < upper; i++) {
				pid = (*ind)[i]; 
				sng.m_parents[pid] = pid;
				prID[pid] = tid;
			}

			#pragma omp barrier

			// Compute core points and neighbors
			for(i = lower; i < upper; i++) {
				pid = (*ind)[i];

				ne.clear();
            	sng.m_kdtree->r_nearest_around_point(pid, 0, sng.m_epsSquare, ne);

				if(ne.size() >= sng.m_minPts) {
					sng.m_corepoint[pid] = 1;
					sng.m_member[pid] = 1;
					
					// Get the root containing pid
					root = pid;

					for (j = 0; j < ne.size(); j++) {
						npid= ne[j].idx;
						if(prID[npid] != tid) {
							merge[tid].push_back(pid);
							merge[tid].push_back(npid);
							continue;
						}

						//Get the root containing npid
						root1 = npid;
						root2 = root;

						if(sng.m_corepoint[npid] == 1 || sng.m_member[npid] == 0) {
							sng.m_member[npid] = 1;
	
							// REMS algorithm to merge the trees
							while(sng.m_parents[root1] != sng.m_parents[root2]) {
								if(sng.m_parents[root1] < sng.m_parents[root2]) {
									if(sng.m_parents[root1] == root1) {
										sng.m_parents[root1] = sng.m_parents[root2];
										root = sng.m_parents[root2];
										break;
									}

				        	        // Splicing
                	        		int z = sng.m_parents[root1];
				            	    sng.m_parents[root1] = sng.m_parents[root2];
                    	       		root1 = z;

								} else {

									if(sng.m_parents[root2] == root2) {
										sng.m_parents[root2] = sng.m_parents[root1];
										root = sng.m_parents[root1];
										break;
									}

					   	       		// Splicing
				        	        int z = sng.m_parents[root2];
            	                	sng.m_parents[root2] = sng.m_parents[root1];					
									root2 = z;
								}
							}
						}
					}
				}
			}
		}

    	// Continue with merging clusters using locks
		int v1, v2, size;
		double stop = omp_get_wtime() ; 
		cout << "Local computation took " << stop - start << " seconds." << endl;

		// Allocate and initiate locks
    	omp_lock_t *nlocks;
		nlocks = (omp_lock_t *) malloc(sng.m_pts->m_i_num_points*sizeof(omp_lock_t));

		//Start = Stop;
		start = omp_get_wtime();

		#pragma omp parallel for private(i) shared(nlocks)
    	for(i = 0; i < sng.m_pts->m_i_num_points; i++) 
      		omp_init_lock(&nlocks[i]); // initialize locks

		#pragma omp parallel for shared(maxthreads, merge, nlocks) private(i, v1, v2, root1, root2, size, tid)
		for(tid = 0; tid < maxthreads; tid++) {
			size = merge[tid].size()/2;

			for(i = 0; i < size; i++) {
                v1 = merge[tid][2 * i];
				v2 = merge[tid][2 * i + 1];
		
				int con = 0;
				if(sng.m_corepoint[v2] == 1)
					con = 1;
				else if(sng.m_member[v2] == 0) {
                	omp_set_lock(&nlocks[v2]);
                    if(sng.m_member[v2] == 0) { // If v2 is not a member yet
                        con = 1;
						sng.m_member[v2] = 1;
                    }
                    	omp_unset_lock(&nlocks[v2]);
				}

				if(con == 1) {
				
					// lLock based approach for merging
					root1 = v1;
					root2 = v2;

					// REMS algorithm with splicing compression techniques
					while (sng.m_parents[root1] != sng.m_parents[root2]) {
						if (sng.m_parents[root1] < sng.m_parents[root2]) {
							
							if(sng.m_parents[root1] == root1) { // root1 is a root
								omp_set_lock(&nlocks[root1]);
								int p_set = false;
								if(sng.m_parents[root1] == root1) { // If root1 is still a root
									sng.m_parents[root1] = sng.m_parents[root2];
									p_set = true;
								}
								omp_unset_lock(&nlocks[root1]);
								if (p_set) // Merge successful
    	      						break;
							}
	
							// splicing
							int z = sng.m_parents[root1];
							sng.m_parents[root1] = sng.m_parents[root2];
							root1 = z;

						} else {
							if(sng.m_parents[root2] == root2) { // root2 is a root			
                	            omp_set_lock(&nlocks[root2]);
                    	        int p_set = false;
                        	    if(sng.m_parents[root2] == root2) { // Check if root2 is a root			
       				                sng.m_parents[root2] = sng.m_parents[root1];
                                	p_set = true;
	                            }
    	                        omp_unset_lock(&nlocks[root2]);
        	                    if (p_set) // Merge Successful
            	                    break;
                	        }
							
							//Splicing
				        	int z = sng.m_parents[root2];
                           	sng.m_parents[root2] = sng.m_parents[root1];
 	                        root2 = z;
						}	
					}
				}
			}
		}

		stop = omp_get_wtime();
		free(nlocks);
		cout << "Merging took " << stop - start << " seconds."<< endl;

		for(tid = 0; tid < maxthreads; tid++)
			merge[tid].clear();
		
		merge.clear();
		ne.clear();
	}
	





	// Run the Sequential DBSCAN clustering algorithm (CLEAN)
	void run_sngcan_algo(ClusteringAlgo& sng) {
		int i, pid, j, k, npid;
		int cid = 1; // cluster id
		vector <int> c;
		c.reserve(sng.m_pts->m_i_num_points);

    	// Initialize clustering parameters
		sng.m_noise.resize(sng.m_pts->m_i_num_points, false);
        sng.m_visited.resize(sng.m_pts->m_i_num_points, false);		
		sng.m_pid_to_cid.resize(sng.m_pts->m_i_num_points, 0);
		sng.m_clusters.clear();

		cout << "SNG SEQUENTIAL ALGORITHM" << endl;

		kdtree2_result_vector ne;
		kdtree2_result_vector ne2;
		//kdtree2_result_vector ne3;
		ne.reserve(sng.m_pts->m_i_num_points);
		ne2.reserve(sng.m_pts->m_i_num_points);

		vector<int>* ind = sng.m_kdtree->getIndex();

		double start = omp_get_wtime() ;		

		// Iterate through points
		for(i = 0; i < sng.m_pts->m_i_num_points; i++) {
			pid = (*ind)[i];

			if (!sng.m_visited[pid]) {
				sng.m_visited[pid] = true;
				ne.clear();
				sng.m_kdtree->r_nearest_around_point(pid, 0, sng.m_epsSquare, ne);
				
				if(ne.size() < sng.m_minPts)
					sng.m_noise[pid] = true;
				else {
					// Start a new cluster
					c.clear();
					c.push_back(pid);
					sng.m_pid_to_cid[pid] = cid;

					// Traverse the neighbors
					for (j = 0; j < ne.size(); j++) {
						npid= ne[j].idx;

						// Not already visited
						if(!sng.m_visited[npid]) {
							sng.m_visited[npid] = true;
	
							// Explore neighbors of neighbors
							ne2.clear();
							sng.m_kdtree->r_nearest_around_point(npid, 0, sng.m_epsSquare, ne2);

							// Sufficient support
							if (ne2.size() >= sng.m_minPts)	{
								// Join Clusters
								for(k = 0; k < ne2.size(); k++)
									ne.push_back(ne2[k]);
							}
						}

						// Not already assigned to a cluster
						if (!sng.m_pid_to_cid[npid]) {
							c.push_back(npid);
							sng.m_pid_to_cid[npid]=cid;
							sng.m_noise[npid] = false;
						}
					}

					sng.m_clusters.push_back(c);
					cid++;
				}	
			}
		}
		
	    double stop = omp_get_wtime();
        cout << "Local computation took " << stop - start << " seconds." << endl;
		cout << "No merging stage in classical SNG"<< endl;
		ind = NULL;
		ne.clear();
		ne2.clear();
	}


	// Run the Sequential SNG clustering algorithm
	void run_sng_algo(ClusteringAlgo& sng) {
		int i, pid, j, k, npid;
		int cid = 1; // cluster id
		vector <int> c;
		c.reserve(sng.m_pts->m_i_num_points);


		// Select m_seeds random points
		vector<int> random_seeds;
		random_seeds.reserve(sng.m_seeds);

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<int> dis(0, sng.m_pts->m_i_num_points - 1);

		cout << "Selected random point(s):" << endl; // Print selected points

		while (random_seeds.size() < sng.m_seeds) {
			int random_index = dis(gen);

			// Check if the point has not been selected before
			if (std::find(random_seeds.begin(), random_seeds.end(), random_index) == random_seeds.end()) {
				random_seeds.push_back(random_index);
				cout << random_index << " "; // Print the selected point
			}
		}

		cout << endl; // Print a newline to separate the list




    	// Initialize clustering parameters
		sng.m_noise.resize(sng.m_pts->m_i_num_points, false);
        sng.m_visited.resize(sng.m_pts->m_i_num_points, false);		
		sng.m_pid_to_cid.resize(sng.m_pts->m_i_num_points, 0);
		sng.m_clusters.clear();

		cout << "SNG SEQUENTIAL ALGORITHM" << endl;
		cout << sng.m_seeds << endl;

		kdtree2_result_vector ne;
		kdtree2_result_vector ne2;
	
		ne.reserve(sng.m_pts->m_i_num_points);
		ne2.reserve(sng.m_pts->m_i_num_points);

		vector<int>* ind = sng.m_kdtree->getIndex();

		double start = omp_get_wtime() ;		

		// Iterate through points
		for (int i = 0; i < sng.m_seeds; i++) {
        	int pid = random_seeds[i];

			if (!sng.m_visited[pid]) {
				sng.m_visited[pid] = true;
				ne.clear();
				sng.m_kdtree->r_nearest_around_point(pid, 0, sng.m_epsSquare, ne);
				
				if(ne.size() < sng.m_minPts)
					sng.m_noise[pid] = true;
				else {
					// Start a new cluster
					c.clear();
					c.push_back(pid);
					sng.m_pid_to_cid[pid] = cid;

					// Traverse the neighbors
					for (j = 0; j < ne.size(); j++) {
						npid= ne[j].idx;

						// Not already visited
						if(!sng.m_visited[npid]) {
							sng.m_visited[npid] = true;
	
							// Explore neighbors of neighbors
							ne2.clear();
							sng.m_kdtree->r_nearest_around_point(npid, 0, sng.m_epsSquare, ne2);

							// Sufficient support
							if (ne2.size() >= sng.m_minPts)	{
								// Join Clusters
								for(k = 0; k < ne2.size(); k++)
									ne.push_back(ne2[k]);
							}
						}

						// Not already assigned to a cluster
						if (!sng.m_pid_to_cid[npid]) {
							c.push_back(npid);
							sng.m_pid_to_cid[npid]=cid;
							sng.m_noise[npid] = false;
						}
					}

					sng.m_clusters.push_back(c);
					cid++;
				}	
			}
		}
		
	    double stop = omp_get_wtime();
        cout << "Local computation took " << stop - start << " seconds." << endl;
		cout << "No merging stage in classical SNG"<< endl;
		ind = NULL;
		ne.clear();
		ne2.clear();
	}
};