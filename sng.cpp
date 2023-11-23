#include "sng.h"


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





// Run the Union-Find version of the Sow-and-Grow (SNG) clustering algorithm
void run_sng_algo_uf(ClusteringAlgo& sng) {	
    cout << endl; 
    cout << "SNG Parallel ALGORITHM" << endl;
    cout << sng.m_seeds << ": seeds" << endl;
    cout << sng.m_pts->m_i_num_points << ": total points" << endl;

    // Select m_seeds random points
    vector<int> random_seeds;
    random_seeds.reserve(sng.m_seeds);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(0, sng.m_pts->m_i_num_points - 1);

    cout << "Selected random point(s):" << endl; 

    while (random_seeds.size() < sng.m_seeds) {
        int random_index = dis(gen);

        // Check if the point has not been selected before
        if (std::find(random_seeds.begin(), random_seeds.end(), random_index) == random_seeds.end()) {
            random_seeds.push_back(random_index);
            cout << random_index << "!" << endl; 
        }
    }
    cout << endl; 

    // Initialize clustering parameters
    sng.m_clusters.clear();
    kdtree2_result_vector neighbors;

    sng.m_parents.resize(sng.m_pts->m_i_num_points, -1);
    sng.m_member.resize(sng.m_pts->m_i_num_points, 0);
    sng.m_corepoint.resize(sng.m_pts->m_i_num_points, 0);

    int thread_id, i, j, k, neighbor_point_id, point_id, root, root1, root2;
    int max_threads = omp_get_max_threads();
    
    vector<vector<int>> merge;
    vector<int> init;
    merge.resize(max_threads, init);

    vector<int> thread_id_map;
    vector<int>* point_indices = sng.m_kdtree->getIndex();
    thread_id_map.resize(sng.m_pts->m_i_num_points, -1); 
    
    //Setting all points to the same thread.
    double start_time = omp_get_wtime();

    // Reserve space for the merge vector
    for (i = 0; i < max_threads; i++) 
    {
        merge[i].reserve(sng.m_pts->m_i_num_points);
    }


	//Loop through every point assign default parent.
	for (int point_id = 0; point_id < sng.m_pts->m_i_num_points; ++point_id) {
    	sng.m_parents[point_id] = point_id;
	}


    // Initialize queue for point expansion
    vector<queue<int>> pointQueues(max_threads); // Create a vector of queues, one for each thread

    #pragma omp parallel private(root, root1, root2, thread_id, neighbors, neighbor_point_id, i, j, point_id) shared(pointQueues, sng, random_seeds)
    {
        thread_id = omp_get_thread_num();
        queue<int>& pointQueue = pointQueues[thread_id];

        // Distribute seeds to threads without overlap
        int seeds_per_thread = sng.m_seeds / max_threads;
        int start_seed = thread_id * seeds_per_thread;
        int end_seed = (thread_id == max_threads - 1) ? sng.m_seeds : start_seed + seeds_per_thread;

        #pragma omp critical
        {
            cout << "Thread " << thread_id << " processes seeds from " << start_seed << " to " << end_seed - 1 << endl;
        }


        for (int seed_idx = start_seed; seed_idx < end_seed; ++seed_idx) {
            pointQueue.push(random_seeds[seed_idx]);
        }

        while (!pointQueue.empty()) {

            #pragma omp critical
            {
                // Debugging code to print the contents of the queue
                cout << "Contents of " << thread_id << " pointQueue: ";
                queue<int> tempQueue = pointQueue; // Create a temporary queue for printing
                while (!tempQueue.empty()) {
                    cout << tempQueue.front() << " ";
                    tempQueue.pop();
                }
                cout << endl;
            }
			

            //Select Front of Queue
			int currentPoint = pointQueue.front();
            pointQueue.pop();
            //cout << "CurrentPoint: " << currentPoint << endl;
            //cout << "m_member: " << sng.m_member[currentPoint] <<endl;
            //cout << "m_corepoint: " << sng.m_corepoint[currentPoint] <<endl;

            // Find core points and neighbors for currentPoint
            neighbors.clear();
            sng.m_kdtree->r_nearest_around_point(currentPoint, 0, sng.m_epsSquare, neighbors);

            if (neighbors.size() >= sng.m_minPts && sng.m_corepoint[currentPoint] != 1) {
                //cout << "Added: " << currentPoint << " as a Core Point" << endl;
				sng.m_corepoint[currentPoint] = 1;
                //cout << "Updated m_corepoint: " << sng.m_corepoint[currentPoint] <<endl;
				sng.m_member[currentPoint] = 1;

                // Get the root containing currentPoint
                root = currentPoint;
            
                for (j = 0; j < neighbors.size(); j++) {
                    
					neighbor_point_id = neighbors[j].idx;
                    //cout << "Neighborhood Point ID: " << neighbor_point_id <<endl;

                    root1 = neighbor_point_id;
                    root2 = root;

                    if (sng.m_member[neighbor_point_id] == 0) {
                        //cout << "Assigning Membership to: " << neighbor_point_id << endl;
                        
                        sng.m_member[neighbor_point_id] = 1;

                        // Union-Find algorithm to merge the trees
                        while (sng.m_parents[root1] != sng.m_parents[root2]) {
							//cout << "Union Find" << endl;
                            if (sng.m_parents[root1] < sng.m_parents[root2]) {
                                if (sng.m_parents[root1] == root1) {
                                    sng.m_parents[root1] = sng.m_parents[root2];
                                    root = sng.m_parents[root2];
                                    //cout << root << "New Root" << endl;
                                    break;
                                }
                                // Splicing
                                int temp = sng.m_parents[root1];
                                sng.m_parents[root1] = sng.m_parents[root2];
                                root1 = temp;
                            } else {
                                if (sng.m_parents[root2] == root2) {
                                    sng.m_parents[root2] = sng.m_parents[root1];
                                    root = sng.m_parents[root1];
                                    break;
                                }
                                // Splicing
                                int temp = sng.m_parents[root2];
                                sng.m_parents[root2] = sng.m_parents[root1];
                                root2 = temp;
                            }
                        }
                    }	
					
                    //If the neighbor hasn't been searched and it's the same point, add it to the queue.
					//cout << "Checking if " << neighbor_point_id << ": should be added to queue" << endl;
					if (sng.m_corepoint[neighbor_point_id] == 0){
                        sng.m_corepoint[neighbor_point_id] = -1;
						pointQueue.push(neighbor_point_id);
						//cout << "Added Point to Queue: " << neighbor_point_id << endl;
					}
                }
            } 
        }
    }

    // Continue with merging clusters using locks
    int vertex1, vertex2, merge_size;
    double stop_time = omp_get_wtime(); 
    cout << "Local computation took " << stop_time - start_time << " seconds." << endl;

    // Allocate and initiate locks
    omp_lock_t *node_locks;
    node_locks = (omp_lock_t *) malloc(sng.m_pts->m_i_num_points * sizeof(omp_lock_t));

    start_time = omp_get_wtime();

    #pragma omp parallel for private(i) shared(node_locks)
    for (i = 0; i < sng.m_pts->m_i_num_points; i++) 
        omp_init_lock(&node_locks[i]); // Initialize locks

    #pragma omp parallel for shared(max_threads, merge, node_locks) private(i, vertex1, vertex2, root1, root2, merge_size, thread_id)
    for (thread_id = 0; thread_id < max_threads; thread_id++) {
        merge_size = merge[thread_id].size() / 2;

        for (i = 0; i < merge_size; i++) {
            vertex1 = merge[thread_id][2 * i];
            vertex2 = merge[thread_id][2 * i + 1];
            int should_merge = 0;
            if (sng.m_corepoint[vertex2] == 1)
                should_merge = 1;
            else if (sng.m_member[vertex2] == 0) {
                omp_set_lock(&node_locks[vertex2]);
                if (sng.m_member[vertex2] == 0) {
                    should_merge = 1;
                    sng.m_member[vertex2] = 1;
                }
                omp_unset_lock(&node_locks[vertex2]);
            }

            if (should_merge == 1) {
                // Union-Find based approach for merging
                root1 = vertex1;
                root2 = vertex2;

                // Union-Find algorithm with splicing compression techniques
                while (sng.m_parents[root1] != sng.m_parents[root2]) {
                    if (sng.m_parents[root1] < sng.m_parents[root2]) {
                        if (sng.m_parents[root1] == root1) {
                            omp_set_lock(&node_locks[root1]);
                            int is_parent_set = 0;
                            if (sng.m_parents[root1] == root1) {
                                sng.m_parents[root1] = sng.m_parents[root2];
                                is_parent_set = 1;
                            }
                            omp_unset_lock(&node_locks[root1]);
                            if (is_parent_set) // Merge successful
                                break;
                        }

                        // Splicing
                        int temp = sng.m_parents[root1];
                        sng.m_parents[root1] = sng.m_parents[root2];
                        root1 = temp;
                    } else {
                        if (sng.m_parents[root2] == root2) {
                            omp_set_lock(&node_locks[root2]);
                            int is_parent_set = 0;
                            if (sng.m_parents[root2] == root2) {
                                sng.m_parents[root2] = sng.m_parents[root1];
                                is_parent_set = 1;
                            }
                            omp_unset_lock(&node_locks[root2]);
                            if (is_parent_set) // Merge Successful
                                break;
                        }

                        // Splicing
                        int temp = sng.m_parents[root2];
                        sng.m_parents[root2] = sng.m_parents[root1];
                        root2 = temp;
                    }
                }
            }
        }
    }
    stop_time = omp_get_wtime();
    free(node_locks);
    cout << "Merging took " << stop_time - start_time << " seconds." << endl;

    for (thread_id = 0; thread_id < max_threads; thread_id++)
        merge[thread_id].clear();

    merge.clear();
    neighbors.clear();
}









// Run the Sequential Sow & Grow (SNG) Clustering Algorithm
void run_sng_algo(ClusteringAlgo& sng) {
		
	cout << endl; 
	cout << "SNG SEQUENTIAL ALGORITHM" << endl;
	cout << sng.m_seeds << endl;
		
	// Select m_seeds random points
	vector<int> random_seeds;
	random_seeds.reserve(sng.m_seeds);

	// Check if there are more points than seeds
    if (sng.m_seeds >= sng.m_pts->m_i_num_points) {
        cout << "Error: Number of seeds is greater than or equal to the number of points. Aborting." << endl;
        exit(-1);
    }

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


		
	int i, pid, j, k, npid;
	int cid = 1; // cluster id
	vector <int> c;
	c.reserve(sng.m_pts->m_i_num_points);

   	// Initialize clustering parameters
	sng.m_noise.resize(sng.m_pts->m_i_num_points, false);
    sng.m_visited.resize(sng.m_pts->m_i_num_points, false);		
	sng.m_pid_to_cid.resize(sng.m_pts->m_i_num_points, 0);
	sng.m_clusters.clear();

		

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