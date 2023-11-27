#include "sng.h"


namespace NWUClustering
{
	// Set SNG algorithm parameters
    void ClusteringAlgo::set_sng_params(double eps, int minPts, int seeds, int seedMethod) {
		m_epsSquare =  eps * eps;
		m_minPts =  minPts;
		m_seeds = seeds;
        m_seedMethod = seedMethod;
	}

	// Destructor to clean up resources
	ClusteringAlgo::~ClusteringAlgo() {
		m_noise.clear();
		m_visited.clear();
		m_parents.clear();
		m_corepoint.clear();
		m_member.clear();
        selected_seeds.clear();
	}
	

    //TODO: Add Seeds to Output file!


	// Write clusters to the output stream
	void ClusteringAlgo::writeClusters(ostream& o) {
		
        // Writing point id and cluster id pairs per line; noise has cluster id 0	
		int iMaxID = m_clusters.size(), id, i;

		int sum_points = 0;
		int noise = 0;
        int unclustered = 0;

        o << "- - - Sequential SNG Clustering Output - - -" << endl;
        o << "Key: -1 = Noise | 0 = Unclustered | > 0 = Cluster ID" << endl;
        o << endl;
        o << "ID" << " | " << "Cluster" << endl; 

		for(i = 0; i < m_pts->m_i_num_points; i++) {

			id = m_pid_to_cid[i];
            if (id == 0 && m_visited[i]) {
                id = -1;
            }

            //Formatting for Output
            if (i < 10){
                o << " " << i << " | " << id << endl;
            } else {
                o << i << " | " << id << endl;
            }
		}

		for(i = 0; i < m_clusters.size(); i++) {
			sum_points += m_clusters[i].size();
		}
	
		for (i = 0; i < m_pts->m_i_num_points; i++) {
			if(m_noise[i]) {
                noise++;
            }	
		}	

        for (int i = 0; i < m_pts->m_i_num_points; i++) {
            if (!m_visited[i]) {
                unclustered++;
            }
        }
    

        cout << "Total points " << unclustered + noise + sum_points << " | pts_cls " << sum_points << " | noise " << noise << " | pts_uncls " << unclustered << endl;
		cout << "Number of clusters: " << m_clusters.size() << endl;

		// Output summary information
        o << endl;
        o << "Total points " << unclustered + noise + sum_points << " | pts_cls " << sum_points << " | noise " << noise << " | pts_uncls " << unclustered << endl;
        o << "Number of clusters: " << m_clusters.size() << endl;

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

        o << "- - - Sequential SNG Clustering Output - - -" << endl;
        o << "Key: 0 = Noise | -1 = Unclustered | > 0 = Cluster ID" << endl;
        o << " " << endl;
        o << "ID" << " | " << "Cluster" << endl; 
        
		// Write point id and cluster ids to the output stream
		for (i = 0; i < m_pts->m_i_num_points; i++) {
			o << i << " | " << clusters[m_parents[i]] << endl;
		}
        
		// Output summary information
        o << endl;
        o << "Total points " << noise + sum_points << " pt_in_cls " << sum_points << " noise " << noise << endl;
        o << "Number of clusters: " << count << endl;

        cout << "Total points " << noise + sum_points << " pt_in_cls " << sum_points << " noise " << noise << endl;
		cout << "Number of clusters: " << count << endl;

		clusters.clear();
	}

    

    void seed_selection(ClusteringAlgo& sng) {
        
        cout << endl;
        cout << "Seed Selection | # of Seeds: " << sng.m_seeds << endl;
        cout << "Seed Method: " << sng.m_seedMethod << endl;
        //TODO: Add other Seed Selection Methodologies 
        //Sections (X-Y) | 1st Fourth

        switch (sng.m_seedMethod) {

            case 0: //Random Seed Selection - m_seeds random points
                {
                    cout << "Random Seed Selection" << endl;
                    sng.selected_seeds.reserve(sng.m_seeds);

                    // Check if there are more points than seeds
                    if (sng.m_seeds >= sng.m_pts->m_i_num_points) {
                        cout << "Error: Number of seeds is greater than or equal to the number of points. Aborting." << endl;
                        exit(-1);
                    }

                    std::random_device rd;
                    std::mt19937 gen(rd());
                    std::uniform_int_distribution<int> dis(0, sng.m_pts->m_i_num_points - 1);

                    cout << "Selected random point(s): "; 

                    while (sng.selected_seeds.size() < sng.m_seeds) {
                        int random_index = dis(gen);

                        //Checks if point hasn't been selected previously.
                        if (std::find(sng.selected_seeds.begin(), sng.selected_seeds.end(), random_index) == sng.selected_seeds.end()) {
                            sng.selected_seeds.push_back(random_index);
                            cout << random_index << " "; 
                        }
                    }

                    cout << endl; // Print a newline to end seed list
                    cout << endl;
                    break;
                }
            case 1: 
                {
                    cout << "Odd Seed Selection" << endl;
                    int half = sng.m_pts->m_i_num_points / 2;
                    sng.selected_seeds.reserve(half);

                    // Use a fixed sequence starting from 1 as the seed
                    std::vector<int> seed_sequence(sng.m_pts->m_i_num_points);
                    std::iota(seed_sequence.begin(), seed_sequence.end(), 1); // Fills the sequence with 1, 2, 3, ..., n

                    // Select every odd number as a seed
                    for (int i = 0; i < half; ++i) {
                        int odd_seed = seed_sequence[i * 2]; // Select every odd index
                        sng.selected_seeds.push_back(odd_seed);
                        cout << odd_seed << " ";
                    }

                    sng.m_seeds = half; //This is so that sequential the system will loop through all seeds.

                    cout << endl; // Print a newline to end seed list
                    cout << endl;
                    break;
                }
            case 2: 
                {
                    cout << "Even Seed Selection" << endl;
                    int half = sng.m_pts->m_i_num_points / 2;
                    sng.selected_seeds.reserve(half);

                    // Use a fixed sequence starting from 0 as the seed
                    std::vector<int> seed_sequence(sng.m_pts->m_i_num_points);
                    std::iota(seed_sequence.begin(), seed_sequence.end(), 0); // Fills the sequence with 0, 1, 2, ..., n-1

                    // Select every even number as a seed
                    for (int i = 0; i < half; ++i) {
                        int even_seed = seed_sequence[i * 2]; // Select every even index
                        if (even_seed != half * 2) { // Exclude half when it's even
                            sng.selected_seeds.push_back(even_seed);
                            cout << even_seed << " ";
                        }
                    }

                    sng.m_seeds = sng.selected_seeds.size(); // Update m_seeds to the actual number of selected seeds

                    cout << endl; // Print a newline to end seed list
                    cout << endl;
                    break;

                }
            case 3: 
                {

                    //TODO: Add Error Checking

                    cout << "Lower Partition Seed Selection" << endl;
                    int num = sng.m_seeds;
                    sng.selected_seeds.reserve(num);

                    // Use a fixed sequence starting from 0 as the seed
                    std::vector<int> seed_sequence(sng.m_pts->m_i_num_points);
                    std::iota(seed_sequence.begin(), seed_sequence.end(), 0); // Fills the sequence with 0, 1, 2, ..., n-1

                    // Add every point below sng.m_seeds to sng.selected_seeds
                    for (int i = 0; i < sng.m_seeds; ++i) {
                        sng.selected_seeds.push_back(seed_sequence[i]); 
                        cout << seed_sequence[i] << " ";
                    }

                    sng.m_seeds = sng.selected_seeds.size(); // Update m_seeds to the actual number of selected seeds
                    cout << endl; // Print a newline to end seed list
                    cout << endl;
                    break;
                }
            case 4: 
                {

                    //TODO: Add Error Checking

                    cout << "Upper Partition Seed Selection" << endl;
                    int num = sng.m_pts->m_i_num_points - sng.m_seeds;
                    sng.selected_seeds.reserve(num);

                    // Use a fixed sequence starting from 0 as the seed
                    std::vector<int> seed_sequence(sng.m_pts->m_i_num_points);
                    std::iota(seed_sequence.begin(), seed_sequence.end(), 0); // Fills the sequence with 0, 1, 2, ..., n-1

                    // Add every point below sng.m_seeds to sng.selected_seeds
                    for (int i = sng.m_pts->m_i_num_points -1; i >= sng.m_seeds; --i) {
                        sng.selected_seeds.push_back(seed_sequence[i]); 
                        cout << seed_sequence[i] << " ";
                    }

                    sng.m_seeds = sng.selected_seeds.size(); // Update m_seeds to the actual number of selected seeds
                    cout << endl; // Print a newline to end seed list
                    cout << endl;
                    break;
                }
            default: 
                {
                    cout << "Not a proper seed method." << endl;
                    exit(1);
                }
        }
    }


// Run the Union-Find version of the Sow-and-Grow (SNG) clustering algorithm
void run_sng_algo_uf(ClusteringAlgo& sng) {	
   
    //Select Seeds for Algorithm
	//seed_selection(sng);

    // Add two numbers to the selected_seeds vector
    sng.selected_seeds.reserve(sng.m_seeds);
    sng.selected_seeds.push_back(9);
    sng.selected_seeds.push_back(11);

    int thread_id, i, j, k, neighbor_point_id, point_id, root, root1, root2;
    int max_threads = omp_get_max_threads();

    // Initialize clustering parameters
    sng.m_clusters.clear();
    kdtree2_result_vector neighbors;

    sng.m_parents.resize(sng.m_pts->m_i_num_points, -1);
    sng.m_member.resize(sng.m_pts->m_i_num_points, 0);
    sng.m_corepoint.resize(sng.m_pts->m_i_num_points, 0);


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


    vector<int> prID; 
    prID.resize(sng.m_pts->m_i_num_points, -1);


    #pragma omp parallel private(root, root1, root2, thread_id, neighbors, neighbor_point_id, i, j, point_id) shared(pointQueues, sng, prID)
    {

        //Get Thread ID
        thread_id = omp_get_thread_num();
        queue<int>& pointQueue = pointQueues[thread_id];

        // Distribute seeds to threads without overlap
        int seeds_per_thread = sng.m_seeds / max_threads;
        int start_seed = thread_id * seeds_per_thread;
        int end_seed = (thread_id == max_threads - 1) ? sng.m_seeds : start_seed + seeds_per_thread;

        #pragma omp critical
        {
            cout << "Thread " << thread_id << " processes seeds from " << start_seed << " to " << end_seed - 1 << endl;
            cout << "Size of prID vector: " << prID.size() << endl;
        }


        for (int seed_idx = start_seed; seed_idx < end_seed; ++seed_idx) {
            pointQueue.push(sng.selected_seeds[seed_idx]);
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

            prID[currentPoint] = thread_id;
    
            // Find core points and neighbors for currentPoint
            neighbors.clear();
            sng.m_kdtree->r_nearest_around_point(currentPoint, 0, sng.m_epsSquare, neighbors);

            if (neighbors.size() >= sng.m_minPts && sng.m_corepoint[currentPoint] != 1) {
                

				sng.m_corepoint[currentPoint] = 1;
				sng.m_member[currentPoint] = 1;
                // Get the root containing currentPoint
                root = currentPoint;
            

                //TODO: Sort Neighbor points by highest.

                for (j = 0; j < neighbors.size(); j++) {
                    
					neighbor_point_id = neighbors[j].idx;

                    //cout << "Requesting to Merge | Neighborhood Point Thread ID: " << prID[neighbor_point_id] << " | Current Point Thread ID: " << prID[neighbor_point_id] << endl;
                    //If point Thread ID isn't in this thread THEN threads should merge.
					if(prID[neighbor_point_id] != -1 &&  prID[neighbor_point_id] != prID[currentPoint]) {
						
                        //cout << "Requesting to Merge | Neighborhood Point Thread ID: " << prID[neighbor_point_id] << " | Current Point Thread ID: " << prID[currentPoint] << endl;
                        //cout << "Neighborhood Point: " << neighbor_point_id << " | Current Point: " << currentPoint << endl;
                        merge[thread_id].push_back(neighbor_point_id);
                        merge[thread_id].push_back(currentPoint);
                        continue; 
				    }

                    root1 = neighbor_point_id;
                    root2 = root;

                    //cout << "In Neighborhood " << neighbor_point_id  << " membership: " <<sng.m_member[neighbor_point_id] << " | Current Point: " <<  currentPoint <<" | Thread ID: " << thread_id << endl;

                    if (sng.m_member[neighbor_point_id] == 0) {
                        
                        sng.m_member[neighbor_point_id] = 1;

                        // Union-Find algorithm to merge the trees
                        while (sng.m_parents[root1] != sng.m_parents[root2]) {
                            if (sng.m_parents[root1] < sng.m_parents[root2]) {
                                if (sng.m_parents[root1] == root1) {
                                    sng.m_parents[root1] = sng.m_parents[root2];
                                    root = sng.m_parents[root2];
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
        cout << merge_size << endl;
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
                //cout << endl;
                //cout << "SHOULD_MERGE!" << endl;
                // Union-Find based approach for merging
                root1 = vertex1;
                root2 = vertex2;

                //cout << "Root1: " << root1 << endl;
                //cout << "Root2: " << root2 << endl;

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
                //cout << "Root1: " << root1 << endl;
                //cout << "Root2: " << root2 << endl;
            }
        }
    }

    for (int i = 0; i < sng.m_pts->m_i_num_points; ++i) {
		cout << "Point " << i << ": Parent = " << sng.m_parents[i] << endl;
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
	
    //Select Seeds for Algorithm
	seed_selection(sng);
	
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

    cout << "Start" << endl;
	double start = omp_get_wtime() ;		

	// Iterate through points
	for (int i = 0; i < sng.m_seeds; i++) {
       	
        int pid = sng.selected_seeds[i];
        cout << "pid: " << pid << endl;
		
        if (!sng.m_visited[pid]) {
			sng.m_visited[pid] = true;
			ne.clear();
			sng.m_kdtree->r_nearest_around_point(pid, 0, sng.m_epsSquare, ne);
			
            //If neighborhood size is too small, assign as noise.
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