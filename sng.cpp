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

		cout << "SNG Parallel ALGORITHM" << endl;
		cout << sng.m_seeds << endl;
		cout << sng.m_pts->m_i_num_points << endl;
		
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



		int thread_id, i, point_id, j, k, neighbor_point_id, root, root1, root2;

		// Initialize clustering parameters
		sng.m_clusters.clear();
		kdtree2_result_vector neighbors;

		// Assign a parent to each test point
		sng.m_parents.resize(sng.m_pts->m_i_num_points, -1);
		sng.m_member.resize(sng.m_pts->m_i_num_points, 0);
		sng.m_corepoint.resize(sng.m_pts->m_i_num_points, 0);

		int chunk_size, max_threads = omp_get_max_threads();

		// Calculate the thread distribution
		//Threads Need to be 
		if (sng.m_pts->m_i_num_points % max_threads == 0)
			chunk_size = sng.m_pts->m_i_num_points / max_threads;
		else
			chunk_size = sng.m_pts->m_i_num_points / max_threads + 1;

		vector<vector<int>> merge;
		vector<int> init;
		merge.resize(max_threads, init);

		// Reserve space for the merge vector
		for (thread_id = 0; thread_id < max_threads; thread_id++)
			merge[thread_id].reserve(sng.m_pts->m_i_num_points);

		vector<int> thread_id_map;
		thread_id_map.resize(sng.m_pts->m_i_num_points, -1);

		//The point_indices needs to be Random Points
		vector<int>* point_indices = sng.m_kdtree->getIndex();

		double start_time = omp_get_wtime();

		#pragma omp parallel private(root, root1, root2, thread_id, neighbors, neighbor_point_id, i, j, point_id) shared(chunk_size, point_indices)
		{
			int lower, upper;
			thread_id = omp_get_thread_num();

			lower = chunk_size * thread_id;
			upper = chunk_size * (thread_id + 1);

			if (upper > sng.m_pts->m_i_num_points)
				upper = sng.m_pts->m_i_num_points;

			for (i = lower; i < upper; i++) {
				point_id = (*point_indices)[i];
				sng.m_parents[point_id] = point_id;
				thread_id_map[point_id] = thread_id;
			}

			#pragma omp barrier

			// Compute core points and neighbors
			for (i = lower; i < upper; i++) {
				point_id = (*point_indices)[i];

				neighbors.clear();
				sng.m_kdtree->r_nearest_around_point(point_id, 0, sng.m_epsSquare, neighbors);

				if (neighbors.size() >= sng.m_minPts) {
					sng.m_corepoint[point_id] = 1;
					sng.m_member[point_id] = 1;

					// Get the root containing point_id
					root = point_id;

					for (j = 0; j < neighbors.size(); j++) {
						neighbor_point_id = neighbors[j].idx;
						if (thread_id_map[neighbor_point_id] != thread_id) {
							merge[thread_id].push_back(point_id);
							merge[thread_id].push_back(neighbor_point_id);
							continue;
						}

						// Get the root containing neighbor_point_id
						root1 = neighbor_point_id;
						root2 = root;

						if (sng.m_corepoint[neighbor_point_id] == 1 || sng.m_member[neighbor_point_id] == 0) {
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
	











	// Run the Sequential SNG clustering algorithm
	void run_sng_algo(ClusteringAlgo& sng) {
		
		
		cout << "SNG SEQUENTIAL ALGORITHM" << endl;
		cout << sng.m_seeds << endl;
		
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