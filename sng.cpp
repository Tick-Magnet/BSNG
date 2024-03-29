#include "sng.h"
#include <thread>
#include <chrono>

void displayScatterPlot2D(NWUClustering::ClusteringAlgo * algorithm, vector<int> * clusters);
void displayScatterPlot2D(NWUClustering::ClusteringAlgo * algorithm);


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
            o << " Value: (";
			//Outputing datapoint 
			for(int j = 0; j < m_pts->m_i_dims; j++)
			{
				o << m_pts->m_points[i][j] << ", ";
			}
			o << ")" << endl;
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
		
        const string directoryPath = "utilities/";
	    string fullPath = directoryPath + csvOutputFilename;	
		if(csvOutputFilename != NULL)
		{
			cout << "Writing csv file\n";

			//Open csv file
			ofstream csvFile;
            csvFile.open(fullPath);
            //Add each point in data set
            //First column is cluster label
            for(int i = 0; i < m_pts->m_i_num_points; i++)
            {
				int clusterID = m_pid_to_cid[i];
				csvFile << clusterID << ',';
				for(int j = 0; j < m_pts->m_i_dims; j++)
				{
					csvFile << m_pts->m_points[i][j];
					if(j < m_pts->m_i_dims - 1)
						csvFile << ',';
				}
				
				csvFile << '\n';
			}
		}
		
	}

	void ClusteringAlgo::writeClusters_uf(ostream& o)
	{
		// Writing point id and cluster id pairs per line; noise has cluster id 0	
		vector <int> clusters;
		clusters.resize(m_pts->m_i_num_points, 0);

		int i, j, sum_points = 0, noise = 0, root, rootcount = 0, tmp;

        o << "- - - Parallel SNG Clustering Output - - -" << endl;
        o << "Key: 0 = Noise | > 0 = Cluster ID" << endl;
        o << endl;
        o << "ID" << " | " << "Cluster" << endl; 

		// Calculate cluster information
		for(i = 0; i < m_pts->m_i_num_points; i++)
		{
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
		for (i = 0; i < m_pts->m_i_num_points; i++) 
		{
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
            //Formatting for Output
            if (i < 10){
                o << i << "  | " << clusters[m_parents[i]] << endl;
                o << " Value: (";
				//Outputing datapoint 
				for(int j = 0; j < m_pts->m_i_dims; j++)
				{
					o << m_pts->m_points[i][j] << ", ";
				}
				o << ")" << endl;
            } else {
                o << i << " | " << clusters[m_parents[i]] << endl;
                o << " Value: (";
				//Outputing datapoint 
				for(int j = 0; j < m_pts->m_i_dims; j++)
				{
					o << m_pts->m_points[i][j] << ", ";
				}
				o << ")" << endl;
            }
		}
        
		// Output summary information
        o << endl;
        o << "Total points " << noise + sum_points << " pt_in_cls " << sum_points << " noise " << noise << endl;
        o << "Number of clusters: " << count << endl;

        cout << "Total points " << noise + sum_points << " pt_in_cls " << sum_points << " noise " << noise << endl;
		cout << "Number of clusters: " << count << endl;
		
	const string directoryPath = "utilities/";
	string fullPath = directoryPath + csvOutputFilename;		
		if(csvOutputFilename != NULL)
		{
			cout << "Writing csv file\n";
			
			//Open csv file
			ofstream csvFile(fullPath);
      
            //Add each point in data set
            //First column is cluster label
            for(int i = 0; i < m_pts->m_i_num_points; i++)
            {

				int clusterID = clusters[m_parents[i]];
				csvFile << clusterID << ',';
				for(int j = 0; j < m_pts->m_i_dims; j++)
				{
					csvFile << m_pts->m_points[i][j];
					if(j < m_pts->m_i_dims - 1)
						csvFile << ',';
				}
				
				csvFile << '\n';
			}
		}
		clusters.clear();
	}

    

    void seed_selection(ClusteringAlgo& sng) {
        
        cout << endl;
        cout << "Starting Seed Selection" << endl;
        cout << "Seed Method: " << sng.m_seedMethod << " | ";
	    double start = omp_get_wtime() ;	

        //TODO: Add other Seed Selection Methodologies 
        //Sections (X-Y) | 1st Fourth

        switch (sng.m_seedMethod) {

            case 0: //Random Seed Selection - m_seeds random points
                {
                    cout << "Random Seed Selection" << endl;
                    sng.selected_seeds.reserve(sng.m_seeds);

                    // Check if there are more points than seeds
                    if (sng.m_seeds > sng.m_pts->m_i_num_points) {
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
                            //cout << random_index << " "; 
                        }
                    }

                    // cout << endl; // Print a newline to end seed list
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

                    sng.m_seeds = half; // This is so that sequential the system will loop through all seeds.

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
                        sng.selected_seeds.push_back(seed_sequence[i]); vector<int> selected_seeds;
                        //cout << seed_sequence[i] << " ";
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
                    std::iota(seed_sequence.begin(), seed_sequence.end(), 0); 

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
            case 5: 
                {
                    vector<int> selected_seeds;
                    break;
                }

            default: 
                {
                    cout << "    -z seedMethod   : select method for seed selection (default 0)\n"
                            "       * 0 random   : random seed values\n"
                            "       * 1 even     : even seed values\n"
                            "       * 2 odd      : odd seed values\n"
                            "       * 3 Lower    : all under input seed\n"
                            "       * 4 Upper    : all over input seed\n" << endl;
                    exit(1);
                }
        }

        double stop = omp_get_wtime();
        cout << "Seed Selection took " << stop - start << " seconds." << endl;
        cout << endl;
    }


// Run the parallel version of the Sow-and-Grow (SNG) clustering algorithm
void run_sng_algo_uf(ClusteringAlgo& sng) {	
   
    // Select Seeds for Algorithm
    seed_selection(sng);

    // Initialize clustering parameters
    int thread_id, i, j, k, neighbor_point_id, point_id, root1, root2;
    int max_threads = omp_get_max_threads();
    
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
    
    double start_time = omp_get_wtime();

    // Reserve space for the merge vector
    for (i = 0; i < max_threads; i++) 
    {
        merge[i].reserve(sng.m_pts->m_i_num_points);
    }
    
    // Initialize stack for point expansion
    vector<stack<int>> pointStacks(max_threads); // Create a stack of vectors, one for each thread

    vector<int> prID; // Point ID
    prID.resize(sng.m_pts->m_i_num_points, -1); // Assign all points to thread: -1
    
    // Initilize load balancing values to -1
    int loadBalancingRequests[max_threads];
    for(int i = 0; i < max_threads; i++)
	{
		loadBalancingRequests[i] = -1;
	}
    
    // Lock for loadBalancingRequests array
	omp_lock_t *loadBalancingLock =(omp_lock_t*) malloc(sizeof(omp_lock_t));
	omp_init_lock(loadBalancingLock);
	

    #pragma omp parallel private(root1, root2, thread_id, neighbors, neighbor_point_id, i, j, point_id) shared(pointStacks, sng, prID, loadBalancingLock, loadBalancingRequests)
    {
        
        // Get Thread ID
        thread_id = omp_get_thread_num();
        stack<int>& pointStack = pointStacks[thread_id];


        // Distribute seeds to threads without overlap
        int seeds_per_thread = sng.m_seeds / max_threads;
        int extra_seeds = sng.m_seeds % max_threads;  // Distribute remaining seeds
        int start_seed = thread_id * seeds_per_thread + min(thread_id, extra_seeds); // First Seed for Thread
        int end_seed = (thread_id + 1) * seeds_per_thread + min(thread_id + 1, extra_seeds); // Last Seed for Thread

        // Assigning Seeds to Stack
        for (int seed_idx = start_seed; seed_idx < end_seed; ++seed_idx) {
            pointStack.push(sng.selected_seeds[seed_idx]);
            prID[seed_idx] = thread_id;
        }


        // Assigning Point's Parent to Itself
        int points_per_thread = sng.m_pts->m_i_num_points / max_threads;
        int extra_points = sng.m_pts->m_i_num_points % max_threads;
        int start_point = thread_id * points_per_thread + min(thread_id, extra_points);
        int end_point = (thread_id + 1) * points_per_thread + min(thread_id + 1, extra_points);

        // Loop through every point assign default parent
        for (int point_id = 0; point_id < sng.m_pts->m_i_num_points; ++point_id) {
            sng.m_parents[point_id] = point_id;
        }

        // Wait for all threads to claim seeds and assign point parents.
        #pragma omp barrier

        // Primary Thread Stack
        while (!pointStack.empty()) 
        {
			if(loadBalancingRequests[thread_id] != -1)
			{
				int callingThread = loadBalancingRequests[thread_id];
				//cout << thread_id << " sending points to " << callingThread << endl;
				//Grant calling thread half of point stack
				int numberToSend = pointStack.size() / 2;
				int tempPoint;
				for(int i = 0; i < numberToSend; i++)
				{
					//cout << "Break\n";
					//take points from back of stack
					tempPoint = pointStack.top();
					pointStack.pop();
					//cout << "Break2\n";

					//Reassign prID values to new thread
					prID[tempPoint] = callingThread;
					
					//Push point onto calling thread's stack
					pointStacks[callingThread].push(tempPoint);
				}
				//cout << "grabbing lock\n";
				//Flip flag back to -1
				omp_set_lock(loadBalancingLock);
					//cout << "flipping\n";
					//Flipping requesting threads flag back to -1
					loadBalancingRequests[callingThread] = -1;
					loadBalancingRequests[thread_id] = -1;
				omp_unset_lock(loadBalancingLock);
			}
			
            int currentPoint = pointStack.top();
            pointStack.pop();
           
            // Find neighborhood for current point.
            neighbors.clear(); 
            sng.m_kdtree->r_nearest_around_point(currentPoint, 0, sng.m_epsSquare, neighbors);

            // If Core Point
            if (neighbors.size() >= sng.m_minPts) 
            {
                sng.m_corepoint[currentPoint] = 1;
                sng.m_member[currentPoint] = 1;

                for (j = 0; j < neighbors.size(); j++) 
                {
                    neighbor_point_id = neighbors[j].idx;

                    //Flag to be Merged Across Partitions
                    if(prID[neighbor_point_id] != -1 && prID[neighbor_point_id] != thread_id)
                    {
                        merge[thread_id].push_back(currentPoint);
                        merge[thread_id].push_back(neighbor_point_id);
                        continue;    
                    } 

                    //Merge Neighborhood
                    else 
                    {                        

                        if (prID[neighbor_point_id] == -1) {
                            //Adding Point to Queue
                            pointStack.push(neighbor_point_id);
                            prID[neighbor_point_id] = thread_id;
                        }
                        
                        //Change Roots
                        if(sng.m_corepoint[neighbor_point_id] == 1 || sng.m_member[neighbor_point_id] == 0) // Border Point Check
                        {
                            // Get the root containing currentPoint
                            root1 = neighbor_point_id;
                            root2 = currentPoint;

                            sng.m_member[neighbor_point_id] = 1;

                            // Union-Find algorithm to merge the trees
                            while (sng.m_parents[root1] != sng.m_parents[root2]) 
                            {
                                if (sng.m_parents[root1] < sng.m_parents[root2]) //If Neighbor is less than Current Point
                                {
                                    if (sng.m_parents[root1] == root1) {
                                        sng.m_parents[root1] = sng.m_parents[root2];
                                        //root = dbs.m_parents[root1];
                                        break;
                                    }

                                    // Splicing 
                                    int temp = sng.m_parents[root1];
                                    sng.m_parents[root1] = sng.m_parents[root2];
                                    root1 = temp;

                                } 
                                else //If Neighbor is greater than Current Point
                                {
                                    if (sng.m_parents[root2] == root2) {
                                        sng.m_parents[root2] = sng.m_parents[root1];
                                        //root = dbs.m_parents[root1];
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
            if(pointStack.empty())
            {
				// Grab other data points if available
				for(int i = 0; i < max_threads; i++)
				{
					omp_set_lock(loadBalancingLock);

                    int threshold = 8; // Should be reviewed for optimization 

					if(thread_id != i && pointStacks[i].size() >= threshold && loadBalancingRequests[i] == -1)
					{
						loadBalancingRequests[i] = thread_id;
						loadBalancingRequests[thread_id] = -2;
						omp_unset_lock(loadBalancingLock);
						// Busy wait until points granted from other thread
						while(loadBalancingRequests[thread_id] == -2)
						{
							std::this_thread::sleep_for(std::chrono::microseconds(1));
						}
						break;
					}
					else
					{
						omp_unset_lock(loadBalancingLock);
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

    /*
    for (int i = 0; i < sng.m_pts->m_i_num_points; ++i) {
        //cout << i << ": Member " << sng.m_member[i] << " | Corepoint: " << sng.m_corepoint[i] << endl;
		cout << "Point " << i << ": Parent = " << sng.m_parents[i] << endl;
	} cout << endl;
    */

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

    cout << "Begin Sequential - Sow & Grow (S-SNG)" << endl;
	double start = omp_get_wtime() ;		

	// Iterate through points
	for (int i = 0; i < sng.m_seeds; i++) {
       	
        int pid = sng.selected_seeds[i];
        //cout << "pid: " << pid << endl;
		
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
	cout << "No merging stage in Sequential SNG"<< endl;
	ind = NULL;
	ne.clear();
	ne2.clear();
}

};
