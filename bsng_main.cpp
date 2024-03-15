#include "sng.h"
#include "dbscan.h"
#include "utils.h"
#include "kdtree2.hpp"
#include <cstdlib>

// Function to Display System Usage 
static void usage(char* argv0) {
    
    const char* params =
    "\n"
    "Usage: %s [flags] -i filename -m minpts -e epsilon -t threads -s seeds -o output \n"
	"Example: ./bsng -i clus50k.bin -b -m 5 -e 25 -s 5000 -t 2 -o output\n"
        "\n"
        "    -b isBinary     : input file is in binary format (default no)\n"
        "    -d dbscan       : run dbscan instead of sow-and-grow (default no)\n"
        "    -e epsilon      : input parameter of BSNG, radius of neighborhood, e.g. 5\n"
        "    -i filename     : file containing input data to be clustered\n"
        "    -m minpts       : input parameter of BSNG, min points to form a cluster, e.g. 3\n"
        "    -o output       : clustering results output file - format (point id, clusterid)\n"
        "    -s seeds        : input parameter of BSNG, input seed(s) (default 1)\n"
        "    -t threads      : number of threads to be employed (default 1)\n"      
        "    -z seedMethod   : select method for seed selection (default 0)\n"
        "       * 0 random   : random seed values\n"
        "       * 1 even     : even seed values\n"
        "       * 2 odd      : odd seed values\n"
        "       * 3 Lower    : all under input seed\n"
        "       * 4 Upper    : all over input seed\n"
        "       * 5 UserList : read from userList\n"
        "    -u userList     : Select the seeds to be used from a list\n"
	    "    -j visualize    : visualize data\n"
        "\n";

    fprintf(stderr, params, argv0);
    exit(-1);
}

// Main Loop
int main(int argc, char** argv) {

    double seconds, eps;
    int minPts, threads, opt, seeds, method;
    char* outfilename = NULL;
    char* infilename = NULL;
    char* csvOutputFilename = "visData.csv";
    char* userfilename = NULL;
    bool classical, isBinaryFile, dbscan;
    bool visualize;
    string seedInputFileName;

    // Initialize Default Values
    minPts = -1;
    eps = -1;
    isBinaryFile = 0;
    outfilename = NULL;
    infilename = NULL;
    userfilename = NULL;
    threads = 1;
    seeds = 1;
    method = 0;
    classical = false;
    isBinaryFile = false;
    dbscan = false;
    visualize = false;

    // Input Flags
    while ((opt = getopt(argc, argv, "i:u:t:p:m:e:s:o:s:z:l:bdxghncuj")) != EOF) {
        switch (opt) {
			case 'l':
				seedInputFileName = string(optarg);
				break;
            case 'i':
                infilename = optarg;
                break;
            case 'u':
                userfilename = optarg;
                break;
            case 't':
                threads = atoi(optarg);
                break;
            case 'm':
                minPts = atoi(optarg);
                break;
            case 'e':
                eps = atof(optarg);
                break;
            case 's':
                seeds = atoi(optarg);
                break;
            case 'z':
                method = atoi(optarg);
                break;
            case 'o':
                outfilename = optarg;
                break;
            case 'c':
                classical = true;
                break;
            case 'd':
                dbscan = true;
                break;
            case 'b':
                isBinaryFile = true;
                break;
	    case 'j':
		visualize = true;
		break;
            case '?':
                usage(argv[0]);
                break;
            default:
                usage(argv[0]);
                break;
        }
        //printf("opt: %c, optarg: %s\n", opt, optarg); //Argument Debugging
    }


    /* 
    Sow & Grow Algorithm - SNG (Runs by Default)
    */

    if (dbscan == false) 
    {

        // Check for provided input arguments
        if (infilename == NULL || minPts < 0 || seeds < 0 || eps < 0 || threads < 0) {
            usage(argv[0]);
            exit(-1);
        }

        cout << endl;
        cout << "|| SOW & GROW ||    " << endl;
        cout << endl;

        if (threads > seeds) {
            threads = seeds; 
            cout << "User provided more threads than seeds. Removing excess threads, remaining: " << threads << endl;
        }
        

        // Set the number of OpenMP threads
        omp_set_num_threads(threads);

        // Create an instance of the ClusteringAlgo class
        NWUClustering::ClusteringAlgo sng;
        sng.set_sng_params(eps, minPts, seeds, method);
        sng.seedInputFileName = seedInputFileName;
		sng.csvOutputFilename = csvOutputFilename;
        cout << "Input parameters " << " minPts: " << minPts << ", eps: " << eps << ", seeds: " << seeds << endl;

        // Measure the time to Read File Data
        cout << "Reading points from file: " << infilename << endl;
        double start = omp_get_wtime();

        // Read input data from the file
        if (sng.read_file(infilename, isBinaryFile) == -1) {
            cout << "File Reading Failed." << endl;
            exit(-1);
        }
        cout << "Reading input data file took: " << omp_get_wtime() - start << " seconds." << endl;

        // Build a k-d tree for the points
        cout << "Begin Building kdtree . . ." << endl;
        start = omp_get_wtime();
        sng.build_kdtree();
        cout << "Building kdtree took: " << omp_get_wtime() - start << " seconds." << endl;

        // Run the SNG clustering algorithms
        start = omp_get_wtime();
        cout << endl; 

        //Selecting Algorithm 
        if (classical == true) {
            cout << "- SNG SEQUENTIAL ALGORITHM -" << endl;
            run_sng_algo(sng); // Sequential 
        } else {
            cout << "- SNG Parallel ALGORITHM -" << endl;
            run_sng_algo_uf(sng); // Parallel
        }
            
        cout << "S&G (total) took " << omp_get_wtime() - start << " seconds." << endl;

        // If an output fiilename is provided, write clustering results to the file
        if (outfilename != NULL) {
            ofstream outfile;
            outfile.open(outfilename);

                if (classical == true) {
                sng.writeClusters(outfile); // Sequential 
            } else {
                sng.writeClusters_uf(outfile); // Parallel
            }

            outfile.close();

            // Running visualiztion code
            if(visualize == true){
                cout << "Beginning Visualization\n";

                // Declaring variables used in command for visualization
                const char* cdcommand = "cd ./utilities";
                const char* pycommand = "python3 ClusterVisualizer.py -i ./visData.csv";
                std::string command = std::string(cdcommand) + " && " + pycommand;
                
                // Running visualization code from command line
                    int returnVal = system(command.c_str());
                // Checking status of visualization system call
                if (returnVal != 0){
                    cout << "Error Encountered during Visualization\n";
                }	
            }
        }

    // DBSCAN

    } else {

        // Check if required options are provided
        if (infilename == NULL || minPts < 0 || eps < 0 || threads < 1) {
            usage(argv[0]);
            exit(-1);
        }

        cout << endl;
        cout << "   || DBSCAN ||    " << endl;
        cout << endl;

        // Set the number of OpenMP threads
        omp_set_num_threads(threads);

        // Create an instance of the ClusteringAlgo class
        NWUClustering::ClusteringAlgoDBS dbs;
        dbs.set_dbs_params(eps, minPts);
		dbs.csvOutputFilename = csvOutputFilename;
        cout << "Input parameters " << " minPts " << minPts << " eps " << eps << endl;

        // Measure the time to Read File Data
        cout << "Reading points from file: " << infilename << endl;
        double start = omp_get_wtime();

        // Read input data from the file
        if (dbs.read_file(infilename, isBinaryFile) == -1) {
            cout << "File Reading Failed." << endl;
            exit(-1);
        }

        cout << "Reading input data file took: " << omp_get_wtime() - start << " seconds." << endl;

        // Build a k-d tree for the points
        cout << "Begin Building kdtree." << endl;
        start = omp_get_wtime();
        dbs.build_kdtree();
        cout << "Building kdtree took: " << omp_get_wtime() - start << " seconds." << endl;

        // Run the SNG clustering algorithms
        start = omp_get_wtime();

        //Selecting Algorithm 
        cout << endl; 
        if (classical == true) {
            cout << "- DBSCAN SEQUENTIAL ALGORITHM -" << endl;
            cout << endl;
            run_dbs_algo(dbs); // Sequential 
        } else {
            cout << "- DBSCAN Parallel ALGORITHM -" << endl;
            cout << endl;
            run_dbs_algo_uf(dbs); // Parallel
        }
            
        cout << "DBSCAN (total) took " << omp_get_wtime() - start << " seconds." << endl;

        // If an output filename is provided, write clustering results to the file
        if (outfilename != NULL) {
            ofstream outfile;
            outfile.open(outfilename);
                if (classical == true) {
                dbs.writeClustersDBS(outfile); //Sequential 
            } else {
                dbs.writeClusters_ufDBS(outfile); // Parallel
            }

            outfile.close();
        }
    }

    return 0;
}
