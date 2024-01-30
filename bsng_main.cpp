#include "sng.h"
#include "dbscan.h"
#include "utils.h"
#include "kdtree2.hpp"

// Function to display usage information
static void usage(char* argv0) {
    const char* params =
    "\n"
    "Usage: %s [switches] -i filename -b -m minpts -e epsilon -t threads -s seeds -o output \n"
	"Example: ./bsng -i random_points.bin -b -m 3 -e 5 -s 2 -t 2 -o output\n"
        "\n"
        "    -i filename     : file containing input data to be clustered\n"
        "    -d dbscan       : run dbscan instead of sow-and-grow (default no)\n"
        "    -b isBinary     : input file is in binary format (default no)\n"
        "    -m minpts       : input parameter of BSNG, min points to form a cluster, e.g. 3\n"
        "    -e epsilon      : input parameter of BSNG, radius of neighborhood, e.g. 5\n"
        "    -s seeds        : input parameter of BSNG, input seed(s) (default 1)\n"
        "    -t threads      : number of threads to be employed (default 1)\n"      
        "    -o output       : clustering results output file - format (point id, clusterid)\n"
        "    -z seedMethod   : select method for seed selection (default 0)\n"
        "       * 0 random   : random seed values\n"
        "       * 1 even     : even seed values\n"
        "       * 2 odd      : odd seed values\n"
        "       * 3 Lower    : all under input seed\n"
        "       * 4 Upper    : all over input seed\n"
        "       * 5 UserList :                                "
        "    -u userList     : Select the seeds to be used from a list"
  
    "\n";

    fprintf(stderr, params, argv0);
    exit(-1);
}


int main(int argc, char** argv) {
    double seconds, eps;
    int minPts, threads, opt, seeds, method;
    char* outfilename = NULL;
    char* infilename = NULL;
    char* csvOutputFilename = NULL;
    char* userfilename = NULL;
    bool classical, isBinaryFile, dbscan;

    // Initialize default values
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

    while ((opt = getopt(argc, argv, "i:u:t:p:m:e:s:o:v:s:z:bdxghncul")) != EOF) {
        switch (opt) {
			case 'v':
				csvOutputFilename = optarg;
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
            case '?':
                usage(argv[0]);
                break;
            default:
                usage(argv[0]);
                break;
        }

        //printf("opt: %c, optarg: %s\n", opt, optarg); //Debugging for Arguments
    }


    // SNG (Runs by Default)

    if (dbscan == false) 
    {

        // Check if required options are provided
        if (infilename == NULL || minPts < 0 || seeds < 0 || eps < 0 || threads < 0) {
            usage(argv[0]);
            exit(-1);
        }

        if (threads > seeds) {

            cout << "More Threads than Seeds." << endl;
            threads = seeds; 
        }
        
        cout << endl;
        cout << "   || SOW & GROW ||    " << endl;
        cout << endl;

        // Set the number of OpenMP threads
        omp_set_num_threads(threads);

        // Create an instance of the ClusteringAlgo class
        NWUClustering::ClusteringAlgo sng;
        sng.set_sng_params(eps, minPts, seeds, method);
		sng.csvOutputFilename = csvOutputFilename;
        cout << "Input parameters " << " minPts " << minPts << " eps " << eps << " seeds " << seeds << endl;

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
        cout << "Begin Building kdtree." << endl;
        start = omp_get_wtime();
        sng.build_kdtree();
        cout << "Building kdtree took: " << omp_get_wtime() - start << " seconds." << endl;

        // Run the SNG clustering algorithms
        start = omp_get_wtime();

        //Selecting Algorithm 
        cout << endl; 
        if (classical == true) {
            cout << "- SNG SEQUENTIAL ALGORITHM -" << endl;
            run_sng_algo(sng); // Sequential 
        } else {
            cout << "- SNG Parallel ALGORITHM -" << endl;
            run_sng_algo_uf(sng); // Parallel
        }
            
        cout << "S&G (total) took " << omp_get_wtime() - start << " seconds." << endl;

        // If an output filename is provided, write clustering results to the file
        if (outfilename != NULL) {
            ofstream outfile;
            outfile.open(outfilename);

                if (classical == true) {
                sng.writeClusters(outfile); // Sequential 
            } else {
                sng.writeClusters_uf(outfile); // Parallel
            }

            outfile.close();
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
