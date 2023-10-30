#include "sng.h"
#include "utils.h"
#include "kdtree2.hpp"

// Function to display usage information
static void usage(char* argv0) {
    const char* params =
        "Usage: %s [switches] -i filename -b -m minpts -e epsilon -o output -t threads\n"
	"Example: ./bsng -i clus50k.bin -b -m 5 -e 25 -t 8 -o output\n"
        "    -i filename     : file containing input data to be clustered\n"
        "    -b isBinary     : input file is in binary format (default no)\n"
        "    -m minpts       : input parameter of BSNG, min points to form a cluster, e.g. 2\n"
        "    -e epsilon      : input parameter of BSNG, radius or threshold on neighborhoods retrieved, e.g. 0.8\n"
        "    -s seeds        : input parameter of Sow-n-Grow, number of seeds to be chosen from each partition of data\n"
        "    -o output       : clustering results, format, (each line, point id, clusterid)\n"
        "    -t threads      : number of threads to be employed\n\n";

    fprintf(stderr, params, argv0);
    exit(-1);
}


int main(int argc, char** argv) {
    double seconds, eps;
    int minPts, threads, opt, seeds;
    char* outfilename = NULL;
    char* infilename = NULL;
    bool classical, isBinaryFile;

    // Initialize default values
    minPts = -1;
    eps = -1;
    isBinaryFile = 0;
    outfilename = NULL;
    infilename = NULL;
    threads = -1;
    seeds = -1;
    classical = false;
    isBinaryFile = false;

    while ((opt = getopt(argc, argv, "i:t:d:p:m:e:s:o:v:z:bxghncul")) != EOF) {
        switch (opt) {
            case 'i':
                infilename = optarg;
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
            case 'o':
                outfilename = optarg;
                break;
            case 'c':
                classical = true;
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
    }

	// Check if required options are provided
    if (infilename == NULL || minPts < 0 || seeds < 0 || eps < 0 || threads < 1) {
        usage(argv[0]);
        exit(-1);
    }

	// Set the number of OpenMP threads
    omp_set_num_threads(threads);

    // Create an instance of the ClusteringAlgo class
    NWUClustering::ClusteringAlgo sng;
    sng.set_sng_params(eps, minPts, seeds);

    cout << "Input parameters " << " minPts " << minPts << " eps " << eps << " seeds " << seeds << endl;

    // Measure the start time
    double start = omp_get_wtime();
    cout << "Reading points from file: " << infilename << endl;

    // Read input data from the file
    if (sng.read_file(infilename, isBinaryFile) == -1)
        exit(-1);

    cout << "Reading input data file took " << omp_get_wtime() - start << " seconds." << endl;

    // Build a k-d tree for the points
    start = omp_get_wtime();
    sng.build_kdtree();
    cout << "Build kdtree took " << omp_get_wtime() - start << " seconds." << endl;

    // Run the SNG clustering algorithms
    start = omp_get_wtime();

    if (classical == true) {
        run_sng_algo(sng); // Sequential 
    } else {
        run_sng_algo_uf(sng); // Parallel
    }
        
    cout << "S&G (total) took " << omp_get_wtime() - start << " seconds." << endl;

	
	// If an output filename is provided, write clustering results to the file
    if (outfilename != NULL) {
        ofstream outfile;
        outfile.open(outfilename);

            if (classical == true) {
            sng.writeClusters(outfile); //Sequential 
        } else {
            sng.writeClusters_uf(outfile); // Parallel
        }

        outfile.close();
    }

    return 0;
}