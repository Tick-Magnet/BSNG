#pragma once

#include <matplot/matplot.h>
#include <vector>
#include "../sng.h"

void displayScatterPlot2D(NWUClustering::ClusteringAlgo * algorithm, vector<int> * clusters);
void displayScatterPlot2D(NWUClustering::ClusteringAlgo * algorithm);


void displayScatterPlot3D(NWUClustering::ClusteringAlgo * alogrithm);
