#include "Visualization.h"


#include <random>
/*
// Write point id and cluster ids to the output stream
		for (i = 0; i < m_pts->m_i_num_points; i++) {
            //Formatting for Output
            if (i < 10){
                o << i << "  | " << clusters[m_parents[i]] << endl;
                o << " Value: (";
//---------------------DATA OUPUT CHANGE----------------------                
				//Outputing datapoint 
				for(int j = 0; j < m_pts->m_i_dims; j++)
				{
					o << m_pts->m_points[i][j] << ", ";
				}
				o << ")" << endl;
				
*/

void displayScatterPlot2D(NWUClustering::ClusteringAlgo * algorithm, vector<int> * clusters)
{
	//Fill x and y vectors
	
	
	vector<double> x;
	vector<double> y;
	vector<double> c;
	//Iterate through each point
	for(int i = 0; i < algorithm->m_pts->m_i_num_points; i++)
	{
		x.push_back(algorithm->m_pts->m_points[i][2]);
		y.push_back(algorithm->m_pts->m_points[i][3]);
		c.push_back(clusters->at(algorithm->m_parents[i]) * 1000);
		//ls
		cout << "X: " << algorithm->m_pts->m_points[i][2] << "Y: " << algorithm->m_pts->m_points[i][3] << endl;
		cout << clusters->at(algorithm->m_parents[i]) << endl;
		//c.push_back(0);
		//c.push_back(i * 20);
	}
	
	auto plot = matplot::scatter(x,y,vector<double>{}, c);
	plot->marker_face(true);
	
	matplot::show();
}
void displayScatterPlot2D(NWUClustering::ClusteringAlgo * algorithm)
{
	//Fill x and y vectors
	
	
	vector<double> x;
	vector<double> y;
	vector<double> c;
	vector<int> clusters = algorithm->m_pid_to_cid;
	//Iterate through each point
	for(int i = 0; i < algorithm->m_pts->m_i_num_points; i++)
	{
		x.push_back(algorithm->m_pts->m_points[i][2]);
		y.push_back(algorithm->m_pts->m_points[i][3]);
		c.push_back(clusters.at(i) * 1000);
		//ls
		//cout << "X: " << algorithm->m_pts->m_points[i][2] << "Y: " << algorithm->m_pts->m_points[i][3] << endl;
		//cout << clusters.at(algorithm->m_parents[i]) << endl;
		//c.push_back(0);
		//c.push_back(i * 20);
	}
	
	auto plot = matplot::scatter(x,y,vector<double>{}, c);
	plot->marker_face(true);
	
	matplot::show();
}

void displayScatterPlot3D(NWUClustering::ClusteringAlgo * alogrithm)
{
	
}
