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

void displayScatterPlot2D(NWUClustering::Clusters * algorithm)
{
	//Fill x and y vectors
	
	
	vector<double> x;
	vector<double> y;
	vector<double> c;
	//Iterate through each point
	for(int i = 0; i < algorithm->m_pts->m_i_num_points; i++)
	{
		x.push_back(algorithm->m_pts->m_points[i][1]);
		y.push_back(algorithm->m_pts->m_points[i][2]);
		c.push_back(algorithm->m_pts->m_points[i][0] * (1000));
		//c.push_back(i * 20);
	}
	
	matplot::scatter(x,y,vector<double>{}, c);
	matplot::show();
}

void displayScatterPlot3D(NWUClustering::Clusters * alogrithm)
{
	
}
