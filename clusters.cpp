//Clusters.cpp prepares the data for clustering, by reading the file and creating the kdtree.

#include "clusters.h"

namespace NWUClustering {

	// Destructor for the Clusters class
	Clusters::~Clusters() {
       
	    // Clean up allocated memory for points and k-d tree
		if(m_pts){
			m_pts -> m_points.clear();
			delete m_pts;
			m_pts = nullptr;
		}

		if(m_kdtree)
		{
			delete m_kdtree;
			m_kdtree = nullptr;
		}
	}

	// Function to read data from a file
	int Clusters::read_file(char* infilename, int isBinaryFile) {
		
		ssize_t numBytesRead;
		int i, j;
		int num_points, dims;

		if(isBinaryFile == 1) {
			
			///		Binary Text Files
			cout << " - - -  Binary Text File Inputted - - - " << endl ;

			ifstream file (infilename, ios::in|ios::binary);
			if(file.is_open())
  			{
				// Read the number of points and dimensions
				file.read((char*)&num_points, sizeof(int));
				file.read((char*)&dims, sizeof(int));
				cout << "Points " << num_points << " dims " << dims << endl;

				//Allocate memory for points
				m_pts = new Points;				
				m_pts->m_i_dims = dims;
                m_pts->m_i_num_points = num_points;
				
                m_pts->m_points.resize(num_points);
                for(int ll = 0; ll < num_points; ll++)
                    m_pts->m_points[ll].resize(dims);

				point_coord_type* pt;					
				
                pt = new point_coord_type[dims];
                        
                for (i = 0; i < num_points; i++){
                    file.read((char*)pt, dims*sizeof(point_coord_type));
                                
                    for (j = 0; j < dims; j++)
                        m_pts->m_points[i][j] = pt[j];
                }
			
				delete [] pt;	
				file.close();

  			} else {
				
				cout << "Error: no such file: " << infilename << endl;
				return -1;
			}

		} else {
			
			///		Non-Binary Text Files
			cout << " - - -  Non-Binary Text File Inputted - - - " << endl;

			// Process text file
			string line, line2, buf;
			ifstream file(infilename);
			stringstream ss;

			if (file.is_open())
  			{
				// Get the dimensions from the first line
				getline(file, line);
				line2 = line;
				ss.clear();				
				ss << line2;
			
				dims = 0;
				while(ss >> buf) // get the corordinate of the points
					dims++;

				// Count the number of points
				num_points = 0;
                while (!file.eof()) {
                    if (line.length() == 0)
                        continue;
                    num_points++;
                    getline(file, line);
                }
				
				cout << "Points " << num_points << " dimensions " << dims << endl;
                               
				// Allocate memory for points
				m_pts = new Points;
				m_pts->m_points.resize(num_points);
                for(int ll = 0; ll < num_points; ll++)
                    m_pts->m_points[ll].resize(dims);
				
				file.clear();
				file.seekg (0, ios::beg);
				
				getline(file, line);

				i = 0;
    			while (!file.eof()){
					if(line.length() == 0)
						continue;

					ss.clear();
					ss << line;

					j = 0;
					while(ss >> buf && j < dims) //Get the corordinate of the points
					{
						m_pts->m_points[i][j] = atof(buf.c_str());
						j++;
					}
					
					i++;
					getline(file, line);
    			}

    			file.close();
  			
                m_pts->m_i_dims = dims;
            	m_pts->m_i_num_points = num_points;

			} else {
                cout << "Error: no such file: " << infilename << endl;
                return -1;
			}			
		}
		
		return 0;		
	}

	// Function to build a k-d tree
    int Clusters::build_kdtree() {
        if (m_pts == nullptr) {
            cout << "Point set is empty" << endl;
            return -1;
        }

        m_kdtree = new kdtree2(m_pts->m_points, false); //kdtree constructor

        if (m_kdtree == nullptr) {
            cout << "Failed to allocate a new kd tree" << endl;
            return -1;
        }
		
		cout << "kdtree successfully built" << endl ;
        return 0;
    }
}