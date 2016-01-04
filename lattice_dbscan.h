/*

TODO: Short Intro. Algorithm implementation inspired by Jono Brogan: 
http://codereview.stackexchange.com/questions/23966/density-based-clustering-of-image-keypoints

*/

#include <iostream>
#include <memory>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TTree.h"


#define DIM 256
#define NUM_TYPE short int

// #define MATRIX_LOOP(OUTER_STATEMENT, INNER_STATEMENT) \
// for (int i=0; i < DIM; i++) \
// { \
// 	OUTER_STATEMENT	\
// 	for(int j = 0; j<DIM; j++) \
// 	{ \
// 		INNER_STATEMENT \
// 	} \
// } \

typedef std::unique_ptr<TFile> TFilePtr;
typedef std::unique_ptr<TTree> TTreePtr;

using namespace std;

typedef vector<NUM_TYPE> data_point;
typedef vector<data_point> point_set;

class LatticeDBScan
{
public:
	LatticeDBScan();
	LatticeDBScan(std::string input_filename);
	//TODO:
	LatticeDBScan(const LatticeDBScan &src_LatticeDBScan);
	~LatticeDBScan();

	void FlushClustersToFile(const char * output_tfile, const char * title);
	void ClusterProjectedPoints(double epsilon, int min_pts);

	void CreateProjectedPoints();
	point_set GetProjectedPoints();
	vector<point_set> GetClusters();
	NUM_TYPE ** GetOriginalDataMatrix();

private:
	point_set R2_projected_points;
	vector<point_set> clusters;
	//Need to do dynamic allocation for this matrix
	NUM_TYPE ** original_data_matrix; // matrix DIM x DIM
	double GetEuclideanDistance(const data_point p1, const data_point p2);
	double ComputeAverageEnergy();
	vector<int> RegionQuery(vector<NUM_TYPE> point, double epsilon);
};

//TODO: PUT EVERYTHING IN A UNIQUE POINTER

LatticeDBScan::LatticeDBScan()
{
	//Initialize the data matrix
	original_data_matrix = new NUM_TYPE * [DIM];

	for(int i=0; i<DIM; i++)
	{
		original_data_matrix[i] = new NUM_TYPE[DIM];
	}
}

LatticeDBScan::~LatticeDBScan()
{
	for (int i =0; i < DIM; i++)
	{
		delete [] original_data_matrix[i];
	}

	delete [] original_data_matrix;
};

LatticeDBScan::LatticeDBScan(const LatticeDBScan &src_LatticeDBScan)
{
	//Initialize the data matrix
	original_data_matrix = new NUM_TYPE * [DIM];

	for(int i=0; i<DIM; i++)
	{
		original_data_matrix[i] = new NUM_TYPE[DIM];
	}

};

void LatticeDBScan::FlushClustersToFile(const char * output_tfile, const char * title)
{
	TFilePtr output_file(new TFile(output_tfile, "RECREATE"));
	TTreePtr tree(new TTree(title, ""));

	// ID of the current cluster
	Int_t clusterID;

	/* Calculate the maximum size of one of our clusters so we can initialize
	 * arrays to hold the data points
	 */
	Int_t max_size_of_cluster=0;
	for (int i =0; i < clusters.size(); i++)
	{
		int current_size = clusters[i].size();
		if (current_size > max_size_of_cluster)
			max_size_of_cluster = current_size;
	}

	// X-coordinate of cluster pixels
	NUM_TYPE x_coordinates[max_size_of_cluster];
	// Y-coordinates of cluster pixels
	NUM_TYPE y_coordinates[max_size_of_cluster];
	// "Density" - the amount of energy (data points) in a pixel
	Int_t pix_density[max_size_of_cluster];

	tree->Branch("MaxClusterSize", &max_size_of_cluster, "max_size_of_cluster/I");
	tree->Branch("ClusterID", &clusterID, "clusterID/I");
	tree->Branch("x_pixels", x_coordinates, "x_coordinates[max_size_of_cluster]/S");
	tree->Branch("y_pixels", y_coordinates, "y_coordinates[max_size_of_cluster]/S");
	tree->Branch("Pixel_density", pix_density, "pix_density[max_size_of_cluster]/I");

	for (int i =0; i < clusters.size(); i++)
	{	
		/* First clear the arrays from the contents of the last cluster
		 * (where -1 implies that it's an empty slot)
		 */
		for (int j=0; j < max_size_of_cluster; j++)
		{
			x_coordinates[j] = -1;
		}
		for (int j=0; j < max_size_of_cluster; j++)
		{
			y_coordinates[j] = -1;
		}
		/* Initialize the ensity to 0. */
		for (int j=0; j < max_size_of_cluster; j++)
		{
			pix_density[j] = 0;
		}

		// Add the current cluster
		clusterID = i;

		int cluster_size = 1;
		for (int j=0; j < clusters[i].size(); j++)
		{
			NUM_TYPE current_x_coordinate = clusters[i][j][0];
			NUM_TYPE current_y_coordinate = clusters[i][j][1];
			bool exists = false;

			// Check if element already exists, and if so bump the density
			for (int k=1; k < cluster_size; k++)
			{
				if ( (x_coordinates[k] == current_x_coordinate) &&
					(y_coordinates[k] == current_y_coordinate) )
					{
						exists = true;
						pix_density[k]++;
						break;
					} 
			}

			// If the element doesn't exist, add it to the list of pixels in the cluster
			if (!exists)
			{
				x_coordinates[cluster_size] = current_x_coordinate;
				y_coordinates[cluster_size] = current_y_coordinate;
				pix_density[cluster_size] = 0;
				cluster_size++;
			}
		}

		// Append the size of the cluster as the first element, for easy access
		x_coordinates[0] = cluster_size;
		y_coordinates[0] = cluster_size;
		pix_density[0] = cluster_size;

		// Flush the cluster to file
		tree->Fill();
		tree->Write();
	}
};

double
LatticeDBScan::GetEuclideanDistance(const data_point p1, const data_point p2)
{
	//FIX THIS
	vector<NUM_TYPE>::const_iterator itr;
	vector<NUM_TYPE>::const_iterator itr2 = p2.begin();
	double distance = 0;

	for (itr = p1.begin(); itr!=p1.end(); ++itr)
	{
		distance += pow((*itr) - (*itr2),2);
		++itr2;
	}

	return sqrt(distance);
}

LatticeDBScan::LatticeDBScan(std::string input_filename)
{
	//Initialize the data matrix
	original_data_matrix = new NUM_TYPE * [DIM];

	for(int i=0; i<DIM; i++)
	{
		original_data_matrix[i] = new NUM_TYPE[DIM];
	}

	std::ifstream input_file;
	input_file.open(input_filename.c_str(), ios::in);

	std::string value;
	std::string line;

	for (int i=0; i < DIM; i++)
	{
		getline(input_file, line);
		istringstream string_stream(line);

		for(int j = 0; j<DIM; j++) \
		{
			string_stream >> value;
			original_data_matrix[i][j] = std::stoi(value);
		}
	}
}

double LatticeDBScan::ComputeAverageEnergy()
{
	double average =0; 

	for (int i=0; i < DIM; i++)
	{
		for(int j = 0; j<DIM; j++) \
		{
			average += original_data_matrix[i][j];
		}
	}

	average /= (DIM * DIM);
	return average;
}

void LatticeDBScan::CreateProjectedPoints()
{
	double energy_threshold = ComputeAverageEnergy();

	for (int i=0; i < DIM; i++)
	{
		for(int j = 0; j<DIM; j++) \
		{
			if (original_data_matrix[i][j] != 0)
			{
				if (original_data_matrix[i][j] < energy_threshold)
				{
					data_point vec; 
					vec.push_back(i);
					vec.push_back(j);
					R2_projected_points.push_back(vec);
				}
				else
				{
					data_point vec; 
					vec.push_back(i);
					vec.push_back(j);
					//3 times, since high energy
					R2_projected_points.push_back(vec);
					R2_projected_points.push_back(vec);
					R2_projected_points.push_back(vec);
				}
			}
		}
	}

	// MATRIX_LOOP(,
	// 	if (original_data_matrix[i,j] != 0)
	// 	{
	// 		if (original_data_matrix[i,j] < energy_threshold)
	// 		{
	// 			data_point vec; 
	// 			vec.push_back(i);
	// 			vec.push_back(j);
	// 			point_set.push_back(vec);
	// 		}
	// 		else
	// 		{
	// 			data_point vec; 
	// 			vec.push_back(i);
	// 			vec.push_back(j);
	// 			//3 times, since high energy
	// 			point_set.push_back(vec);
	// 			point_set.push_back(vec);
	// 			point_set.push_back(vec);
	// 		}
	// 	}
		// )
}

void LatticeDBScan::ClusterProjectedPoints(double epsilon, int min_pts)
{

	vector<bool> clustered;
	vector<bool> visited;
	vector<int> noise;
	vector<int> neighborPts;
	vector<int> neighborPts_;
	int c = 0;

	int total_points = R2_projected_points.size();

	//init clustered and visited
	for(int k = 0; k < total_points; k++)
	{
	    clustered.push_back(false);
	    visited.push_back(false);
	}

	int cluster_size = 0;
	// clusters.push_back(point_set()); //will stay empty?

	//for each unvisted point P in dataset keypoints
	for(int i = 0; i < total_points; i++)
	{
	    if(!visited[i])
	    {
	        visited[i] = true;
	        neighborPts = RegionQuery(R2_projected_points.at(i),epsilon);

	        if(neighborPts.size() < min_pts)
	            //Mark as noise
	            noise.push_back(i);

	        else
	        {
	            clusters.push_back(point_set());

	            //Start expanding cluster
	            cluster_size++;
	            // add original core point to cluster 
	            clusters[c].push_back(R2_projected_points.at(i));

	            //Add each point in neighborPts
	            for(int j = 0; j < neighborPts.size(); j++)
	            {
	                //if P' is not visited
	                if(!visited[neighborPts[j]])
	                {
	                    //Mark P' as visited
	                    visited[neighborPts[j]] = true;
	                    neighborPts_ = RegionQuery(R2_projected_points.at(neighborPts[j]),
	                    	epsilon);

	                    //Add more points to process
	                    if(neighborPts_.size() >= min_pts)
	                    {
	                        neighborPts.insert(neighborPts.end(),neighborPts_.begin(),neighborPts_.end());
	                    }
	                }
	                // if P' is not yet a member of any cluster
	                // (I.E. IF IT IS VISITED, BUT HAS NOT YET BEEN CLUSTERED)
	                // add P' to cluster c
	                if(!clustered[neighborPts[j]])
	                {
	                    clusters[c].push_back(R2_projected_points.at(neighborPts[j]));
	                    cluster_size++;
	                }
	            }

	            //We added one cluster
	            c++;
	        }

	    }
	}
}

vector<int> LatticeDBScan::RegionQuery(vector<NUM_TYPE> point, double epsilon)
{
	double dist;
	vector<int> close_points;

	for(int i = 0; i< R2_projected_points.size(); i++)
	{
	    dist = GetEuclideanDistance(point, R2_projected_points[i]);
	  
	    if(dist <= epsilon && dist != 0.0f)
	    {
	        close_points.push_back(i);
	    }
	}
	return close_points;
}

point_set LatticeDBScan::GetProjectedPoints()
{
	return R2_projected_points;
}

vector<point_set> LatticeDBScan::GetClusters()
{
	return clusters;
}

NUM_TYPE ** LatticeDBScan::GetOriginalDataMatrix()
{
	return original_data_matrix;
}
