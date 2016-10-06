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
#include <algorithm>
#include "TTree.h"
#include "TFile.h"

#include "TPaveStats.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphPainter.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TPad.h"
#include "TError.h"
#include "TH1I.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TVirtualFitter.h"
#include "TFormula.h"
#include "TList.h"
#include "TLatex.h"
#include "TH2I.h"
#include "TSystem.h"
#include "TF1.h"
#include "Riostream.h"


#define DIM 256
#define NUM_TYPE short int

TCanvas * c0;

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

	void FlushClustersToFile(const char * output_tfile, const char * title, bool root_file = true);
	void ClusterProjectedPoints(double epsilon, int min_pts);
	void PlotClusters(std::string cuts="", std::string _title = "", std::string outname = "");

	void CreateProjectedPoints();
	point_set GetProjectedPoints();
	vector<point_set> GetClusters();
	NUM_TYPE ** GetOriginalDataMatrix();

private:
	point_set R2_projected_points;
	vector<point_set> clusters;
	NUM_TYPE ** original_data_matrix; // matrix DIM x DIM
	double GetEuclideanDistance(const data_point p1, const data_point p2);
	NUM_TYPE ComputePercentileEnergy(const double percentile);
	vector<int> RegionQuery(vector<NUM_TYPE> point, double epsilon);

	void set_plot_style();
	int put_clusters_in_frame(std::vector<int> &Frame_Density, std::vector<int> &cluster_x_coords,
	std::vector<int> &cluster_y_coords, int cluster_index);
	void defineHistogram(TH1* h, string x, string y, string title, bool adjust_range = false, double titleSize = 0.075, bool rebin = false, double average_entries_per_bin = 10.);
};

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
}

// TODO:
LatticeDBScan::LatticeDBScan(const LatticeDBScan &src_LatticeDBScan)
{
	//Initialize the data matrix
	original_data_matrix = new NUM_TYPE * [DIM];

	for(int i=0; i<DIM; i++)
	{
		original_data_matrix[i] = new NUM_TYPE[DIM];
	}

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

		for(int j = 0; j<DIM; j++)
		{
			string_stream >> value;
			original_data_matrix[i][j] = std::stoi(value);
		}
	}
}


void LatticeDBScan::FlushClustersToFile(const char * output_tfile, const char * title, bool root_file)
{
	if (!root_file)
	{
		for (int i=0; i < clusters.size(); i++)
		{

		}
	}
	
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
		/* Initialize the density to 0. */
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
				pix_density[cluster_size] = 1;
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
}

void 
LatticeDBScan::PlotClusters(std::string cuts, std::string _title, std::string outname)
{
	std::cout << " Starting to plot.. " << std::endl;

	cuts = (cuts).c_str();

	set_plot_style();

	TH2F* pix = new TH2F("pixelmatrix", "", 256, -0.5,255.5,256,-0.5,255.5);

	std::cout << "Got pixel matrix" << std::endl;

	pix->SetTitle("");
	pix->GetXaxis()->SetTitle("Pixel No. X");
	pix->GetYaxis()->SetTitle("Pixel No. Y");

	std::cout << "Set pixel matrix details" << std::endl;

	for(int i = 0; i < clusters.size(); i++)
	{
		// Holder of frames of current cluster to plot
		std::vector<int> frame_density;
		// Holder of coordinates of current cluster to plot
		std::vector<int> x_coordinates;
		std::vector<int> y_coordinates;
		// Fill the vectors with appropriate data
		int cluster_size = put_clusters_in_frame(frame_density, x_coordinates,
		y_coordinates, i);

		// Then fill the display frame with the cluster data
		for(int j = 0; j < cluster_size; j++)
		{
			pix->Fill(y_coordinates[j],x_coordinates[j],frame_density[j]);
			// std::cout << x_coordinates[j] << " " << y_coordinates[j] << " " << frame_density[j] << std::endl;
		}
	}		

	std::cout << "Finished filling out the frame " << std::endl;

	pix->SetMinimum(0);
	//pix->SetMaximum(78);
	std::cout << "Error right before getting drawable object" << std::endl;
	c0 = new TCanvas("New Canvas", "", 665, 662);
	std::cout << "Error right after getting drawable object" << std::endl;
	//c0->SetRightMargin(0.23);
	c0->SetLeftMargin(0.23);
	c0->SetRightMargin(0.12);
	c0->SetBottomMargin(0.23);
	

	TPaveText* pt = new TPaveText(0.75, 0.91 ,0.97, 0.99, "brNDC");

	std::cout << "Added text" << std::endl;

	pix->SetMinimum(0.9);
	pt->AddText("integral"); pt->AddText("clock counts");


	
	defineHistogram(pix, "Pixel No. X", "Pixel No. Y", "");
	//c0->SetGrid();
	pt->SetFillColor(kWhite);
	pt->SetFillStyle(0);
	pt->SetBorderSize(0);
	pt->SetShadowColor(0);
	pt->SetTextFont(42);
	pt->SetTextSize(0.05);

	std::cout << " Midway point" << std::endl;

	//pt->SetTextAttributes(42);
	c0->SetFixedAspectRatio(1);
	c0->SetLogz();
	pix->GetYaxis()->SetTitleOffset(1.5);
	pix->SetStats(kFALSE);
	pix->SetTitle(_title.c_str());
	pix->Draw("COLZ9");
	pt->Draw();

}

int
LatticeDBScan::put_clusters_in_frame(std::vector<int> &Frame_Density, std::vector<int> &cluster_x_coords, 
	std::vector<int> &cluster_y_coords, int cluster_index)
{	
	int cluster_size = clusters[cluster_index].size();

	int current_size = 0;
	for (int j=0; j < cluster_size; j++)
	{
		NUM_TYPE current_x_coordinate = clusters[cluster_index][j][0];
		NUM_TYPE current_y_coordinate = clusters[cluster_index][j][1];
		bool exists = false;

		// Check if element already exists, and if so bump the density
		for (int k=0; k < current_size; k++)
		{
			if ( (cluster_x_coords[k] == current_x_coordinate) &&
				(cluster_y_coords[k] == current_y_coordinate) )
			{
 				exists = true;
				Frame_Density[k]++;
				break;
			} 
		}

		// If the element doesn't exist, add it to the list of pixels in the cluster
		if (!exists)
		{
			cluster_x_coords.push_back(current_x_coordinate);
			cluster_y_coords.push_back(current_y_coordinate);
			Frame_Density.push_back(1);
			current_size++;
		}
	}
	// std::cout << "Finished putting clusters in frame!" << std::endl;
	// Return the size of the vectors of frames and coordinates
	return current_size;
}

NUM_TYPE 
LatticeDBScan::ComputePercentileEnergy(const double percentile)
{
	// unsigned int percentile_size;

	// std::vector<NUM_TYPE> lower_heap;
	// std::vector<NUM_TYPE> upper_heap;
	// std::make_heap(lower_heap.begin(), lower_heap.end());
	// std::make_heap(upper_heap.begin(), upper_heap.end());

	// for (int i=0; i < DIM; i++)
	// {
	// 	for(int j = 0; j<DIM; j++)
	// 	{
	// 		/* If we have an empty heap always push it to that heap so
	// 		   that the front() operation is defined */
	// 		if (lower_heap.empty())
	// 			lower_heap.push_back(original_data_matrix[i][j]);
	// 		else if (upper_heap.empty())
	// 			upper_heap.push_back(original_data_matrix[i][j]);
	// 		else
	// 		{
	// 			if (original_data_matrix[i][j] <= lower_heap.front())
	// 			{
	// 				lower_heap.push_back(original_data_matrix[i][j]);
	// 				std::push_heap(lower_heap.begin(), lower_heap.end(), std::less<NUM_TYPE>());
	// 			}
	// 			else
	// 			{
	// 				upper_heap.push_back(original_data_matrix[i][j]);
	// 				std::push_heap(upper_heap.begin(), upper_heap.end(), std::greater<NUM_TYPE>());
	// 			}
	// 		}

	// 		// Reorganize so that percentiles are maintained
	// 		percentile_size = ((lower_heap.size() + upper_heap.size()) * percentile) + 1;
	// 		// std::cout << "lh: " << lower_heap.size() << std::endl;
	// 		// std::cout << "uh: " << upper_heap.size() << std::endl;
	// 		// std::cout << percentile_size << std::endl;

	// 		if (lower_heap.size() > percentile_size)
	// 		{
	// 			/*  First move the maximal element to the end of the heap, and reorganize
	// 			  the rest of the heap, O(logn) heapify operation */
	// 			std::pop_heap(lower_heap.begin(), lower_heap.end());
	// 			/* Push the maximal element from the lower heap to the upper heap,
	// 			 O(1) time */
	// 			upper_heap.push_back(lower_heap.back());
	// 			/* Reorganize the upper heap so that the heap properties are maintained
	// 			 with the new added element, O(logn) heapify ooperation*/
	// 			std::push_heap(upper_heap.begin(), upper_heap.end(), std::greater<NUM_TYPE>());
	// 			/* Now remove the maximal element from the lower heap, it is now located as the
	// 			 last element of the array, therefore O(1) time*/
	// 			lower_heap.pop_back();
	// 		}

	// 		if (lower_heap.size() < percentile_size)
	// 		{
	// 			/*  First move the maximal element to the end of the heap, and reorganize
	// 			  the rest of the heap, O(logn) heapify operation */
	// 			std::pop_heap(upper_heap.begin(), upper_heap.end(), std::greater<NUM_TYPE>());
	// 			/* Push the maximal element from the lower heap to the upper heap,
	// 			 O(1) time */
	// 			lower_heap.push_back(upper_heap.back());
	// 			/* Reorganize the upper heap so that the heap properties are maintained
	// 			 with the new added element, O(logn) heapify ooperation*/
	// 			std::push_heap(lower_heap.begin(), lower_heap.end(), std::less<NUM_TYPE>());
	// 			/* Now remove the maximal element from the lower heap, it is now located as the
	// 			 last element of the array, therefore O(1) time*/
	// 			upper_heap.pop_back();
	// 			percentile_size = ((lower_heap.size() + upper_heap.size()) * percentile) + 1;
	// 		}
	// 	}
	// }

	// Finally, take the maximal element from the lower heap as the desired percentile
	// std::cout << "Lower heap" << lower_heap.back() << std::endl;
	// return lower_heap.front();

	std::vector<NUM_TYPE> pix_values;
	for (int i=0; i < DIM; i++) 
	{
		for (int j=0; j < DIM; j++) 
		{
			pix_values.push_back(original_data_matrix[i][j]);
		}
	}

	NUM_TYPE max_el = *(std::max_element(pix_values.begin(), pix_values.end()));
	NUM_TYPE min_el = *(std::min_element(pix_values.begin(), pix_values.end()));

	return ( min_el + ( (max_el - min_el) * percentile) );
}


void LatticeDBScan::CreateProjectedPoints()
{
	// double energy_threshold = ComputePercentileEnergy(0);
	NUM_TYPE energy_threshold = ComputePercentileEnergy(.45);

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

void 
LatticeDBScan::
defineHistogram(TH1* h, string x, string y, string title, bool adjust_range, double titleSize, bool rebin, double average_entries_per_bin)
{
	h->GetXaxis()->SetTitle(x.c_str());
	h->GetYaxis()->SetTitle(y.c_str());
	h->SetStats(kFALSE);
	h->SetLineWidth(2);
	h->SetTitle(title.c_str());
	h->GetYaxis()->SetDecimals();
	h->GetYaxis()->SetTitleOffset(1.2);
	h->GetXaxis()->SetTitleSize(titleSize);
	h->GetYaxis()->SetTitleSize(titleSize);
	h->GetXaxis()->SetLabelSize(titleSize);
	h->GetYaxis()->SetLabelSize(titleSize);
	h->SetTitleSize(titleSize);
	
	if(x.compare("#varphi (Deg)")==0)
		h->GetXaxis()->SetNdivisions(505, true);
	else
		h->GetXaxis()->SetNdivisions(503, true);
	
	h->GetYaxis()->SetNdivisions(505, true);

	int lastBin = h->GetNbinsX();
	if(adjust_range == true) {
		for(int i = 1; i <= h->GetNbinsX(); i++) {
			if(h->GetBinContent(i) > 0)
				lastBin = i;
		}
		lastBin= (int) lastBin * 1.2;
		//if(lastBin > 3 * h->GetMean())
		//	lastBin = 3 * h->GetMean();

		if(rebin == true) {
			int bin_factor = (average_entries_per_bin * lastBin / h->GetEntries());
			if(bin_factor > 1) {
				h->Rebin(bin_factor);
				lastBin = lastBin / bin_factor * 1.2;
			}
		}

		h->GetXaxis()->SetRange(0,lastBin);
	}
 }


void 
LatticeDBScan::set_plot_style()
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    gStyle->SetNumberContours(NCont);
}
