#include "lattice_dbscan.h"

#ifdef CREATEANDFLUSH
int
main(int argc, char * argv[])
{
	LatticeDBScan myscan("../output");
	myscan.CreateProjectedPoints();
	// point_set points = myscan.GetProjectedPoints();
	// for (int i =0; i < points.size(); i++)
	// 	std::cout << points[i][0] << "," << points[i][1] << std::endl;
	myscan.ClusterProjectedPoints(1,4);
	vector<point_set> clusters = myscan.GetClusters();
	std::cout << clusters.size() << std::endl;

	// for (int i =0; i < clusters.size(); i++)
	// {
	// 	point_set cur_cluster = clusters[i];
	// 	std::cout << cur_cluster.size() << std::endl;
	// 	std::cout << "Cluster " << i << ": ";
	// 	for (int j=0; j < cur_cluster.size(); j++)
	// 	{
	// 		std::cout << "(" << cur_cluster[j][0] << "," << cur_cluster[j][1] << ")" << " ";
	// 	}
	// 	std::cout << std::endl;
	// }

	myscan.FlushClustersToFile("cluster_file.root", "11_75degStableBeam.root, clusterType=4, cluster_num = 3");
}
#endif

#ifdef VIEW_CLUSTERS
int
main(int argc, char * argv[])
{
	TFilePtr file = TFilePtr(new TFile("cluster_file.root", "READ"));
	if (file == NULL)
	{
		std::cout << "Error opening file!" << std::endl;
		return 1;
	}

	TTree * tree = (TTree *)file->Get("11_75degStableBeam.root, clusterType=4, cluster_num = 3");
	if (tree == NULL)
	{
		std::cout << "Error getting tree from file" << std::endl;
		return 1;
	}

	// Get the max size of the cluster to initialize our holding arrays
	Int_t max_size_of_cluster;
	tree->GetBranch("MaxClusterSize")->SetAddress(&max_size_of_cluster);
	tree->GetEntry(0);

	Int_t clusterID;
	NUM_TYPE x_coordinates[max_size_of_cluster];
	NUM_TYPE y_coordinates[max_size_of_cluster];
	Int_t density[max_size_of_cluster];

	tree->GetBranch("ClusterID")->SetAddress(&clusterID);
	tree->GetBranch("x_pixels")->SetAddress(x_coordinates);
	tree->GetBranch("y_pixels")->SetAddress(y_coordinates);
	tree->GetBranch("Pixel_density")->SetAddress(density);

	int index = 0;
	int more_clusters=1;
	while ( (more_clusters = tree->GetEntry(index)) != 0 )
	{
		// Print cluster id and size
		std::cout << "Cluster ID: " << clusterID << " Size: " << y_coordinates[0] << std::endl;

		for (int k=1; k < x_coordinates[0]; k++)
		{
			// Print only heavy pixels (just for testing)
			if (density[k] == 0)
				continue;
			std::cout << "( " << x_coordinates[k] << "," << y_coordinates[k]
				<< " ): " << density[k] << std::endl; 
		}

		index++;
	}

	return 0;
}
#endif