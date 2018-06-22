#include <cstdlib>			// C standard includes
#include <iostream>			// C++ I/O
#include <string>			// C++ strings
#include <vector>
#include <map>
#include <list>

#include "k_means_wrapper.h"
#include "Kmeans/KMlocal.h"			// k-means algorithms

CKmeansWrapper::CKmeansWrapper()
{
}


CKmeansWrapper::~CKmeansWrapper()
{
}

//CKmeansWrapper::CKmeansWrapper(int k, int dim, int nPts)
//{
//	k_ = k;
//	dim_ = dim;
//	nPts_ = nPts;
//}

void CKmeansWrapper::setTermParameters()
{
	// under construction
	std::cout << "@CKmeansWrapper::setTermParameters: under construction" << std::endl;
}

void CKmeansWrapper::executeKMeans(
	int k, 
	std::vector<std::vector<double> > f_vec, 
	std::vector<int> &ClusterId,
	std::vector<std::vector<double> >& ClusCenters)
{
	std::vector<std::vector<int>> Labels2Clusters;
	std::vector<double> cluster_skewness;
	executeKMeans(k, f_vec, ClusterId, ClusCenters, Labels2Clusters);

	return;
}

void CKmeansWrapper::executeKMeans(
	int k, 
	std::vector<std::vector<double>> f_vec, 
	std::vector<int>& ClusterId, 
	std::vector<std::vector<double>>& ClusCenters, 
	std::vector<std::vector<int>>& Labels2Clusters)
{
	assert(k > 0);
	//std::cout << "#Clusters " << k << std::endl;
	assert(f_vec.size() > k);
	int nPts = f_vec.size();
	//std::cout << "#feature vectors " << nPts << std::endl;
	assert(f_vec[0].size() > 0);
	int dim = f_vec[0].size();
	//std::cout << "Dim of a feature vector " << dim << std::endl;

	// termination parameters
	KMterm  term(100, 0, 0, 0,						// run for 100 stages
		0.10, 0.10, 3,								// other typical parameter values 
		0.50, 10, 0.95);

	KMdata dataPts(dim, nPts);						// allocate data storage
	setVectorToData(f_vec, dataPts);
	//std::cout << "build KcTree... " << std::endl;
	dataPts.buildKcTree();								// build filtering structure
	KMfilterCenters	ctrs(k, dataPts);				// allocate centers

	//KMlocalLloyds	kmAlg(ctrs, term);				// repeated Lloyd's
													// you can use the following methods:
													// KMlocalLloyds, KMlocalSwap, KMlocalEZ_Hybrid, KMlocalHybrid
	KMlocalHybrid	kmAlg(ctrs, term);

	//std::cout << "Compute Kmeans... " << std::endl;
	ctrs = kmAlg.execute();							// execute

													// print number of stages
	//std::cout << "Number of stages: " << kmAlg.getTotalStages() << "\n";

	// print average distortion
	//cout << "Average distortion: " << ctrs.getDist() / nPts << "\n";

	//ctrs.print();				// print final centers
								// get/print final cluster assignments


	KMctrIdxArray closeCtr = new KMctrIdx[dataPts.getNPts()];
	double* sqDist = new double[dataPts.getNPts()];
	ctrs.getAssignments(closeCtr, sqDist);

	double* dists = ctrs.getDists();

	ClusCenters.clear();
	for (int i = 0; i < k; i++)
	{
		std::vector<double> tmpCtr;
		for (int j = 0; j < dim; j++)
		{
			tmpCtr.push_back(ctrs[i][j]);
		}
		ClusCenters.push_back(tmpCtr);
	}

	// print closeCtr
	Labels2Clusters.clear(); Labels2Clusters.resize(k);
	ClusterId.clear();
	for (int i = 0; i < dataPts.getNPts(); i++)
	{
		ClusterId.push_back(closeCtr[i]);
		Labels2Clusters[closeCtr[i]].push_back(i);
	}
	return;

}

void CKmeansWrapper::setVectorToData(std::vector<std::vector<double> > f_vec, KMdata & dataPts)
{
	// feature_vec.size = nPts
	int dim = f_vec[0].size();
	for (int i = 0; i < f_vec.size(); i++)
	{
		for (int j = 0; j < dim; j++)
		{
			// feature_ is int, while dataPts[i][j] is double
			dataPts[i][j] = f_vec[i][j];
		}
	}
}