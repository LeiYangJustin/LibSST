#ifndef CKMEANS_CLUSTERING_CLASS_H
#define CKMEANS_CLUSTERING_CLASS_H

#include "general_tool_prereq.h"

class KMdata;

class GENERAL_TOOLS_CLASS CKmeansWrapper
{
public:
	CKmeansWrapper();
	~CKmeansWrapper();

	//CKmeansWrapper(int k, int dim, int nPts);

	void executeKMeans(
		int k, 
		std::vector<std::vector<double>> feature_vec, 
		std::vector<int> &ClusterId, 
		std::vector<std::vector<double> > &ClusCenters);
	
	void executeKMeans(
		int k,
		std::vector<std::vector<double>> feature_vec,
		std::vector<int> &ClusterId,
		std::vector<std::vector<double> > &ClusCenters,
		std::vector<std::vector<int> > &Labels2Clusters);

	void setTermParameters();

private:
	void setVectorToData(std::vector<std::vector<double>> feature_vec, KMdata & dataPts);
	

};

#endif