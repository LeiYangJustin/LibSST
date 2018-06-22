#ifndef C_SST_PARAMETERIZATION_H
#define C_SST_PARAMETERIZATION_H

#include "def_alg_prereq.h"
#include "def_alg_typedef.h"
#include "../DataColle/sst_object.h"
#include "../DataColle/mesh_object.h"
#include "../GeneralTools/geo_calculator.h"

class DEF_ALGCOLLE_CLASS CSSTparameterization
{
public:
	CSSTparameterization();
	~CSSTparameterization();

	// set up the parameterization
	void SetUpAndEncode(COpenMeshT &mesh, std::vector<COpenMeshT::Point> skel_pts);

	// set and get data
	void SetSkeleton(CSkeleton s) {
		skeleton_ = s;
	}
	CSkeleton GetSkeleton() {
		return skeleton_;
	}
	void SetCrossSectionList(std::vector<CCrossSection> cs_list) {
		cross_sections_ = cs_list;
	}
	std::vector<CCrossSection> GetAllCrossSections() {
		return cross_sections_;
	}
	void PushBackCrossSection(CCrossSection cs) {
		cross_sections_.push_back(cs);
	}
	void GetCrossSectionWithId(int id, CCrossSection & cs) {
		if (id < cross_sections_.size())
			cs = cross_sections_[id];
	}

	//void SkelDeformation(DenseMatrixXd & new_pts, const DenseMatrixXd & def_skel_pts);
	void GetDeformedCrossSectionViaSkelDeformation(
		std::vector<DenseMatrixXd> & def_vecs_list,
		std::vector<std::vector<COpenMeshT::Point>> & def_cs_pts_list,
		const std::vector<std::vector<COpenMeshT::Point>> & emb_cs_pts_list,
		const DenseMatrixXd & def_skel_pts,
		const std::vector<int> cs_id_list);


	// DECODE DATA FROM THE SKELETON-INDUCED EMBEDDING
	void DecodeVectorField(DenseMatrixXd & U, const DenseMatrixXd & V);
	void DecodeVectorField(DenseMatrixXd & U, const DenseMatrixXd & V, std::vector<int> ids);
	void DecodePositions(DenseMatrixXd & new_pts);
	
	void DecodePosition(Vector3d & new_pt, 
		const Vector3d & old_pt, const int id_min);

	// EXTRACT CROSS-SECTIONS
	void ExtractSingleCrossSection(std::vector<COpenMeshT::Point> &pt_list,
		COpenMeshT::Point center, COpenMeshT::Point normal, bool is_encoded =false);

	// GET EMBEDDED PTS FOR DEFORMATION
	void GetEmbeddedPts(DenseMatrixXd & emb_pts) { emb_pts = emb_pts_; };
	void GetEmbeddedMesh(COpenMeshT & embMesh) { embMesh = embMesh_; };
	void GetPidSkeidMap(std::map<int, int> & pid_skeid_map) { pid_skeid_map = pid_skeid_map_; };
	void GetPidVidMap(std::map<int, int> & pid_vid_map) { pid_vid_map = pid_vid_map_; };
	void GetUnitTangentVectors(std::vector<COpenMeshT::Point> &tang_vecs);

	void computeSkeletalRMF(const DenseMatrixXd & skeletal_pts, std::vector<Mat3d> &rmf_list);

private:
	// mesh
	COpenMeshT triMesh_;
	COpenMeshT embMesh_;

	CSkeleton skeleton_;
	std::vector<CCrossSection> cross_sections_;

	// pts
	DenseMatrixXd ori_pts_;
	DenseMatrixXd emb_pts_;
	
	// skel
	DenseMatrixXd skel_pts_;
	std::vector<Mat3d> rmf_list_;

	// map
	std::map<int, int> vid_pid_map_;
	std::map<int, int> pid_vid_map_;
	std::map<int, int> pid_skeid_map_;

	bool is_encoded_;
private:
	// Compute space embedding
	
	void doubleReflection(const DenseMatrixXd & skeletal_pts,
		const DenseMatrixXd & T,
		DenseMatrixXd & R,
		DenseMatrixXd & S);

	// ENCODE MESH TO THE SKELETON-INDUCED EMBEDDING
	void Encode();

	//void getEmbeddedCrossSections();
};

#endif // !C_SST_PARAMETERIZATION_H



