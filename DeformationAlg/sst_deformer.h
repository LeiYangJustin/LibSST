#ifndef C_SST_DEFORM_H
#define C_SST_DEFORM_H

#include "def_alg_prereq.h"
#include "../DataColle/mesh_object.h"
//#include "sst_object.h"
#include "skeleton_object.h"
#include "mesh_deformation.h"

enum DeformType{Global, Local, NotGiven};

class DEF_ALGCOLLE_CLASS CSSTDeformer
{
public:
	CSSTDeformer();
	~CSSTDeformer();

	void SetUseLocalToSolveGlobal(bool b = true) {
		use_local_to_solve_global = b;
	}
	bool GetUseLocalToSolveGlobal() const {
		return use_local_to_solve_global;
	}

	// local setup and solve
	void LocalDeformationSetUp(
		const std::vector<CCrossSection> def_cs_list);
	bool LocalDeformationSolve(
		const std::vector<CCrossSection> def_cs_list,
		COpenMeshT *mesh,
		COpenMeshT *def_mesh);

	// global setup and solve
	void GlobalDeformationSetup(std::vector<CCrossSection> def_cs_list);
	bool GlobalDeformationSolve(
		const std::vector<CCrossSection> def_cs_list,
		COpenMeshT *mesh,
		COpenMeshT *def_mesh);

	void Update()
	{
		p_deformer_map_.clear();
	}

private:
	std::map<std::pair<int, int>, CMeshDeformation*> p_deformer_map_;
	std::vector<std::pair<int, int>> cs_pair_list_;
	bool use_local_to_solve_global;
	
	DeformType deform_type_; // 1-global; 2-local; 0-not initialized

private:
	void get_deformer(const std::vector<COpenMeshT::Point> left_cs_pts,
		const std::vector<COpenMeshT::Point> right_cs_pts,
		CMeshDeformation * deformer);
	void get_deformer(const std::vector<COpenMeshT::Point> all_cs_pts,
		CMeshDeformation * deformer);

	void local_deformation_solve(
		const std::vector<COpenMeshT::Point> left_handles_vecs,
		const std::vector<COpenMeshT::Point> right_handles_vecs,
		const std::vector<COpenMeshT::Point> roi_pts,
		CMeshDeformation* p_deformer,
		std::vector<COpenMeshT::Point> &roi_vecs);

	void global_deformation_solve(
		const std::vector<COpenMeshT::Point> all_handles_vecs,
		const std::vector<COpenMeshT::Point> roi_pts,
		CMeshDeformation* p_deformer,
		std::vector<COpenMeshT::Point> &roi_vecs);

	//void get_ROI_pids(std::pair<int, int> csid_pair, std::vector<int> roi_ids);


};
#endif
