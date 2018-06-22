#ifndef C_SST_DEFORM_WRAPPER_H
#define C_SST_DEFORM_WRAPPER_H

#include "def_alg_prereq.h"
#include "def_alg_typedef.h"
#include "../DataColle/custom_openmesh_type.h"
#include "../DataColle/mesh_object.h"
#include "../GeneralTools/geo_calculator.h"
#include "../GeneralTools/fe_meshIO.h"
#include "sst_parameterization.h"
#include "mesh_deformation.h"

class DEF_ALGCOLLE_CLASS CSSTDeformWrapper
{
public:
	CSSTDeformWrapper();
	~CSSTDeformWrapper();

	void SetUp(COpenMeshT triMesh, std::vector<COpenMeshT::Point> skelPts);
	void InitSST();
	void DeformWithSST(COpenMeshT & out_mesh,
		std::vector<std::vector<COpenMeshT::Point>> & ori_cs_pts_list,
		std::vector<int> & cs_type_list);

	void DeformationLoop(COpenMeshT & out_mesh,
		std::vector<std::vector<COpenMeshT::Point>> & ori_cs_pts_list,
		std::vector<int> & cs_type_list,
		const std::string &fpath,
		double scale = 1.0);

	void GlobalDeformationLoop(
		COpenMeshT & out_mesh, 
		std::vector<std::vector<COpenMeshT::Point>> & def_cs_pts_list, // plot
		const std::string &fpath,
		double scale = 1.0);

	void outputEmbeddingMesh(std::string fpath, double scale = 1.0);

	

private:
	// read-in
	COpenMeshT triMesh_; 
	std::vector<COpenMeshT::Point> skelPts_;

	// obj
	CSSTparameterization * sst_obj_; // construct
	std::vector<CMeshDeformation*> p_deformer_list_;

	// meshes
	COpenMeshT embMesh_; // get from sst

	// maps: 
	// this is not needed if we can ensure there is no topo-change in the original mesh
	std::map<int, int> pid_vid_map_; // get from sst
	std::map<int, int> pid_skeid_map_; // get from sst

private:
	void DeformationSetup(const DenseMatrixXd & cspair_map,
		const std::vector<std::vector<COpenMeshT::Point>> & cs_pts_list);

	void GlobalDeformationSetup(
		const std::vector<std::vector<COpenMeshT::Point>> & cs_pts_list);

	void DeformationSolveAndOutput(int sample_id,
		COpenMeshT & out_mesh,
		const DenseMatrixXd & cspair_map,
		const std::vector<std::vector<int>> & rid_keeper_list,
		const std::vector<int> &fix_mask,
		const std::string &fpath);

	void GlobalDeformationSolveAndOutput(int sample_id,
		COpenMeshT & out_mesh,
		DenseMatrixXd & vecs_handles,
		const std::string &fpath);

	void GlobalDeformationSolveAndOutput(int sample_id,
		COpenMeshT & out_mesh,
		const DenseMatrixXd & cspair_map,
		const std::vector<DenseMatrixXd> & vecs_array,
		const std::vector<std::vector<int>> & rid_keeper_list,
		const std::string &fpath);
	
	void getCrossSections(
		std::vector<std::vector<COpenMeshT::Point>> & emb_cs_pts_list,
		std::vector<std::vector<COpenMeshT::Point>> & cs_pts_list,
		std::vector<int> & id_cs);

	void getROIPts(int numROI, 
		std::vector<int> id_array, 
		DenseMatrixXd & cspair_map, 
		std::vector<std::vector<int>> &rid_keeper_list);

	bool readCSdisplacement(
		const std::string filename,
		std::vector<COpenMeshT::Point> &con_pts,
		std::vector<COpenMeshT::Point> &con_vecs);

	bool readCSData(
		const std::string filename,
		std::vector<COpenMeshT::Point> &vec3d);
	bool readCSData(
		const std::string filename,
		std::vector<double> &data);

	bool readSkelCtrlPts(
		const std::string filename, 
		std::vector<double> &data);

	bool writeMeshVertsToTxt(COpenMeshT & out_mesh, int sample_id);
};

#endif // !C_SST_DEFORM_WRAPPER_H


