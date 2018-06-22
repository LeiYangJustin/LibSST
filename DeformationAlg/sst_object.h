#ifndef SST_OBJECT_H
#define SST_OBJECT_H

#include "../DataColle/custom_openmesh_type.h"
#include "def_alg_prereq.h"
#include "def_alg_typedef.h"
#include "skeleton_object.h"
#include "sst_deformer.h"

class DEF_ALGCOLLE_CLASS CSstObject
{
public:
	CSstObject(CSkeleton s, COpenMeshT& in_mesh) :
		has_global_deformation_(false), has_local_deformation_(false),
		is_encoded_(false), is_deformed_(false)
	{ 
		skeleton_ = s; trimesh_ = in_mesh;
		deformer_ = new CSSTDeformer;
		set_parameters();
	};
	CSstObject() :
		has_global_deformation_(false), has_local_deformation_(false),
		is_encoded_(false), is_deformed_(false)
	{
		deformer_ = new CSSTDeformer;
		set_parameters();
	};
	~CSstObject() 
	{ 
		if (deformer_ != NULL) 
			delete deformer_; 
	};
	void Encode() 
	{ 
		encoding_mesh(); 
		is_encoded_ = true;
	};
	void InitCrossSections(std::vector<int> sid_list) 
	{
		extracting_cross_sections(sid_list);
	};
	void SetDefSkeleton(CSkeleton def_skel)
	{
		has_global_deformation_ = true;
		has_local_deformation_ = !has_global_deformation_;
		def_skeleton_.CopyFrom(def_skel);

		// deform the cross-sections
		cross_section_transformed_due_to_skeleton_change();
	};
	void RenewCSInfo(int csid, std::vector<COpenMeshT::Point> cs_pts)
	{
		cross_sections_[csid].SetEmbProfPts(cs_pts);
	};
	// return false if #def_cs_pts != #cs_pts
	bool SetDefCSInfo(int csid, std::vector<COpenMeshT::Point> def_cs_pts)
	{
		has_local_deformation_ = true;
		has_global_deformation_ = !has_local_deformation_;
		cross_sections_[csid].SetDeformed(true);
		cross_sections_[csid].SetDefEmbProfPts(def_cs_pts);
		if (def_cs_pts.size() == cross_sections_[csid].GetEmbProfPts().size())
			return true;
		else
			return false;
	};
	void OutputDefMesh(COpenMeshT &mesh) 
	{
		if (is_deformed_)
			mesh = def_trimesh_;
	};

	// deformation call
	bool LocalDeformSetup(/*std::vector<std::pair<int, int>> cs_pair_list*/);
	bool GlobalDeformSetup(bool use_local_setup = false);

	bool LocalDeformSolve();
	bool GlobalDeformSolve();

	//void SetParameters(std::string fname_src_skeleton,
	//	std::string fname_dst_skeleton,
	//	std::string fname_src_cross_sections,
	//	std::string fname_dst_cross_sections,
	//	std::string fname_src_mesh,
	//	std::string fname_src_emb_mesh,
	//	std::string fname_dst_mesh,
	//	std::string fname_dst_emb_mesh);

	// print
	bool PrintSSTandMesh();
	bool PrintMesh(std::string fname_original, std::string fname_embedded, bool is_deformed = false);
	bool PrintMesh(std::string fname, const COpenMeshT &m, bool is_emb = false);
	bool PrintSkeleton(std::string fname, const CSkeleton &s);
	bool PrintCrossSectionList(std::string fname, const std::vector<CCrossSection> &cs_list, bool is_emb = false);

private:
	// deformer
	CSSTDeformer* deformer_;
	//
	CSkeleton skeleton_;
	CSkeleton def_skeleton_;
	//
	std::vector<CCrossSection> cross_sections_;
	std::vector<CCrossSection> def_cross_sections_;
	//
	COpenMeshT trimesh_;
	COpenMeshT def_trimesh_;
	//
	bool is_encoded_;
	bool is_deformed_;
	bool has_global_deformation_;
	bool has_local_deformation_;

	// path
	std::string fname_src_skeleton_;
	std::string fname_dst_skeleton_;
	std::string fname_src_cross_sections_;
	std::string fname_dst_cross_sections_;
	std::string fname_src_emb_cross_sections_;
	std::string fname_dst_emb_cross_sections_;
	std::string fname_src_mesh_;
	std::string fname_dst_mesh_;
	std::string fname_src_emb_mesh_;
	std::string fname_dst_emb_mesh_;

private:
	// functions
	//void binding_skeleton_mesh();
	void encoding_mesh();
	void extracting_cross_sections(std::vector<int> sid_list);
	void extracting_single_cross_section(
		COpenMeshT::Point c, COpenMeshT::Point d, 
		std::vector<COpenMeshT::Point> &cs_pts);
	void decoding_vector_field(DenseMatrixXd & U,
		const DenseMatrixXd & V, const std::vector<int> ids);
	void decoding_mesh();
	void cross_section_transformed_due_to_skeleton_change();

	void set_parameters();


};
#endif
