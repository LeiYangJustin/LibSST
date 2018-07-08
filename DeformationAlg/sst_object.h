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
	//
	void SetSkeleton(CSkeleton &skel) 
	{
		skeleton_.CopyFrom(skel);
	};
	CSkeleton GetSkeleton() const
	{
		return skeleton_;
	};
	void SetMesh(COpenMeshT &mesh)
	{
		trimesh_ = mesh;
	};
	COpenMeshT& GetMesh()
	{
		return trimesh_;
	};
	bool OutputDefMesh(COpenMeshT &mesh) 
	{
		if (is_deformed_) {
			mesh = def_trimesh_;
			return true;
		}
		else
			return false;
	};

	// 
	void Encode() 
	{ 
		if (is_encoded_ == false)
		{
			encoding_mesh();
			is_encoded_ = true;
		}
	};
	bool IsEncode()
	{
		return is_encoded_;
	};
	
	//
	void InitCrossSections(std::vector<int> sid_list) 
	{
		extracting_cross_sections(sid_list);
	};
	void GetCrossSections(std::vector<CCrossSection> &cs_list)
	{
		for (int i = 0; i < cross_sections_.size(); i++)
		{
			CCrossSection cs(cross_sections_[i]);
			cs_list.push_back(cs);
		}
	}
	// 
	void InsertCrossSections(std::vector<int> sid_list)
	{
		std::cout << "insert cross sections function is under construction" << std::endl;
	}

	// 
	bool IsDeformed()
	{
		return is_deformed_;
	}

	// Update
	void SetUpdate(bool b)
	{
		is_updated_ = b;
	}
	void Update()
	{
		deformer_->Update();
	}

	// Set deformation
	void SetDefSkeleton(CSkeleton def_skel)
	{
		has_global_deformation_ = true;
		has_local_deformation_ = !has_global_deformation_;
		def_skeleton_.CopyFrom(def_skel);

		// deform the cross-sections
		skeleton_driven_cross_section_transformation();
	};
	CSkeleton GetDefSkeleton()
	{
		return def_skeleton_;
	};

	void SetDefCSList(std::vector<int> def_csid_list) {
		has_local_deformation_ = true;
		has_global_deformation_ = !has_local_deformation_;
		def_csid_list_ = def_csid_list;
	};
	bool RenewCSInfo(int csid, std::vector<COpenMeshT::Point> cs_pts, std::vector<COpenMeshT::Point> def_cs_pts)
	{
		cross_sections_[csid].SetDeformed(true);
		cross_sections_[csid].SetEmbProfPts(cs_pts);
		cross_sections_[csid].SetDefEmbProfPts(def_cs_pts);
		if (def_cs_pts.size() == cs_pts.size())
			return true;
		else
			return false;
	};

	// deformation call
	void LocalDeformSetup(/*std::vector<std::pair<int, int>> cs_pair_list*/);
	void GlobalDeformSetup(bool use_local_setup = false);

	void LocalDeformSolve();
	void GlobalDeformSolve();

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
	void ResetCrossSectionsAfterDeformation();
private:
	// deformer
	CSSTDeformer* deformer_;
	//
	CSkeleton skeleton_;
	CSkeleton def_skeleton_;
	//
	std::vector<CCrossSection> cross_sections_;
	std::vector<int> def_csid_list_;
	//std::vector<CCrossSection> def_cross_sections_;
	//
	COpenMeshT trimesh_;
	COpenMeshT def_trimesh_;
	//
	bool is_encoded_;
	bool is_deformed_;
	bool is_updated_;
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
		std::vector<COpenMeshT::Point> &cs_pts, bool is_emb = false);
	void decoding_vector_field(DenseMatrixXd & U,
		const DenseMatrixXd & V, const std::vector<int> ids);
	void decoding_mesh();
	void skeleton_driven_cross_section_transformation();

	void set_parameters();

};
#endif
