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
	};
	CSstObject() :
		has_global_deformation_(false), has_local_deformation_(false),
		is_encoded_(false), is_deformed_(false)
	{
		deformer_ = new CSSTDeformer;
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
		is_encoded_ = false;
	};
	CSkeleton GetSkeleton() const
	{
		return skeleton_;
	};
	void SetMesh(COpenMeshT &mesh)
	{
		trimesh_ = mesh;
		is_encoded_ = false;
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
		for (auto map_iter = map_id_cross_sections_.begin(); map_iter != map_id_cross_sections_.end(); ++map_iter)
		{
			CCrossSection cs(map_iter->second);
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

	void SetDefCSList(std::vector<int> def_sid_list) {
		has_local_deformation_ = true;
		has_global_deformation_ = !has_local_deformation_;
		def_sid_list_ = def_sid_list;
	};
	bool RenewCSInfo(int sid, std::vector<COpenMeshT::Point> cs_pts, std::vector<COpenMeshT::Point> def_cs_pts)
	{
		map_id_cross_sections_[sid].SetDeformed(true);
		map_id_cross_sections_[sid].SetEmbProfPts(cs_pts);
		map_id_cross_sections_[sid].SetDefEmbProfPts(def_cs_pts);
		if (def_cs_pts.size() == cs_pts.size())
			return true;
		else
			return false;
	};

	// deformation call
	void LocalDeformSetup(/*std::vector<std::pair<int, int>> cs_pair_list*/);
	void GlobalDeformSetup(bool use_local_setup = false);
	bool LocalDeformSolve();
	bool GlobalDeformSolve();

	void ResetCrossSectionsAfterDeformation();

private:
	// deformer
	CSSTDeformer* deformer_;
	//
	CSkeleton skeleton_;
	CSkeleton def_skeleton_;
	//
	std::vector<int> def_sid_list_;
	std::map<int, CCrossSection> map_id_cross_sections_;
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
};
#endif
