#ifndef C_SKELETON_H
#define C_SKELETON_H

#include "../DataColle/custom_openmesh_type.h"
#include "def_alg_prereq.h"
#include "def_alg_typedef.h"


class DEF_ALGCOLLE_CLASS CCrossSection
{
public:
	CCrossSection() { };
	~CCrossSection() { };
	CCrossSection(const CCrossSection &b) {
		prof_pts_ = b.GetProfPts();
		emb_prof_pts_ = b.GetEmbProfPts();
		def_emb_prof_pts_ = b.GetDefEmbProfPts();
		skel_id_ = b.GetSid();
		is_closed_ = b.IsClosed();
		is_deformed_ = b.IsDeformed();
	};

	void CopyFrom(const CCrossSection &b) {
		prof_pts_ = b.GetProfPts();
		emb_prof_pts_ = b.GetEmbProfPts();
		skel_id_ = b.GetSid();
		is_closed_ = b.IsClosed();
	}

	void SetProfPts(const std::vector<COpenMeshT::Point> right) {
		prof_pts_ = right;
		//update_emb_prof_pts_with_new_prof_pts();
	};
	std::vector<COpenMeshT::Point> GetProfPts() const {
		return prof_pts_;
	};

	void SetEmbProfPts(const std::vector<COpenMeshT::Point> right) {
		emb_prof_pts_ = right;
		//update_prof_pts_with_new_emb_prof_pts();
	};
	std::vector<COpenMeshT::Point> GetEmbProfPts() const {
		return emb_prof_pts_;
	};

	void SetDefEmbProfPts(const std::vector<COpenMeshT::Point> right) {
		def_emb_prof_pts_ = right;
		//update_prof_pts_with_new_emb_prof_pts();
	};
	std::vector<COpenMeshT::Point> GetDefEmbProfPts() const {
		return def_emb_prof_pts_;
	};

	void SetSid(int id) {
		skel_id_ = id;
	};
	int GetSid() const {
		return skel_id_;
	}

	void SetClosed(bool b = true) {
		is_closed_ = b;
	}
	bool IsClosed() const {
		return is_closed_;
	}

	void SetDeformed(bool b = true) {
		is_deformed_ = b;
	}
	bool IsDeformed() const {
		return is_deformed_;
	}
	//void Print();

private:
	std::vector<COpenMeshT::Point> prof_pts_;
	std::vector<COpenMeshT::Point> emb_prof_pts_;
	std::vector<COpenMeshT::Point> def_emb_prof_pts_;
	int skel_id_;
	//COpenMeshT::Point normal_;
	bool is_closed_;
	bool is_deformed_;

private:
	void update_prof_pts_with_new_emb_prof_pts();
	void update_emb_prof_pts_with_new_prof_pts();
};

class DEF_ALGCOLLE_CLASS CSkeleton
{
public:
	CSkeleton(const std::vector<COpenMeshT::Point> skeletal_pts) {
		skeletal_pts_ = skeletal_pts;
		compute_rotation_minimizing_frames();
		compute_accumulated_arc_length();
	};
	CSkeleton() { };
	~CSkeleton() { };
	CSkeleton(const CSkeleton &b) {
		skeletal_pts_.clear();
		skeletal_pts_ = b.GetSkeletalPts();
		accum_arclength_ = b.GetAccumArcLength();
		b.GetRMF(RMF_list_);
	};

	void CopyFrom(CSkeleton &b) {
		skeletal_pts_.clear();
		skeletal_pts_ = b.GetSkeletalPts();
		accum_arclength_ = b.GetAccumArcLength();
		b.GetRMF(RMF_list_);
	}

	void InitData(const std::vector<COpenMeshT::Point> skeletal_pts) {
		skeletal_pts_.clear();
		skeletal_pts_ = skeletal_pts;
		compute_rotation_minimizing_frames();
		compute_accumulated_arc_length();
	};
	std::vector<COpenMeshT::Point> GetSkeletalPts() const {
		return skeletal_pts_;
	};
	std::vector<double> GetAccumArcLength() const {
		return accum_arclength_;
	}
	void GetRMF(std::vector<Mat3d> &RMF_list) const {
		for (int i = 0; i < RMF_list_.size(); i++)
			RMF_list.push_back(RMF_list_[i]);
	}
	//void Print();

private:
	std::vector<COpenMeshT::Point> skeletal_pts_;
	std::vector<Mat3d> RMF_list_;
	std::vector<double> accum_arclength_;

	// functions
	void compute_accumulated_arc_length();
	void compute_rotation_minimizing_frames();
	void double_reflection(const DenseMatrixXd & skeletal_pts,
		const DenseMatrixXd & T,
		DenseMatrixXd & R,
		DenseMatrixXd & S);
};

#endif // !C_SKELETON_H
