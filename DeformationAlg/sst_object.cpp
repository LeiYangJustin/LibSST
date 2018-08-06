#include "sst_object.h"
#include "../GeneralTools/geo_calculator.h"
#include "../GeneralTools/euclidean_MST.h"
//#include "../GeneralTools/ransac_wrapper.h"

void CSstObject::LocalDeformSetup()
{
	assert(has_local_deformation_);
	
	auto t1 = std::chrono::high_resolution_clock::now();

	// DEFORMATION CODES
	std::vector<CCrossSection> def_cross_sections;
	for (int i = 0; i < def_sid_list_.size(); i++)
	{
		int id = def_sid_list_[i];
		def_cross_sections.push_back(map_id_cross_sections_[id]);
	}
	deformer_->LocalDeformationSetUp(def_cross_sections);
	// DEFORMATION CODES

	auto t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	std::cout << "Time spent for Local Setup: " << time_span.count() << " seconds" << std::endl;
}

void CSstObject::GlobalDeformSetup(bool use_local_setup)
{
	assert(has_global_deformation_);
	if (use_local_setup) {
		std::cerr << "The locally supported RBF is not available now" << std::endl;
	}
	else {
		auto t1 = std::chrono::high_resolution_clock::now();

		// DEFORMATION CODES
		std::vector<CCrossSection> def_cross_sections;
		for (auto miter = map_id_cross_sections_.begin(); miter != map_id_cross_sections_.end(); ++miter)
		{
			def_cross_sections.push_back(miter->second);
		}
		deformer_->GlobalDeformationSetup(def_cross_sections);
		// DEFORMATION CODES
		auto t2 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
		std::cout << "Time spent for Global Setup: " << time_span.count() << " seconds" << std::endl;
	}
}

bool CSstObject::LocalDeformSolve()
{
	assert(has_local_deformation_);
	std::vector<CCrossSection> def_cross_sections;
	for (int i = 0; i < def_sid_list_.size(); i++)
	{
		int id = def_sid_list_[i];
		def_cross_sections.push_back(map_id_cross_sections_[id]);
	}
	this->def_trimesh_ = this->trimesh_;
	//
	auto t1 = std::chrono::high_resolution_clock::now();
	//
	if (deformer_->LocalDeformationSolve(def_cross_sections, &this->trimesh_, &this->def_trimesh_)) {
		decoding_mesh(); // have the deformed mesh in the original space
		//
		auto t2 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
		std::cout << "Time spent for Local Solve: " << time_span.count() << " seconds" << std::endl;
		//
		is_deformed_ = true;
		return true;
	}
	else {
		std::cout << "local deformation failed" << std::endl;
		return false;
	}
		
}

bool CSstObject::GlobalDeformSolve()
{
	assert(has_global_deformation_);
	// use global support
	if (!deformer_->GetUseLocalToSolveGlobal())
	{
		def_trimesh_ = trimesh_;
		std::vector<CCrossSection> def_cross_sections;
		for (auto miter = map_id_cross_sections_.begin(); miter != map_id_cross_sections_.end(); ++miter)
		{
			def_cross_sections.push_back(miter->second);
		}
		//
		auto t1 = std::chrono::high_resolution_clock::now();
		//
		if (deformer_->GlobalDeformationSolve(def_cross_sections, &this->trimesh_, &this->def_trimesh_))
		{
			//
			auto t2 = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
			std::cout << "Time spent for Global Solve: " << time_span.count() << " seconds" << std::endl;
			//
			is_deformed_ = true;
			return true;
		}
		else {
			std::cout << "global deformation failed" << std::endl;
			return false;
		}
	}
	return false;
}

void CSstObject::encoding_mesh()
{
	std::vector<COpenMeshT::Point> skeletal_pts = skeleton_.GetSkeletalPts();
	std::vector<double> accu_arc_length = skeleton_.GetAccumArcLength();
	std::vector<Mat3d> rmf_list; skeleton_.GetRMF(rmf_list);

	// encoding
	for (auto viter = trimesh_.vertices_begin(); viter!= trimesh_.vertices_end(); viter++)
	{
		// roi label
		double dmin = std::numeric_limits<double>::max();
		int id_min = -1;
		for (int j = 0; j < skeletal_pts.size(); j++)
		{
			double d = (trimesh_.point(*viter) - skeletal_pts[j]).norm();
			if (d < dmin)
			{
				dmin = d;
				id_min = j;
			}
		}
		trimesh_.data(*viter).set_vlabel(id_min);
		// encode
		Vector3d ptmp;
		Vector3d v(trimesh_.point(*viter)[0], trimesh_.point(*viter)[1], trimesh_.point(*viter)[2]);
		Vector3d s(skeletal_pts[id_min][0], skeletal_pts[id_min][1], skeletal_pts[id_min][2]);
		ptmp = rmf_list[id_min] * (v - s);
		//ptmp(0) = accu_arc_length[id_min];
		// or ?
		ptmp(0) += accu_arc_length[id_min];
		trimesh_.data(*viter).set_emb_coord(COpenMeshT::Point(ptmp(0), ptmp(1), ptmp(2)));
	}
	is_encoded_ = true;
}

void CSstObject::extracting_cross_sections(std::vector<int> sid_list)
{
	//std::vector<int> sid_list;
	std::vector<Mat3d> rmf_list; 
	skeleton_.GetRMF(rmf_list);
	for (int i = 0; i < sid_list.size(); i++) {
		int sid = sid_list[i];
		COpenMeshT::Point cemb(skeleton_.GetAccumArcLength()[sid], 0.0, 0.0);
		std::vector<COpenMeshT::Point> cs_pts, emb_cs_pts;
		//COpenMeshT::Point tang_vec(rmf_list[sid](0, 0), rmf_list[sid](0, 1), rmf_list[sid](0, 2));
		std::cout << "cs_extraction: " << sid << "; " << cemb << std::endl;
		//extracting_single_cross_section(skeleton_.GetSkeletalPts()[sid], tang_vec, cs_pts);

		if (extracting_single_cross_section(cemb, COpenMeshT::Point(1.0, 0.0, 0.0), emb_cs_pts, true))
		{
			Vector3d skelPt(skeleton_.GetSkeletalPts()[sid][0],
				skeleton_.GetSkeletalPts()[sid][1],
				skeleton_.GetSkeletalPts()[sid][2]);
			Vector3d cembPt(cemb[0], cemb[1], cemb[2]);
			for (int j = 0; j < emb_cs_pts.size(); j++)
			{
				Vector3d ept(emb_cs_pts[j][0], emb_cs_pts[j][1], emb_cs_pts[j][2]);
				ept(0) = ept(0) - skeleton_.GetAccumArcLength()[sid];
				Vector3d pt = rmf_list[sid].transpose() * ept + skelPt;
				cs_pts.push_back(COpenMeshT::Point(pt(0), pt(1), pt(2)));
			}
			// init cs
			CCrossSection cs;
			cs.SetSid(sid);
			cs.SetProfPts(cs_pts);
			cs.SetEmbProfPts(emb_cs_pts);
			cs.SetClosed(false);
			cs.SetDeformed(false);
			map_id_cross_sections_[sid] = cs;
		}
	}
}

bool CSstObject::extracting_single_cross_section(COpenMeshT::Point center, 
	COpenMeshT::Point normal, std::vector<COpenMeshT::Point>& cs_pts, bool is_emb)
{
	cs_pts.clear();
	assert(is_encoded_);

	// circulating the mesh to find edges with two endpoints having difference sign to the plane
	int cnt = 0;
	for (auto eiter = trimesh_.edges_begin(); eiter != trimesh_.edges_end(); ++eiter)
	{
		COpenMeshT::HalfedgeHandle heh = trimesh_.halfedge_handle(*eiter, 0);
		COpenMeshT::VertexHandle vh_from = trimesh_.from_vertex_handle(heh);
		COpenMeshT::VertexHandle vh_to = trimesh_.to_vertex_handle(heh);

		COpenMeshT::Point p_from;
		COpenMeshT::Point p_to;
		if (is_emb) {
			p_from = trimesh_.data(vh_from).get_emb_coord();
			p_to = trimesh_.data(vh_to).get_emb_coord();
		}
		else {
			p_from = trimesh_.point(vh_from);
			p_to = trimesh_.point(vh_to);
		}

		double d_from = CGeoCalculator::ComputeDistFromPoint2Plane(p_from, center, normal);
		double d_to = CGeoCalculator::ComputeDistFromPoint2Plane(p_to, center, normal);

		if (d_from*d_to <= 0)
		{
			//cs_pts.push_back(p_from);
			//cs_pts.push_back(p_to);
			/*cnt++;*/
			double d_deno = fabs(d_from) + fabs(d_to);
			if (d_deno == 0 && 0 == fabs(d_to))
			{
				cs_pts.push_back(p_from);
				cs_pts.push_back(p_to);
			}
			else {
				double t = fabs(d_from) / d_deno;
				COpenMeshT::Point pt = p_from + t*(p_to - p_from);
				cs_pts.push_back(pt);
			}
		}
	}

	//std::ofstream cspts_file("cspts.txt");
	//for (int i = 0; i < cs_pts.size(); i++) {
	//	cspts_file << cs_pts[i] << std::endl;
	//}
	//cspts_file.close();

	std::vector<COpenMeshT::Point> test_pts;
	//test_pts = cs_pts;
	CfindChainLoopUsingEMST fcl_emst_obj;
	fcl_emst_obj.SetPoints(cs_pts);
	fcl_emst_obj.Compute();
	fcl_emst_obj.GetSortedPoints(test_pts);

	//
	std::cout << "extracted: " << cs_pts.size()
		<< "prunned: " << cs_pts.size() << std::endl;
	cs_pts.clear();
	cs_pts = test_pts;
	CGeoCalculator::simplify_polygon(cs_pts);

	////
	//bool use_sorting = true;
	//if (use_sorting) {
	//	if (cs_pts.size() == 0)
	//		return false;

	//	CGeoCalculator::pts_sorting_alg(cs_pts);
	//	CGeoCalculator::simplify_polygon(cs_pts);
	//	CGeoCalculator::sample_polygon(cs_pts, 0.01, false); 
	//}
	////
	return true;
}

void CSstObject::decoding_vector_field(DenseMatrixXd & U, 
	const DenseMatrixXd & V, const std::vector<int> ids)
{
	assert(V.rows() == ids.size());
	if (!is_encoded_) {
		std::cerr << "no encoding data!" << std::endl;
		return;
	}

	std::vector<Mat3d> rmf_list; skeleton_.GetRMF(rmf_list);
	U.resize(V.rows(), 3);
	for (int i = 0; i < V.rows(); i++)
	{
		Eigen::RowVectorXd v = V.row(i);
		U.row(i) = v*rmf_list[ids[i]];
	}
}

void CSstObject::decoding_mesh()
{
	std::vector<COpenMeshT::Point> skeletal_pts = skeleton_.GetSkeletalPts();
	std::vector<double> accu_arc_length = skeleton_.GetAccumArcLength();
	std::vector<Mat3d> rmf_list; 
	skeleton_.GetRMF(rmf_list);

	// decoding
	for (auto viter = def_trimesh_.vertices_begin(); viter != def_trimesh_.vertices_end(); viter++)
	{
		int sid = def_trimesh_.data(*viter).get_vlabel();
		COpenMeshT::VertexHandle vh0 = trimesh_.vertex_handle(viter->idx());
		COpenMeshT::Point p0 = trimesh_.data(vh0).get_emb_coord();
		COpenMeshT::Point p1 = def_trimesh_.data(*viter).get_emb_coord();
		Vector3d d0(p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]);
		Vector3d d1 = rmf_list[sid].transpose()*d0;
		COpenMeshT::Point ptmp = def_trimesh_.point(*viter) + COpenMeshT::Point(d1(0), d1(1), d1(2));
		def_trimesh_.set_point(*viter, ptmp);
	}
}

void CSstObject::skeleton_driven_cross_section_transformation()
{
	//std::cout << "no cross-sections are deformed according to the skeleton's change" << std::endl;

	// src skeleton
	std::vector<COpenMeshT::Point> skeletal_pts = skeleton_.GetSkeletalPts();
	std::vector<Mat3d> rmf_list;
	skeleton_.GetRMF(rmf_list);

	// def skeleton
	std::vector<COpenMeshT::Point> def_skeletal_pts = def_skeleton_.GetSkeletalPts();
	std::vector<Mat3d> def_rmf_list;
	def_skeleton_.GetRMF(def_rmf_list);

	for (auto miter = map_id_cross_sections_.begin(); miter != map_id_cross_sections_.end(); ++miter)
	{
		CCrossSection cs = miter->second;
		
		int sid = cs.GetSid();
		std::vector<COpenMeshT::Point> cspts = cs.GetProfPts();
		Vector3d psrc(skeletal_pts[sid][0], skeletal_pts[sid][1], skeletal_pts[sid][2]);
		Vector3d pdef(def_skeletal_pts[sid][0], def_skeletal_pts[sid][1], def_skeletal_pts[sid][2]);
		Mat3d Tsrc = rmf_list[sid];
		Mat3d Tdeftranspose = def_rmf_list[sid].transpose();
		for (int j = 0; j < cspts.size(); j++)
		{
			Vector3d csrc(cspts[j][0], cspts[j][1], cspts[j][2]);
			Vector3d cdef = Tdeftranspose*Tsrc * csrc - Tdeftranspose*Tsrc*psrc + pdef;
			cspts[j] = COpenMeshT::Point(cdef(0), cdef(1), cdef(2));
		}
		
		// update cs;
		cs.SetDefEmbProfPts(cspts);
		cs.SetDeformed(true);
		miter->second = cs;
	}
}

void CSstObject::ResetCrossSectionsAfterDeformation()
{
	for (auto miter = map_id_cross_sections_.begin(); miter != map_id_cross_sections_.end(); ++miter)
		miter->second.SetDeformed(false);
}

// 
void CSstObject::InsertCrossSections(std::vector<int> spid_list)
{
	// add to map
	extracting_cross_sections(spid_list);
}

void CSstObject::DeleteCrossSections(std::vector<int> spid_list)
{
	// remove from map
	for (int i = 0; i < spid_list.size(); i++)
	{
		if (map_id_cross_sections_.find(spid_list[i]) != map_id_cross_sections_.end())
		{
			map_id_cross_sections_.erase(spid_list[i]);
		}
	}
}