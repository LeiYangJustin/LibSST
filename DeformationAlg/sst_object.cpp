#include "sst_object.h"
#include "../GeneralTools/sorting_2d_pts.h"

void CSstObject::LocalDeformSetup()
{
	assert(has_local_deformation_);
	std::cout << "Time spent for Local Setup: " << std::endl;
	auto t1 = std::chrono::high_resolution_clock::now();
	// SOME CODES
	std::vector<CCrossSection> def_cross_sections;
	for (int i = 0; i < def_csid_list_.size(); i++)
	{
		int id = def_csid_list_[i];
		def_cross_sections.push_back(cross_sections_[id]);
	}
	deformer_->LocalDeformationSetUp(def_cross_sections);
	// SOME CODES
	auto t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	std::cout << time_span.count() << " seconds" << std::endl;

	//std::vector<CCrossSection> def_cross_sections;
	//for (int i = 0; i < def_csid_list_.size(); i++)
	//{
	//	int id = def_csid_list_[i];
	//	def_cross_sections.push_back(cross_sections_[id]);
	//}
	//deformer_->LocalDeformationSetUp(def_cross_sections);
}

void CSstObject::GlobalDeformSetup(bool use_local_setup)
{
	assert(has_global_deformation_);
	if (use_local_setup) {
		std::cerr << "The locally supported RBF is not available now" << std::endl;
		//// local support
		//deformer_->SetUseLocalToSolveGlobal(use_local_setup);
		//deformer_->LocalDeformationSetUp(def_cross_sections);
	}
	else {
		//deformer_->GlobalDeformationSetup(cross_sections_);

		std::cout << "Time spent for Global Setup: " << std::endl;
		auto t1 = std::chrono::high_resolution_clock::now();
		// SOME CODES
		deformer_->GlobalDeformationSetup(cross_sections_);
		auto t2 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
		std::cout << time_span.count() << " seconds" << std::endl;
	}
}

void CSstObject::LocalDeformSolve()
{
	assert(has_local_deformation_);
	std::vector<CCrossSection> def_cross_sections;
	for (int i = 0; i < def_csid_list_.size(); i++)
	{
		int id = def_csid_list_[i];
		def_cross_sections.push_back(cross_sections_[id]);
	}
	this->def_trimesh_ = this->trimesh_;
	deformer_->LocalDeformationSolve(def_cross_sections, &this->trimesh_, &this->def_trimesh_);
	decoding_mesh(); // have the deformed mesh in the original space
	is_deformed_ = true;
}


void CSstObject::GlobalDeformSolve()
{
	assert(has_global_deformation_);
	// use global support
	if (!deformer_->GetUseLocalToSolveGlobal())
	{
		//deformer_->GlobalDeformationSolve(cross_sections_, &this->trimesh_, &this->def_trimesh_);

		std::cout << "Time spent for Global Solve: " << std::endl;
		auto t1 = std::chrono::high_resolution_clock::now();
		// SOME CODES
		def_trimesh_ = trimesh_;
		std::cout << trimesh_.n_vertices() << std::endl;
		std::cout << def_trimesh_.n_vertices() << std::endl;
		deformer_->GlobalDeformationSolve(cross_sections_, &this->trimesh_, &this->def_trimesh_);
		auto t2 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
		std::cout << time_span.count() << " seconds" << std::endl;
	}
	// use local support
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
		COpenMeshT::Point tang_vec(rmf_list[sid](0, 0), rmf_list[sid](0, 1), rmf_list[sid](0, 2));
		extracting_single_cross_section(skeleton_.GetSkeletalPts()[sid], tang_vec, cs_pts);
		extracting_single_cross_section(cemb, COpenMeshT::Point(1.0, 0.0, 0.0), emb_cs_pts, true);
		// init cs
		CCrossSection cs;
		cs.SetSid(sid);
		cs.SetProfPts(cs_pts);
		cs.SetEmbProfPts(emb_cs_pts);
		cs.SetClosed(false);
		cs.SetDeformed(false);
		//cs.SetProfPts();
		cross_sections_.push_back(cs);
	}
}

void CSstObject::extracting_single_cross_section(COpenMeshT::Point center, 
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

		COpenMeshT::Point p_from, p_to;
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
			cnt++;
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

	// try CGAL::Optimal Transportation Curve Reconstruction
	// sorting the vertices
	CGeoCalculator::pts_sorting_alg(cs_pts);
	CGeoCalculator::simplify_polygon(cs_pts);
	if (has_global_deformation_)
		CGeoCalculator::sample_polygon(cs_pts, 0.02, true);
	else
		CGeoCalculator::sample_polygon(cs_pts, 0.01, true);
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

	std::vector<CCrossSection> def_cs_list;
	for (int i = 0; i < cross_sections_.size(); i++)
	{
		CCrossSection cs = cross_sections_[i];
		
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
		cross_sections_[i] = cs;
	}
}

void CSstObject::set_parameters()
{
	fname_src_skeleton_ = "fDataRepo/src_skeleton.txt";
	fname_dst_skeleton_ = "fDataRepo/dst_skeleton.txt";
	
	fname_src_cross_sections_ = "fDataRepo/CSdata/src_cross_sections";
	fname_dst_cross_sections_ = "fDataRepo/CSdata/dst_cross_sections";

	fname_src_emb_cross_sections_ = "fDataRepo/src_emb_cross_sections.txt";
	fname_dst_emb_cross_sections_ = "fDataRepo/dst_emb_cross_sections.txt";

	fname_src_mesh_ = "fDataRepo/src_mesh.txt";
	fname_dst_mesh_ = "fDataRepo/dst_mesh.txt";

	fname_src_emb_mesh_ = "fDataRepo/src_emb_mesh.txt";
	fname_dst_emb_mesh_ = "fDataRepo/dst_emb_mesh.txt";
}

void CSstObject::ResetCrossSectionsAfterDeformation()
{
	for (int i = 0; i < cross_sections_.size(); i++)
		cross_sections_[i].SetDeformed(false);
}

bool CSstObject::PrintSSTandMesh()
{
	// print src
	bool is_emb_cs = true;
	if (!PrintSkeleton(fname_src_skeleton_, skeleton_))
		return false;
	if (!PrintCrossSectionList(fname_src_cross_sections_, cross_sections_, is_emb_cs))
		return false;
	if (!PrintMesh(fname_src_mesh_, fname_src_emb_mesh_))
		return false;

	std::vector<CCrossSection> def_cross_sections;
	for (int i = 0; i < cross_sections_.size(); i++)
	{
		if (cross_sections_[i].IsDeformed()) {
			def_cross_sections.push_back(cross_sections_[i]);
		}
	}

	// print def
	if (is_deformed_) {
		if (!PrintSkeleton(fname_dst_skeleton_, def_skeleton_))
			return false;
		if (!PrintCrossSectionList(fname_dst_cross_sections_, def_cross_sections))
			return false;
		if (!PrintMesh(fname_dst_mesh_, fname_dst_emb_mesh_))
			return false;
	}

	return true;
}

bool CSstObject::PrintMesh(std::string fname_original, std::string fname_embedded, bool is_deformed)
{
	bool is_embedded = true;
	if (!is_deformed) {
		if (!PrintMesh(fname_original, trimesh_, !is_embedded))
			return false;
		if (!PrintMesh(fname_embedded, trimesh_, is_embedded))
			return false;
	}
	else {
		if (!PrintMesh(fname_original, def_trimesh_, !is_embedded))
			return false;
		if (!PrintMesh(fname_embedded, def_trimesh_, is_embedded))
			return false;
	}
	
	return true;
}

bool CSstObject::PrintMesh(std::string fname, const COpenMeshT &m, bool is_emb)
{
	// open file
	std::ofstream  fout_file;
	fout_file.open(fname);
	if (!fout_file.is_open()) {
		std::cerr << "PrintMesh(): cannot open the file" << std::endl;
		return false;
	}
	// check if there is any embedded mesh
	if (is_emb && !is_encoded_)
	{
		std::cerr << "PrintMesh(): no embedded mesh" << std::endl;
		return false;
	}
	// write
	fout_file << "##Mesh" << std::endl;
	fout_file << "#Vertices" << std::endl;
	if (!is_emb) {
		for (auto viter = m.vertices_begin();
			viter != m.vertices_end(); ++viter)
		{
			fout_file << viter->idx() << " " 
				<< m.point(*viter)
				<< std::endl;
		}
	}
	else {
		for (auto viter = m.vertices_begin();
			viter != m.vertices_end(); ++viter)
		{
			fout_file << viter->idx() << " " 
				<< m.data(*viter).get_emb_coord()
				<< std::endl;
		}
	}
	fout_file << "#Faces" << std::endl;
	for (auto fiter = m.faces_begin();
		fiter != m.faces_end(); ++fiter)
	{
		for (auto ccwfv_iter = m.cfv_ccwbegin(*fiter);
			ccwfv_iter != m.cfv_ccwend(*fiter); ++ccwfv_iter)
		{
			fout_file << ccwfv_iter->idx() << " ";
		}
		fout_file << std::endl;
	}
	fout_file.close();
	return true;
}

bool CSstObject::PrintSkeleton(std::string fname, const CSkeleton &s)
{
	// open file
	std::ofstream  fout_file;
	fout_file.open(fname);
	if (!fout_file.is_open()) {
		std::cerr << "PrintSkeleton(): cannot open the file" << std::endl;
		return false;
	}

	// data
	std::vector<COpenMeshT::Point> skel_pts = s.GetSkeletalPts();
	std::vector<double> alength = s.GetAccumArcLength();
	std::vector<Mat3d> rmf_list; s.GetRMF(rmf_list);
	// assert
	if (skel_pts.size() != alength.size() || skel_pts.size() != rmf_list.size())
	{
		std::cerr << "PrintSkeleton(): data size do not match" << std::endl;
		return false;
	}
	

	// write
	for (int i = 0; i < skel_pts.size(); i++) {
		fout_file << skel_pts[i] <<  std::endl;
	}
	//fout_file << "##Skeleton" << std::endl;
	//fout_file << "#Points" << std::endl;
	//for (int i = 0; i < skel_pts.size(); i++) {
	//	fout_file << skel_pts[i] << ", " << alength[i] <<  std::endl;
	//}
	//fout_file << "#RMF" << std::endl;
	//for (int i = 0; i < rmf_list.size(); i++) {
	//	fout_file << rmf_list[i].row(0) << ", "
	//		<< rmf_list[i].row(1) << ", "
	//		<< rmf_list[i].row(2) << ", "
	//		<< std::endl;
	//}
	fout_file.close();
	return true;
}

bool CSstObject::PrintCrossSectionList(std::string fname, const std::vector<CCrossSection> &cs_list, bool is_emb)
{
	// write
	// open file
	for (int i = 0; i < cs_list.size(); i++) {
		CCrossSection cs = cs_list[i];

		std::string fname_cs = fname;
		std::ostringstream id_str;
		//id_str << cs.GetSid();
		id_str << i;
		fname_cs.append(id_str.str());
		fname_cs.append(".txt");

		std::ofstream fout_file;
		fout_file.open(fname_cs);
		std::vector<COpenMeshT::Point> cs_pts;
		if (is_emb) {
			cs_pts = cs_list[i].GetEmbProfPts();
		}
		else {
			cs_pts = cs_list[i].GetProfPts();
		}
		for (int j = 0; j < cs_pts.size(); j++) {
			fout_file << cs_pts[j] << std::endl;
		}
		fout_file.close();
	}
	return true;
	
	//// make data
	//std::map<int, COpenMeshT::Point> id_pt_map;
	//std::vector<std::vector<int>> cs_ids_list;
	//int id = 0;
	//for (int i = 0; i < cs_list.size(); i++)
	//{
	//	std::vector<COpenMeshT::Point> cs_pts;
	//	if (is_emb) {
	//		cs_pts = cs_list[i].GetEmbProfPts();
	//	}
	//	else {
	//		cs_pts = cs_list[i].GetProfPts();
	//	}
	//	std::vector<int> cs_ids;
	//	for (int j = 0; j < cs_pts.size(); j++) {
	//		id_pt_map.insert(std::make_pair(id, cs_pts[j]));
	//		cs_ids.push_back(id);
	//		id++;
	//	}
	//	cs_ids_list.push_back(cs_ids);
	//}
	//std::ofstream  fout_file;
	//fout_file.open(fname);
	//if (!fout_file.is_open()) {
	//	std::cerr << "PrintCrossSectionList(): cannot open the file" << std::endl;
	//	return false;
	//}
	//fout_file << "##CrossSection" << std::endl;
	//fout_file << "#Points" << std::endl;
	//for (auto itmap = id_pt_map.begin(); itmap != id_pt_map.end(); ++itmap) {
	//	fout_file << itmap->first << ", " << itmap->second << std::endl;
	//}
	//for (int i = 0; i < cs_ids_list.size(); i++) {
	//	fout_file << "#CS " << i << std::endl;
	//	for (int j = 0; j < cs_ids_list[i].size(); j++) {
	//		fout_file << cs_ids_list[i][j] << " ";
	//	}
	//	fout_file << std::endl;
	//}
	//fout_file.close();
	//return true;
}
