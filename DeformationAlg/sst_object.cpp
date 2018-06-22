#include "sst_object.h"
#include "../GeneralTools/sorting_2d_pts.h"

bool CSstObject::LocalDeformSetup()
{
	// factorize the matrix on the left hand side
	if (has_local_deformation_)
	{
		// compute using the deformation information
		def_cross_sections_.clear();
		for (int i = 0; i < cross_sections_.size(); i++)
		{
			if (cross_sections_[i].IsDeformed()) {
				def_cross_sections_.push_back(cross_sections_[i]);
			}
		}
		deformer_->LocalDeformationSetUp(def_cross_sections_);
		return true;
	}
	else
	{
		std::cerr << "Cannot peform deformation as no information is supplied" << std::endl;
		return false;
	}
}

bool CSstObject::GlobalDeformSetup(bool use_local_setup)
{
	if (has_global_deformation_)
	{
		def_cross_sections_.clear();
		for (int i = 0; i < cross_sections_.size(); i++)
		{
			if (cross_sections_[i].IsDeformed()) {
				def_cross_sections_.push_back(cross_sections_[i]);
			}
		}
		if (use_local_setup) {
			// locally supported RBF (truncated)
			deformer_->SetUseLocalToSolveGlobal(use_local_setup);
			deformer_->LocalDeformationSetUp(this->def_cross_sections_);
		}
		else {
			// globally supported RBF
			deformer_->GlobalDeformationSetup(this->def_cross_sections_);
		}
		return true;
	}
	else
	{
		std::cerr << "Cannot peform deformation as no information is supplied" << std::endl;
		return false;
	}
}

bool CSstObject::LocalDeformSolve()
{
	// provide right hand side
	// factorize the matrix on the left hand side
	if (has_local_deformation_)
	{
		// do some deformation
		this->def_trimesh_ = this->trimesh_;
		is_deformed_ = true;
		deformer_->LocalDeformationSolve(this->def_cross_sections_, &this->trimesh_, &this->def_trimesh_);
		decoding_mesh(); // have the deformed mesh in the original space
		return true;
	}
	else
	{
		std::cerr << "Cannot peform deformation as no information is supplied" << std::endl;
		return false;
	}
}


bool CSstObject::GlobalDeformSolve()
{
	// provide right hand side
	// factorize the matrix on the left hand side
	if (has_global_deformation_ && !deformer_->GetUseLocalToSolveGlobal())
	{
		// do some deformation
		deformer_->GlobalDeformationSolve(this->def_cross_sections_, &this->trimesh_, &this->def_trimesh_);
		return true;
	}
	else if (has_global_deformation_ && deformer_->GetUseLocalToSolveGlobal())
	{
		deformer_->LocalDeformationSolve(this->def_cross_sections_, &this->trimesh_, &this->def_trimesh_);
		return true;
	}
	else
	{
		std::cerr << "Cannot peform deformation as no information is supplied" << std::endl;
		return false;
	}
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
		ptmp(0) = accu_arc_length[id_min];
		// or ?
		// ptmp(0) += accu_arc_length[id_min];
		trimesh_.data(*viter).set_emb_coord(COpenMeshT::Point(ptmp(0), ptmp(1), ptmp(2)));
	}
	is_encoded_ = true;
}

void CSstObject::extracting_cross_sections(std::vector<int> sid_list)
{
	//std::vector<int> sid_list;
	for (int i = 0; i < sid_list.size(); i++) {
		int sid = sid_list[i];
		COpenMeshT::Point cemb(skeleton_.GetAccumArcLength()[sid], 0.0, 0.0);
		std::vector<COpenMeshT::Point> cs_pts;
		extracting_single_cross_section(cemb, COpenMeshT::Point(1.0, 0.0, 0.0), cs_pts);
		// init cs
		CCrossSection cs;
		cs.SetSid(sid);
		cs.SetEmbProfPts(cs_pts);
		cs.SetClosed(false);
		cs.SetDeformed(false);
		//cs.SetProfPts();
		cross_sections_.push_back(cs);
	}
}

void CSstObject::extracting_single_cross_section(COpenMeshT::Point center, 
	COpenMeshT::Point normal, std::vector<COpenMeshT::Point>& cs_pts)
{
	cs_pts.clear();

	// circulating the mesh to find edges with two endpoints having difference sign to the plane
	int cnt = 0;
	for (auto eiter = trimesh_.edges_begin(); eiter != trimesh_.edges_end(); ++eiter)
	{
		COpenMeshT::HalfedgeHandle heh = trimesh_.halfedge_handle(*eiter, 0);
		COpenMeshT::VertexHandle vh_from = trimesh_.from_vertex_handle(heh);
		COpenMeshT::VertexHandle vh_to = trimesh_.to_vertex_handle(heh);

		COpenMeshT::Point p_from, p_to;
		if (is_encoded_) {
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
	// sorting the vertices
	CGeoCalculator::pts_sorting_alg(cs_pts);
	CGeoCalculator::simplify_polygon(cs_pts);
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
		COpenMeshT::Point ptmp = def_trimesh_.data(*viter).get_emb_coord();
		Vector3d p(ptmp[0], ptmp[1], ptmp[2]);
		Vector3d s(skeletal_pts[sid][0], skeletal_pts[sid][1], skeletal_pts[sid][2]);
		Vector3d v;
		v =  rmf_list[sid].transpose()*p - s;
		def_trimesh_.set_point(*viter, COpenMeshT::Point(v(0), v(1), v(2)));
	}
}

void CSstObject::cross_section_transformed_due_to_skeleton_change()
{
	std::cout << "no cross-sections are deformed according to the skeleton's change" << std::endl;

	
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

	// print def
	if (is_deformed_) {
		if (!PrintSkeleton(fname_dst_skeleton_, def_skeleton_))
			return false;
		if (!PrintCrossSectionList(fname_dst_cross_sections_, def_cross_sections_))
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
	fout_file << "##Skeleton" << std::endl;
	fout_file << "#Points" << std::endl;
	for (int i = 0; i < skel_pts.size(); i++) {
		fout_file << skel_pts[i] << ", " << alength[i] <<  std::endl;
	}
	fout_file << "#RMF" << std::endl;
	for (int i = 0; i < rmf_list.size(); i++) {
		fout_file << rmf_list[i].row(0) << ", "
			<< rmf_list[i].row(1) << ", "
			<< rmf_list[i].row(2) << ", "
			<< std::endl;
	}
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
