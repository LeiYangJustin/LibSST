#include "sst_deformer.h"
#include "../GeneralTools/geo_calculator.h"
//#include "../GeneralTools/visibility_graph.h"
//#include "../GeneralTools/sorting_2d_pts.h"

CSSTDeformer::CSSTDeformer()
	: use_local_to_solve_global(false), deform_type_(DeformType::NotGiven)
{
}

CSSTDeformer::~CSSTDeformer()
{
	for (auto map_iter = p_deformer_map_.begin(); map_iter != p_deformer_map_.end(); ++map_iter)
	{
		delete map_iter->second;
	}
}

void CSSTDeformer::LocalDeformationSetUp(
	const std::vector<CCrossSection> cs_list)
{
	deform_type_ = DeformType::Local;
	for (int i = 0; i < cs_list.size()-1; i++)
	{
		// in
		std::vector<COpenMeshT::Point> left_cs_pts = cs_list[i].GetEmbProfPts();
		std::vector<COpenMeshT::Point> right_cs_pts = cs_list[i+1].GetEmbProfPts();
		// out
		std::pair<int, int> cs_pair(cs_list[i].GetSid(), cs_list[i + 1].GetSid());
		// new deformer
		if (p_deformer_map_.find(cs_pair) == p_deformer_map_.end()) {
			CMeshDeformation* p_deformer = new CMeshDeformation;
			// compute
			get_deformer(left_cs_pts, right_cs_pts, p_deformer);
			p_deformer_map_[cs_pair] = p_deformer;
		}
		//else {
		//	std::cout << "we have the deformer, so we don't need to set up a new one" << std::endl;
		//}
	}
}

bool CSSTDeformer::LocalDeformationSolve(
	const std::vector<CCrossSection> def_cs_list,
	COpenMeshT *p_mesh,
	COpenMeshT *p_mesh_def)
{
	std::map<int, std::pair<COpenMeshT::Point, int>> vh_disp_cnt_map;

	// compute
	//for (int i = 0; i < cs_pair_list_.size(); i++)
	for (int i = 0; i < def_cs_list.size()-1; i++)
	{
		// prepare
		std::pair<int, int> cs_pair(def_cs_list[i].GetSid(), def_cs_list[i + 1].GetSid());

		std::vector<COpenMeshT::Point> left_cs_pts = def_cs_list[i].GetEmbProfPts();
		std::vector<COpenMeshT::Point> left_def_cs_pts = def_cs_list[i].GetDefEmbProfPts();

		std::vector<COpenMeshT::Point> right_cs_pts = def_cs_list[i+1].GetEmbProfPts();
		std::vector<COpenMeshT::Point> right_def_cs_pts = def_cs_list[i+1].GetDefEmbProfPts();

		std::pair<double, double> bound_xCoords;
		bound_xCoords.first = left_def_cs_pts[0][0];
		bound_xCoords.second = right_def_cs_pts[0][0];

		std::cout << cs_pair.first << " " << cs_pair.second << std::endl;
		std::cout << bound_xCoords.first << " " << bound_xCoords.second << std::endl;

		// check
		if (left_cs_pts.size() != left_def_cs_pts.size()) {
			std::cerr << "sizes of the original and deformed left CSs do not match" << std::endl;
			return false;
		}
		if (right_cs_pts.size() != right_def_cs_pts.size()) {
			std::cerr << "sizes of the original and deformed right CSs do not match" << std::endl;
			return false;
		}

		// in
		std::vector<COpenMeshT::Point> left_handles_vecs, right_handles_vecs;
		for (int j = 0; j < left_cs_pts.size(); j++) {
			left_handles_vecs.push_back(left_def_cs_pts[j] - left_cs_pts[j]);
		}
		for (int j = 0; j < right_cs_pts.size(); j++) {
			right_handles_vecs.push_back(right_def_cs_pts[j] - right_cs_pts[j]);
		}

		CMeshDeformation* p_deformer = p_deformer_map_[cs_pair];

		// how to get this?
		// get roi_pts
		std::vector<int> roi_ids;
		std::vector<COpenMeshT::Point> roi_pts;
		for (auto viter = p_mesh->vertices_begin(); viter != p_mesh->vertices_end(); ++viter)
		{
			COpenMeshT::Point pt = p_mesh->data(*viter).get_emb_coord();
			if (pt[0] > bound_xCoords.first && pt[0] < bound_xCoords.second) {
				roi_ids.push_back(viter->idx());
				roi_pts.push_back(pt);
			}
		}
		// out
		std::vector<COpenMeshT::Point> roi_vecs;
		local_deformation_solve(left_handles_vecs, right_handles_vecs, roi_pts, p_deformer, roi_vecs);
		for (int j = 0; j < roi_ids.size(); j++)
		{
			COpenMeshT::Point pt = roi_vecs[j] + roi_pts[j];
			COpenMeshT::VertexHandle def_vh = p_mesh_def->vertex_handle(roi_ids[j]);
			p_mesh_def->data(def_vh).set_emb_coord(pt);
		}
	}

	return true;
}

void CSSTDeformer::GlobalDeformationSetup(std::vector<CCrossSection> cs_list)
{
	std::pair<int, int> cs_pair(0, cs_list.size() - 1);
	// renew deformer if it is not already computed
	if (p_deformer_map_.find(cs_pair) == p_deformer_map_.end()) {
		CMeshDeformation* p_deformer = new CMeshDeformation;
		std::vector<COpenMeshT::Point> all_cs_pts;
		for (int i = 0; i < cs_list.size(); i++) {
			std::vector<COpenMeshT::Point> cs_pts = cs_list[i].GetProfPts();
			all_cs_pts.insert(all_cs_pts.end(), cs_pts.begin(), cs_pts.end());
		}
		get_deformer(all_cs_pts, p_deformer);
		p_deformer_map_[cs_pair] = p_deformer;
	}
}

bool CSSTDeformer::GlobalDeformationSolve(
	const std::vector<CCrossSection> def_cs_list,
	COpenMeshT *p_mesh,
	COpenMeshT *p_mesh_def)
{
	// prepare
	std::vector<COpenMeshT::Point> all_cs_pts, all_def_cs_pts;
	for (int i = 0; i < def_cs_list.size(); i++) {
		std::vector<COpenMeshT::Point> cs_pts = def_cs_list[i].GetProfPts();
		std::vector<COpenMeshT::Point> def_cs_pts = def_cs_list[i].GetDefEmbProfPts();

		all_cs_pts.insert(all_cs_pts.end(), cs_pts.begin(), cs_pts.end());
		all_def_cs_pts.insert(all_def_cs_pts.end(), def_cs_pts.begin(), def_cs_pts.end());
	}
	// check
	if (all_cs_pts.size() != all_def_cs_pts.size()) {
		std::cerr << "sizes of the original and deformed CSs do not match" << std::endl;
		return false;
	}
	// in
	std::vector<COpenMeshT::Point> all_handles_vecs;
	for (int j = 0; j < all_cs_pts.size(); j++) {
		all_handles_vecs.push_back(all_def_cs_pts[j] - all_cs_pts[j]);
		//std::cout << all_def_cs_pts[j] - all_cs_pts[j] << std::endl;
	}
	CMeshDeformation* p_deformer = p_deformer_map_.begin()->second;
	// how to get this?
	std::vector<int> roi_ids;
	std::vector<COpenMeshT::Point> roi_pts;
	for (auto viter = p_mesh->vertices_begin(); viter != p_mesh->vertices_end(); ++viter)
	{
		COpenMeshT::Point pt = p_mesh->point(*viter);
		roi_ids.push_back(viter->idx());
		roi_pts.push_back(pt);
	}
	
	// solve
	std::vector<COpenMeshT::Point> roi_vecs;
	global_deformation_solve(all_handles_vecs, roi_pts, p_deformer, roi_vecs);

	// out
	for (int j = 0; j < roi_ids.size(); j++)
	{
		//std::cout << j << std::endl;
		COpenMeshT::Point pt = roi_vecs[j] + roi_pts[j];
		//std::cout << roi_vecs[j].norm() << std::endl;
		COpenMeshT::VertexHandle def_vh = p_mesh_def->vertex_handle(roi_ids[j]);
		p_mesh_def->set_point(def_vh, pt);
	}

	return true;
}

void CSSTDeformer::get_deformer(
	const std::vector<COpenMeshT::Point> left_cs_pts, 
	const std::vector<COpenMeshT::Point> right_cs_pts,
	CMeshDeformation * deformer)
{
	std::vector<double> handlepts;
	for (int j = 0; j < left_cs_pts.size(); j++)
	{
		handlepts.push_back(left_cs_pts[j][0]);
		handlepts.push_back(left_cs_pts[j][1]);
		handlepts.push_back(left_cs_pts[j][2]);
	}
	for (int j = 0; j < right_cs_pts.size(); j++)
	{
		handlepts.push_back(right_cs_pts[j][0]);
		handlepts.push_back(right_cs_pts[j][1]);
		handlepts.push_back(right_cs_pts[j][2]);
	}
	DenseMatrixXd pts_handles;
	pts_handles = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>(handlepts.data(), 3, handlepts.size() / 3);
	pts_handles.transposeInPlace();
	deformer->SetUp(pts_handles);
}

void CSSTDeformer::get_deformer(const std::vector<COpenMeshT::Point> all_cs_pts, CMeshDeformation * deformer)
{
	std::vector<double> handlepts;
	for (int j = 0; j < all_cs_pts.size(); j++)
	{
		handlepts.push_back(all_cs_pts[j][0]);
		handlepts.push_back(all_cs_pts[j][1]);
		handlepts.push_back(all_cs_pts[j][2]);
	}
	DenseMatrixXd pts_handles;
	pts_handles = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>(handlepts.data(), 3, handlepts.size() / 3);
	pts_handles.transposeInPlace();

	deformer->SetUp(pts_handles);

}


void CSSTDeformer::local_deformation_solve(
	const std::vector<COpenMeshT::Point> left_handles_vecs,
	const std::vector<COpenMeshT::Point> right_handles_vecs,
	const std::vector<COpenMeshT::Point> roi_pts,
	CMeshDeformation* p_deformer,
	std::vector<COpenMeshT::Point> &roi_vecs)
{
	// convert std::vector<Point> to vector<double>
	std::vector<double> left_handles_vecs_data, right_handles_vecs_data, roi_pts_data;
	for (int i = 0; i < left_handles_vecs.size(); i++) {
		left_handles_vecs_data.push_back(left_handles_vecs[i][0]);
		left_handles_vecs_data.push_back(left_handles_vecs[i][1]);
		left_handles_vecs_data.push_back(left_handles_vecs[i][2]);
	}
	for (int i = 0; i < right_handles_vecs.size(); i++) {
		right_handles_vecs_data.push_back(right_handles_vecs[i][0]);
		right_handles_vecs_data.push_back(right_handles_vecs[i][1]);
		right_handles_vecs_data.push_back(right_handles_vecs[i][2]);
	}
	for (int i = 0; i < roi_pts.size(); i++) {
		roi_pts_data.push_back(roi_pts[i][0]);
		roi_pts_data.push_back(roi_pts[i][1]);
		roi_pts_data.push_back(roi_pts[i][2]);
	}

	// initialize the boundary conditions; eigen
	DenseMatrixXd vecs_h_left, vecs_h_right, pts_roi;
	vecs_h_left = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>
		(left_handles_vecs_data.data(), 3, left_handles_vecs_data.size() / 3);
	vecs_h_left.transposeInPlace();
	vecs_h_right = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>
		(right_handles_vecs_data.data(), 3, right_handles_vecs_data.size() / 3);
	vecs_h_right.transposeInPlace();
	pts_roi = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>
		(roi_pts_data.data(), 3, roi_pts_data.size()/3);
	pts_roi.transposeInPlace();

	// solve
	DenseMatrixXd vecs_roi;
	p_deformer->Solve(vecs_h_left, vecs_h_right, pts_roi, vecs_roi);

	roi_vecs.clear();
	for (int i = 0; i < vecs_roi.rows(); i++) {
		Eigen::RowVector3d rv = vecs_roi.row(i);
		roi_vecs.push_back(COpenMeshT::Point(rv(0), rv(1), rv(2)));
	}
}

void CSSTDeformer::global_deformation_solve(
	const std::vector<COpenMeshT::Point> all_handles_vecs, 
	const std::vector<COpenMeshT::Point> roi_pts, 
	CMeshDeformation * p_deformer, 
	std::vector<COpenMeshT::Point> &roi_vecs)
{
	// convert std::vector<Point> to vector<double>
	std::vector<double> all_handles_vecs_data, roi_pts_data;
	for (int i = 0; i < all_handles_vecs.size(); i++) {
		all_handles_vecs_data.push_back(all_handles_vecs[i][0]);
		all_handles_vecs_data.push_back(all_handles_vecs[i][1]);
		all_handles_vecs_data.push_back(all_handles_vecs[i][2]);
	}
	for (int i = 0; i < roi_pts.size(); i++) {
		roi_pts_data.push_back(roi_pts[i][0]);
		roi_pts_data.push_back(roi_pts[i][1]);
		roi_pts_data.push_back(roi_pts[i][2]);
	}

	// initialize the boundary conditions; eigen
	DenseMatrixXd vecs_handles, pts_roi;
	vecs_handles = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>
		(all_handles_vecs_data.data(), 3, all_handles_vecs_data.size() / 3);
	vecs_handles.transposeInPlace();
	pts_roi = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>
		(roi_pts_data.data(), 3, roi_pts_data.size() / 3);
	pts_roi.transposeInPlace();

	// solve
	DenseMatrixXd vecs_roi;
	p_deformer->Solve(vecs_handles, pts_roi, vecs_roi);

	roi_vecs.clear();
	for (int i = 0; i < vecs_roi.rows(); i++) {
		Eigen::RowVector3d rv = vecs_roi.row(i);
		roi_vecs.push_back(COpenMeshT::Point(rv(0), rv(1), rv(2)));
	}
}
