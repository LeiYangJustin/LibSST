#include "sst_parameterization.h"

#include "../DataColle/data_io.h"
#include "../GeneralTools/sorting_2d_pts.h"

CSSTparameterization::CSSTparameterization()
{
	is_encoded_ = false;
}

//CSSTparameterization::CSSTparameterization(CSSTparameterization * sst)
//{
//	this->triMesh_ = sst->triMesh_;
//
//}


CSSTparameterization::~CSSTparameterization()
{
}

void CSSTparameterization::SetUpAndEncode(COpenMeshT & mesh, std::vector<COpenMeshT::Point> skel_pts)
{
	// Set skeletal points
	std::vector<double> skelCoordStride3;
	for (int i = 0; i < skel_pts.size(); i++) {
		skelCoordStride3.push_back(skel_pts[i][0]);
		skelCoordStride3.push_back(skel_pts[i][1]);
		skelCoordStride3.push_back(skel_pts[i][2]);
	}
	skel_pts_ = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>(skelCoordStride3.data(), 3, skelCoordStride3.size() / 3);
	skel_pts_.transposeInPlace();

	// Set mesh
	triMesh_ = mesh;
	std::vector<double> meshCoordStride3;
	int cnt = 0;
	for (auto viter = triMesh_.vertices_begin(); viter != triMesh_.vertices_end(); ++viter) {
		COpenMeshT::Point p = triMesh_.point(*viter);
		meshCoordStride3.push_back(p[0]);
		meshCoordStride3.push_back(p[1]);
		meshCoordStride3.push_back(p[2]);
		vid_pid_map_.insert(std::make_pair(viter->idx(), cnt));
		pid_vid_map_.insert(std::make_pair(cnt, viter->idx()));
		cnt++;
	}

	ori_pts_ = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>(meshCoordStride3.data(), 3, meshCoordStride3.size() / 3);
	ori_pts_.transposeInPlace();

	Encode();
}

void CSSTparameterization::GetDeformedCrossSectionViaSkelDeformation(
	std::vector<DenseMatrixXd>& def_vecs_list, 
	std::vector<std::vector<COpenMeshT::Point>> & def_cs_pts_list,
	const std::vector<std::vector<COpenMeshT::Point>> & emb_cs_pts_list,
	const DenseMatrixXd & def_skel_pts, 
	const std::vector<int> cs_id_list)
{
	def_vecs_list.clear();
	def_cs_pts_list.clear();

	// encoding or no?
	if (!is_encoded_) {
		std::cerr << "no encoding data!" << std::endl;
		return;
	}

	if (cs_id_list.size() != emb_cs_pts_list.size())
	{
		std::cerr << "inconsistent dimension between cs_id list and cs_pts list!" << std::endl;
		return;
	}
	
	std::vector<Mat3d> rmf_list;
	computeSkeletalRMF(def_skel_pts, rmf_list);
	
	for (int k = 0; k < cs_id_list.size(); k++)
	{
		DenseMatrixXd emb_cs_pts, ori_cs_pts, def_cs_pts, def_cs_vecs;
		std::vector<double> d_cs_pts;
		for (int j = 0; j < emb_cs_pts_list[k].size(); j++)
		{
			d_cs_pts.push_back(emb_cs_pts_list[k][j][0]);
			d_cs_pts.push_back(emb_cs_pts_list[k][j][1]);
			d_cs_pts.push_back(emb_cs_pts_list[k][j][2]);
		}
		emb_cs_pts = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>(d_cs_pts.data(), 3, d_cs_pts.size() / 3);
		emb_cs_pts.transposeInPlace();

		//std::cout << d_cs_pts.size() << std::endl;
		//std::cout << emb_cs_pts_list[k].size() << std::endl;
		//std::cout << emb_cs_pts.rows() << " " << emb_cs_pts.cols() << std::endl;
		//std::cout << std::endl;

		ori_cs_pts.resize(emb_cs_pts.rows(), 3);
		def_cs_pts.resize(emb_cs_pts.rows(), 3);
		int csid = cs_id_list[k];
		//std::cout << csid << std::endl;

		for (int i = 0; i < emb_cs_pts.rows(); i++)
		{
			Vector3d p_emb = emb_cs_pts.row(i);
			if (csid > 0) {
				for (int j = 1; j < csid + 1; j++) {
					p_emb(0) -= (skel_pts_.row(j) - skel_pts_.row(j - 1)).norm();
				}
			}
			// original
			ori_cs_pts.row(i) = p_emb.transpose() * rmf_list_[csid];
			ori_cs_pts.row(i) += skel_pts_.row(csid);
			def_cs_pts.row(i) = p_emb.transpose() * rmf_list[csid];
			def_cs_pts.row(i) += def_skel_pts.row(csid);
		}
		def_cs_vecs = def_cs_pts - ori_cs_pts;
		
		def_vecs_list.push_back(def_cs_vecs);

		// output 
		std::vector<COpenMeshT::Point> def_pts;
		for (int i = 0; i < def_cs_pts.rows(); i++)
		{
			def_pts.push_back(COpenMeshT::Point(def_cs_pts(i, 0), def_cs_pts(i, 1), def_cs_pts(i, 2)));
		}
		def_cs_pts_list.push_back(def_pts);
	}
}

//void CSSTparameterization::SetDeformedSkeleton(std::vector<COpenMeshT::Point> def_skel_pt_list, 
//	DenseMatrixXd & def_skel_pts, std::vector<Mat3d> & rmf_list)
//{
//	// Set skeletal points
//	std::vector<double> skelCoordStride3;
//	for (int i = 0; i < def_skel_pt_list.size(); i++) {
//		skelCoordStride3.push_back(def_skel_pt_list[i][0]);
//		skelCoordStride3.push_back(def_skel_pt_list[i][1]);
//		skelCoordStride3.push_back(def_skel_pt_list[i][2]);
//	}
//	def_skel_pts = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>(skelCoordStride3.data(), 3, skelCoordStride3.size() / 3);
//	def_skel_pts.transposeInPlace();
//
//	computeSkeletalRMF(def_skel_pts, rmf_list);
//}

void CSSTparameterization::Encode()
{
	// skeleton-embedding space
	computeSkeletalRMF(skel_pts_, rmf_list_);

	// encoding
	emb_pts_.resize(ori_pts_.rows(), 3);
	emb_pts_.setZero();
	pid_skeid_map_.clear();
	for (int i = 0; i < ori_pts_.rows(); i++)
	{
		double dmin = std::numeric_limits<double>::max();
		int id_min = -1;
		for (int j = 0; j < skel_pts_.rows(); j++)
		{
			double d = (ori_pts_.row(i) - skel_pts_.row(j)).norm();
			if (d < dmin)
			{
				dmin = d;
				id_min = j;
			}
		}
		pid_skeid_map_.insert(std::make_pair(i, id_min));

		// encode
		Vector3d ptmp;
		ptmp = rmf_list_[id_min] * (ori_pts_.row(i) - skel_pts_.row(id_min)).transpose();
		if (id_min != 0) {
			for (int j = 1; j < id_min + 1; j++) {
				ptmp(0) += (skel_pts_.row(j) - skel_pts_.row(j-1)).norm();
			}
		}
		emb_pts_.row(i) = ptmp;
	}

	// embedded mesh
	embMesh_ = triMesh_;
	for (auto viter = embMesh_.vertices_begin(); viter != embMesh_.vertices_end(); ++viter) {
		int pid = vid_pid_map_[viter->idx()];
		embMesh_.set_point(*viter, COpenMeshT::Point(emb_pts_(pid, 0), emb_pts_(pid, 1), emb_pts_(pid, 2)));
	}

	//
	std::ofstream fout_mesh;
	fout_mesh.open("fDataMesh/embed_mesh.txt");
	for (auto viter = embMesh_.vertices_begin(); viter != embMesh_.vertices_end(); ++viter)
		fout_mesh << embMesh_.point(*viter) << std::endl;
	fout_mesh.close();

	is_encoded_ = true;
}

void CSSTparameterization::DecodeVectorField(DenseMatrixXd & U, const DenseMatrixXd & V)
{
	if (!is_encoded_) {
		std::cerr << "no encoding data!" << std::endl;
		return;
	}
	
	U.resize(V.rows(), 3);
	for (int i = 0; i < V.rows(); i++)
	{
		Vector3d v = V.row(i);
		int id_min = pid_skeid_map_[i];
		//U.row(i) = rmf_list_[id_min].transpose()*v;
		U.row(i) = v.transpose()*rmf_list_[id_min];
	}
}

void CSSTparameterization::DecodeVectorField(DenseMatrixXd & U, const DenseMatrixXd & V, std::vector<int> ids)
{
	if (!is_encoded_) {
		std::cerr << "no encoding data!" << std::endl;
		return;
	}

	assert(V.rows() == ids.size());
	U.resize(V.rows(), 3);
	for (int i = 0; i < V.rows(); i++)
	{
		Vector3d v = V.row(i);
		int id_min = pid_skeid_map_[ids[i]];
		U.row(i) = v.transpose()*rmf_list_[id_min];
	}
}

void CSSTparameterization::DecodePositions(DenseMatrixXd & new_pts)
{
	if (!is_encoded_) {
		std::cerr << "no encoding data!" << std::endl;
		return;
	}

	new_pts.resize(emb_pts_.rows(), 3);
	for (int i = 0; i < emb_pts_.rows(); i++)
	{
		Vector3d p_emb = emb_pts_.row(i);
		int id_min = pid_skeid_map_[i];
		if (id_min != 0) {
			for (int j = 1; j < id_min + 1; j++) {
				p_emb(0) -= (skel_pts_.row(j) - skel_pts_.row(j - 1)).norm();
			}
		}
		new_pts.row(i) = p_emb.transpose() * rmf_list_[id_min];
		new_pts.row(i) += skel_pts_.row(id_min);
	}
}
//
//void CSSTparameterization::SkelDeformation(DenseMatrixXd & new_pts, const DenseMatrixXd & def_skel_pts)
//{
//	if (!is_encoded_) {
//		std::cerr << "no encoding data!" << std::endl;
//		return;
//	}
//
//	std::vector<Mat3d> rmf_list;
//	computeSkeletalRMF(def_skel_pts, rmf_list);
//
//	new_pts.resize(emb_pts_.rows(), 3);
//	for (int i = 0; i < emb_pts_.rows(); i++)
//	{
//		Vector3d p_emb = emb_pts_.row(i);
//		int id_min = pid_skeid_map_[i];
//		if (id_min != 0) {
//			for (int j = 1; j < id_min + 1; j++) {
//				p_emb(0) -= (skel_pts_.row(j) - skel_pts_.row(j - 1)).norm();
//			}
//		}
//		new_pts.row(i) = p_emb.transpose() * rmf_list[id_min];
//		new_pts.row(i) += def_skel_pts.row(id_min);
//	}
//}
//


void CSSTparameterization::DecodePosition(Vector3d & new_pt, 
	const Vector3d & old_pt, const int id_min)
{
	if (!is_encoded_) {
		std::cerr << "no encoding data!" << std::endl;
		return;
	}

	Vector3d tmp_pt = old_pt;
	if (id_min != 0) {
		for (int j = 1; j < id_min + 1; j++) {
			tmp_pt(0) -= (skel_pts_.row(j) - skel_pts_.row(j - 1)).norm();
		}
	}
	new_pt = tmp_pt.transpose() * rmf_list_[id_min];
	new_pt += skel_pts_.row(id_min);
}

void CSSTparameterization::ExtractSingleCrossSection(std::vector<COpenMeshT::Point> &pt_list,
	const COpenMeshT::Point center, const COpenMeshT::Point normal, bool is_encoded)
{
	COpenMeshT mesh;
	if (is_encoded)
		mesh = embMesh_;
	else
		mesh = triMesh_;

	// circulating the mesh to find edges with two endpoints having difference sign to the plane
	int cnt = 0;
	for (auto eiter = mesh.edges_begin(); eiter != mesh.edges_end(); ++eiter) 
	{
		COpenMeshT::HalfedgeHandle heh = mesh.halfedge_handle(*eiter, 0);
		COpenMeshT::VertexHandle vh_from = mesh.from_vertex_handle(heh);
		COpenMeshT::VertexHandle vh_to = mesh.to_vertex_handle(heh);

		COpenMeshT::Point p_from = mesh.point(vh_from);
		COpenMeshT::Point p_to = mesh.point(vh_to);

		double d_from = CGeoCalculator::ComputeDistFromPoint2Plane(p_from, center, normal);
		double d_to = CGeoCalculator::ComputeDistFromPoint2Plane(p_to, center, normal);

		if (d_from*d_to <= 0)
		{
			cnt++;
			double d_deno = fabs(d_from) + fabs(d_to);
			if (d_deno == 0 && 0 == fabs(d_to))
			{
				pt_list.push_back(p_from);
				pt_list.push_back(p_to);
			}
			else {
				double t = fabs(d_from) / d_deno;
				COpenMeshT::Point pt = p_from + t*(p_to - p_from);
				pt_list.push_back(pt);
			}
		}
	}

	// sorting the vertices
	CSorting2DPts sorter(pt_list, 0, 0.1);
}

//void CSSTparameterization::ExtracEmbeddedCrossSectionsWithCScomputer(std::vector<COpenMeshT::Point>& pt_list, COpenMeshT::Point center, COpenMeshT::Point normal)
//{
//	cs_computer_->GetSingleCrossSectionWithFixedNormal(center, normal, pt_list);
//}


void CSSTparameterization::GetUnitTangentVectors(std::vector<COpenMeshT::Point> &tang_vecs)
{
	assert(rmf_list_.size() == skel_pts_.rows());

	for (int i = 0; i < rmf_list_.size(); i++)
	{
		COpenMeshT::Point t(rmf_list_[i](0, 0), rmf_list_[i](0, 1), rmf_list_[i](0, 2));
		tang_vecs.push_back(t);
	}
}

void CSSTparameterization::computeSkeletalRMF(const DenseMatrixXd & skel_pts, std::vector<Mat3d> &rmf_list)
{
	int numPts = skel_pts.rows();

	// init T
	DenseMatrixXd T(numPts, 3);
	T.setZero();
	// first
	Vector3d v = skel_pts.row(1) - skel_pts.row(0);
	v.normalize();
	T.row(0) = v.transpose();
	for (int i = 1; i < numPts; i++)
	{
		Vector3d v = skel_pts.row(i) - skel_pts.row(i - 1);
		v.normalize();
		T.row(i) = v.transpose();
	}

	//// init R, S
	//Mat3d M;
	////M << 1, 0, 0,
	////	0, 0, 1,
	////	0, -1, 0;
	//M << 0, -1, 0,
	//	1, 0, 0,
	//	0, 0, 1;
	//Vector3d Tinit = T.row(0);
	//Vector3d Rinter = M.transpose()*Tinit;
	//Vector3d Sinter = Tinit.cross(Rinter);
	//Sinter.normalize();
	//Vector3d Rinit = Tinit.cross(Sinter);
	//Rinit.normalize();
	//Vector3d Sinit = Tinit.cross(Rinit);
	//Sinit.normalize();
	//DenseMatrixXd R(numPts, 3);
	//DenseMatrixXd S(numPts, 3);
	//R.setZero();
	//S.setZero();
	//R.row(0) = Rinit;
	//S.row(0) = Sinit;

	Vector3d Tinit = T.row(0);
	Mat3d Rm, Rtmp;
	Rm.setIdentity();
	Rtmp.setIdentity();
	for (int i = 1; i < 3; i++) {
		std::vector<double> g;
		CGeoCalculator::givensTransform(Tinit(0), Tinit(i), g);
		Rtmp(1, 1) = g[0]; // c
		Rtmp(1, i) = g[1]; // s
		Rtmp(i, 1) = -g[1]; // -s
		Rtmp(i, i) = g[0];
		Rm *= Rtmp;
	}

	DenseMatrixXd R(numPts, 3);
	DenseMatrixXd S(numPts, 3);
	R.setZero();
	S.setZero();
	R.row(0) = Rm.row(1);
	S.row(0) = Rm.row(2);

	// compute
	doubleReflection(skel_pts, T, R, S);
	rmf_list.clear();
	for (int i = 0; i < numPts; i++)
	{
		// row OR col?
		Mat3d rmf;

		rmf.row(0) = T.row(i);
		rmf.row(1) = R.row(i);
		rmf.row(2) = S.row(i);
		rmf_list.push_back(rmf);
	}
}

void CSSTparameterization::doubleReflection(
	const DenseMatrixXd & skel_pts,
	const DenseMatrixXd & T,
	DenseMatrixXd & R,
	DenseMatrixXd & S)
{
	int numPts = skel_pts.rows();
	double th = 0.001;

	for (int i = 0; i < numPts - 1; i++)
	{
		Vector3d ri = R.row(i);
		Vector3d ti = T.row(i);
		Vector3d tii = T.row(i + 1);

		if (tii.norm() < th)
		{
			R.row(i + 1) = R.row(i);
			S.row(i + 1) = S.row(i);
		}
		else {
			Vector3d v1 = skel_pts.row(i + 1) - skel_pts.row(i);
			Vector3d riL = ri - (2 / v1.dot(v1))*(v1.dot(ri))*v1;
			Vector3d tiL = ti - (2 / v1.dot(v1))*(v1.dot(ti))*v1;
			Vector3d v2 = tii - tiL;
			Vector3d rii = riL - (2 / v2.dot(v2))*(v2.dot(riL))*v2;
			R.row(i + 1) = rii;
			S.row(i + 1) = tii.cross(rii);
		}
	}
}

