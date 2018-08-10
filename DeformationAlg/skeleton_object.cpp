#include "skeleton_object.h"
#include "../GeneralTools/geo_calculator.h"


void CSkeleton::compute_accumulated_arc_length()
{
	double alength = 0;
	accum_arclength_.push_back(alength);

	for (int i = 1; i < skeletal_pts_.size(); i++) {
		alength += (skeletal_pts_[i] - skeletal_pts_[i - 1]).norm();
		accum_arclength_.push_back(alength);
	}
}

void CSkeleton::compute_rotation_minimizing_frames()
{
	std::vector<double> skelCoordStride3;
	for (int i = 0; i < skeletal_pts_.size(); i++) {
		skelCoordStride3.push_back(skeletal_pts_[i][0]);
		skelCoordStride3.push_back(skeletal_pts_[i][1]);
		skelCoordStride3.push_back(skeletal_pts_[i][2]);
	}

	int numPts = skeletal_pts_.size();
	DenseMatrixXd skel_pts;
	skel_pts = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>(skelCoordStride3.data(), 3, numPts);
	skel_pts.transposeInPlace();

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

	Vector3d Tinit = T.row(0);
	Mat3d Rm, Rtmp;
	Rm.setIdentity();
	Rtmp.setIdentity();
	for (int i = 1; i < 3; i++) {
		std::vector<double> g;
		CGeoCalculator::givensTransform(Tinit(0), Tinit(i), g);
		Rtmp(0, 0) = g[0]; // c
		Rtmp(0, i) = g[1]; // s
		Rtmp(i, 0) = -g[1]; // -s
		Rtmp(i, i) = g[0];  // c
		Rm *= Rtmp;
		//std::cout << "Rtmp = \n" << Rtmp << std::endl;
	}

	DenseMatrixXd R(numPts, 3);
	DenseMatrixXd S(numPts, 3);
	R.setZero();
	S.setZero();
	R.row(0) = Rm.row(1);
	S.row(0) = Rm.row(2);

	//std::cerr << "Initialization of RMF computation: " << std::endl;
	//std::cerr << T.row(0) << std::endl;
	//std::cerr << R.row(0) << std::endl;
	//std::cerr << S.row(0) << std::endl;


	// compute
	double_reflection(skel_pts, T, R, S);
	RMF_list_.clear();
	for (int i = 0; i < numPts; i++)
	{
		// row OR col?
		Mat3d rmf;
		rmf.row(0) = T.row(i);
		rmf.row(1) = R.row(i);
		rmf.row(2) = S.row(i);
		RMF_list_.push_back(rmf);
	}

	////
	//std::ofstream file_out_rmf("rmf_list.txt");
	//if (file_out_rmf.is_open())
	//{
	//	for (int i = 0; i < RMF_list_.size(); i++)
	//	{
	//		file_out_rmf << RMF_list_[i] << std::endl;
	//	}
	//}
	//file_out_rmf.close();
}

void CSkeleton::double_reflection(
	const DenseMatrixXd & skeletal_pts,
	const DenseMatrixXd & T,
	DenseMatrixXd & R,
	DenseMatrixXd & S)
{
	int numPts = skeletal_pts.rows();
	double th = 0.0000001;

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
			Vector3d v1 = skeletal_pts.row(i + 1) - skeletal_pts.row(i);
			Vector3d riL = ri - (2 / v1.dot(v1))*(v1.dot(ri))*v1;
			Vector3d tiL = ti - (2 / v1.dot(v1))*(v1.dot(ti))*v1;
			Vector3d v2 = tii - tiL;
			Vector3d rii = riL - (2 / v2.dot(v2))*(v2.dot(riL))*v2;
			R.row(i + 1) = rii;
			S.row(i + 1) = tii.cross(rii);
		}
	}
}