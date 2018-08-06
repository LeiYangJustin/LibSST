#ifndef C_MESH_DEFORMATION_H
#define C_MESH_DEFORMATION_H

#include "def_alg_prereq.h"
#include "def_alg_typedef.h"
//#include "../DataColle/mesh_object.h"
//#include "../GeneralTools/geo_calculator.h"


class DEF_ALGCOLLE_CLASS CMeshDeformation
{
public:
	CMeshDeformation();
	CMeshDeformation(double scale);
	~CMeshDeformation();

	// This pre-computes the QR decomp and store it with a scaling metric of diag(0.1, 1, 1)
	void SetUp(const DenseMatrixXd & pts_handles);

	// ONLY applicable to solving the problem given by SetUp function
	void Solve(
		DenseMatrixXd & vecs_left,
		DenseMatrixXd & vecs_right,
		const DenseMatrixXd & roi_pts,
		DenseMatrixXd & roi_vecs);

	void Solve(
		const DenseMatrixXd & vecs_handles,
		const DenseMatrixXd & roi_pts,
		DenseMatrixXd & roi_vecs);

	// using an anisotropic metric but need to compute QR every time
	void CMeshDeformation::DeformWithRBF(
		const DenseMatrixXd & pts_handles_left, const DenseMatrixXd & pts_handles_right,
		const DenseMatrixXd & vecs_handles_left, const DenseMatrixXd & vecs_handles_right,
		const DenseMatrixXd & roi_pts, DenseMatrixXd & roi_vecs);

private:
	//DenseMatrixXd pts_handles_;
	Eigen::ColPivHouseholderQR<DenseMatrixXd> KMQR_; // high accuracy but low efficiency
	Mat3d S_;
	DenseMatrixXd pts_handles_;

private:
	void computeAnisotropicMetric(
		const DenseMatrixXd & pts_at_cs_left,
		const DenseMatrixXd & vec_at_cs_left,
		const DenseMatrixXd & pts_at_cs_right,
		const DenseMatrixXd & vec_at_cs_right,
		Mat3d & metrics);

	// Compute weights
	void computeKernelMatrixRow(VectorXd & v,
		const Point3d & p,
		const DenseMatrixXd & handle_sites, 
		const Mat3d & M = Mat3d::Identity());

	// tools
	void givensTransform(const double x, const double y, Mat2d & m);

};

#endif // !C_MESH_DEFORMATION_H




