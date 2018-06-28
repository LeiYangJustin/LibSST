#include "mesh_deformation.h"
#include <chrono>

CMeshDeformation::CMeshDeformation()
{
	S_.setIdentity();
	//S_(0, 0) = 0.01; // scaling
}


CMeshDeformation::~CMeshDeformation()
{
}

void CMeshDeformation::SetUp(const DenseMatrixXd & pts_handles)
{
	pts_handles_ = pts_handles;
	int numHandles = pts_handles_.rows();
	DenseMatrixXd KerMat(numHandles, numHandles);
	KerMat.setZero();
	for (int i = 0; i < numHandles; i++) {
		VectorXd v;
		Point3d p = pts_handles_.row(i);
		computeKernelMatrixRow(v, p, pts_handles_, S_);
		KerMat.row(i) = v;
	}
	KMQR_ = KerMat.colPivHouseholderQr();
	//std::cout << KerMat.rows() << ", " << KerMat.cols() << std::endl;
}

void CMeshDeformation::Solve(
	DenseMatrixXd & vecs_left,
	DenseMatrixXd & vecs_right,
	const DenseMatrixXd & roi_pts,
	DenseMatrixXd & roi_vecs)
{
	DenseMatrixXd vecs_handles;
	if (vecs_left.rows() == 0 && vecs_right.rows() != 0)
	{
		const int rows_left = pts_handles_.rows() - vecs_right.rows();
		vecs_left.resize(rows_left, 3);
		vecs_left.setZero();
		//
		vecs_handles.resize(vecs_left.rows() + vecs_right.rows(), 3);
		vecs_handles << vecs_left,
			vecs_right;
	}
	else if (vecs_left.rows() != 0 && vecs_right.rows() == 0)
	{
		const int rows_right = pts_handles_.rows() - vecs_left.rows();
		vecs_right.resize(rows_right, 3);
		vecs_right.setZero();
		//
		vecs_handles.resize(vecs_left.rows() + vecs_right.rows(), 3);
		vecs_handles << vecs_left,
			vecs_right;
	}
	else if (vecs_left.rows() != 0 && vecs_right.rows() != 0) {
		vecs_handles.resize(vecs_left.rows() + vecs_right.rows(), 3);
		vecs_handles << vecs_left,
			vecs_right;
		assert(vecs_handles.rows() == pts_handles_.rows());
	}
	else {
		// something wrong
	}
	
	DenseMatrixXd w = KMQR_.solve(vecs_handles);
	// obtain the (embedded) deformation fields
	int numROIPts = roi_pts.rows();
	roi_vecs.resize(numROIPts, 3);
	roi_vecs.setZero();
	for (int i = 0; i < numROIPts; i++) {
		Eigen::VectorXd v;
		computeKernelMatrixRow(v, roi_pts.row(i), pts_handles_, S_);
		roi_vecs.row(i) = w.transpose()*v;
	}
}

void CMeshDeformation::Solve(
	const DenseMatrixXd & vecs_handles, 
	const DenseMatrixXd & roi_pts, 
	DenseMatrixXd & roi_vecs)
{
	// SOME CODES
	DenseMatrixXd w = KMQR_.solve(vecs_handles);
	// obtain the (embedded) deformation fields
	int numROIPts = roi_pts.rows();
	roi_vecs.resize(numROIPts, 3);
	roi_vecs.setZero();

	for (int i = 0; i < numROIPts; i++) {
		//std::cout << "Time spent for decomposition Solve: " << std::endl;
		//auto t1 = std::chrono::high_resolution_clock::now();
		Eigen::VectorXd v;
		computeKernelMatrixRow(v, roi_pts.row(i), pts_handles_, S_);
		roi_vecs.row(i) = w.transpose()*v;
		//auto t2 = std::chrono::high_resolution_clock::now();
		//std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
		//std::cout << time_span.count() << " seconds" << std::endl;
	}
}

void CMeshDeformation::DeformWithRBF(
	const DenseMatrixXd & pts_handles_left, const DenseMatrixXd & pts_handles_right, 
	const DenseMatrixXd & vecs_handles_left, const DenseMatrixXd & vecs_handles_right,
	const DenseMatrixXd & roi_pts, DenseMatrixXd & roi_vecs)
{
	// compute anisotropic metric
	Mat3d M;
	M.setIdentity();
	computeAnisotropicMetric(pts_handles_left, vecs_handles_left, pts_handles_right, vecs_handles_right, M);

	// concatenate matrices
	DenseMatrixXd pts_handle(pts_handles_left.rows() + pts_handles_right.rows(), 3);
	pts_handle << pts_handles_left,
				pts_handles_right;

	DenseMatrixXd vecs_handle(vecs_handles_left.rows() + vecs_handles_right.rows(), 3);
	vecs_handle << vecs_handles_left,
		vecs_handles_right;

	int numHandles = pts_handle.rows();
	DenseMatrixXd w;
	DenseMatrixXd KerMat(numHandles, numHandles);
	KerMat.setZero();
	for (int i = 0; i < numHandles; i++) {
		VectorXd v;
		Point3d p = pts_handle.row(i);
		computeKernelMatrixRow(v, p, pts_handle, M);
		KerMat.row(i) = v;
	}
	// solve for the weights with QR decomposition
	w = KerMat.colPivHouseholderQr().solve(vecs_handle);
	// check if solution exists
	double relative_error = (KerMat*w - vecs_handle).norm() / vecs_handle.norm(); // norm() is L2 norm
	std::cout << "The relative error is:\n" << relative_error << std::endl;

	//std::cout << w.transpose() << std::endl;

	// obtain the (embedded) deformation fields
	int numROIPts = roi_pts.rows();
	roi_vecs.resize(numROIPts, 3);
	for (int i = 0; i < numROIPts; i++) {
		Eigen::VectorXd v;
		computeKernelMatrixRow(v, roi_pts.row(i), pts_handle, M);
		roi_vecs.row(i) = w.transpose()*v;
	}
}

void CMeshDeformation::computeAnisotropicMetric(
	const DenseMatrixXd & pts_at_cs_left,	/*Input; Points at one Cross-section*/
	const DenseMatrixXd & vec_at_cs_left,	/*Input; Vectors at points at one Cross-section*/
	const DenseMatrixXd & pts_at_cs_right, /*Input; Points at the other Cross-section*/
	const DenseMatrixXd & vec_at_cs_right, /*Input; Points at the other Cross-section*/
	Mat3d & M /*Output; metric tensor for this segment*/)
{
	// left
	Point3d p_center_left(0, 0, 0);
	double sum_w = 0.0;
	for (int i = 0; i < pts_at_cs_left.rows(); i++)
	{
		p_center_left += vec_at_cs_left.row(i).norm()*pts_at_cs_left.row(i);
		sum_w += vec_at_cs_left.row(i).norm();
	}
	p_center_left /= sum_w;

	// right
	Point3d p_center_right(0, 0, 0);
	sum_w = 0.0;
	for (int i = 0; i < pts_at_cs_right.rows(); i++)
	{
		p_center_left += vec_at_cs_right.row(i).norm()*pts_at_cs_right.row(i);
		sum_w += vec_at_cs_right.row(i).norm();
	}
	p_center_right /= sum_w;

	// direction
	Vector3d vdir = p_center_left - p_center_right;
	vdir.normalize();

	Mat3d Rm, Rtmp; 
	Rm.setIdentity();
	Rtmp.setIdentity();
	for (int i = 1; i < 3; i++) {
		Mat3d g;
		givensTransform(vdir(0), vdir(i), g);
		Rtmp(1, 1) = g(1, 1);
		Rtmp(1, i) = g(1, 2);
		Rtmp(i, 1) = g(2, 1);
		Rtmp(i, i) = g(2, 2);
		Rm *= Rtmp;
	}

	// this should be changed to use the eigenvalue to decide
	Mat3d Sm; Sm.setIdentity();
	double n = 0; 
	double re = vdir(0) / 10;
	while (re > 1){
		re = re / 10;
		n++;
	}
	Sm(0, 0) = pow(10, n);
	M = Sm*Rm;
}

void CMeshDeformation::computeKernelMatrixRow(VectorXd & v,
	const Point3d & p,
	const DenseMatrixXd & handle_sites,
	const Mat3d & M)
{
	double th = 0.001;
	int numHandles = handle_sites.rows();
	v.resize(numHandles);
	v.setZero();

	// accerlarate
	if (M.isIdentity()) {
		for (int i = 0; i < numHandles; i++) {
			// this is not compact
			Point3d c = handle_sites.row(i);
			Vector3d vtmp = p - c;
			v(i) = vtmp.norm()*vtmp.norm()*vtmp.norm();
		}
	}
	else {
		for (int i = 0; i < numHandles; i++) {
			// this is not compact
			Point3d c = handle_sites.row(i);
			Vector3d vtmp = M*(p - c);
			v(i) = vtmp.norm()*vtmp.norm()*vtmp.norm();
		}
	}

}

void CMeshDeformation::givensTransform(const double x, const double y, Mat2d & m)
{
	double absx = fabs(x);
	double c, s;

	if (absx == 0.0)
	{
		c = 0.0; s = 1.0;
	}
	else {
		double nrm = Eigen::Vector2d(x, y).norm();
		c = absx / nrm;
		s = x / absx*(y / nrm);
	}
	m << c, s,
		-s, c;
}

//bool CMeshDeformation::ReadDeformTestExample(std::string fileName)
//{
//	// check open
//	std::ifstream inFile;
//	inFile.open(fileName);
//	if (!inFile) {
//		std::cerr << "Unable to open file " << fileName << std::endl;
//		return false;
//	}
//	else {
//		std::cerr << "Txt Reader " << fileName << std::endl;
//	}
//
//	// read
//	std::vector<double> pts;
//	std::vector<double> hpts;
//	std::vector<double> hvecs;
//	int flag = 1; // 1 for pts, 2 for handles, 3 for handle displacement
//	int cnt = 0;
//	while (!inFile.eof()) {
//		std::string line;
//		getline(inFile, line);
//		//std::cout << "line: "<< line << std::endl;
//		std::stringstream ss(line);
//		if (line.size() == 0)
//			break;
//		else if (line[0] == 'P')
//			flag = 1;
//		else if (line[0] == 'H')
//			flag = 2;
//		else if (line[0] == 'D')
//			flag = 3;
//		else {
//			std::vector<double> values;
//			for (int i = 0; i < 3; i++) {
//				double d;
//				ss >> d;
//				values.push_back(d);
//			}
//			if (flag == 1) {
//				// add to pts
//				pts.insert(pts.end(), values.begin(), values.end());
//			}
//			else if (flag == 2){
//				// add to handles
//				hpts.insert(hpts.end(), values.begin(), values.end());
//			}
//			else if (flag == 3) {
//				// add to displacement
//				hvecs.insert(hvecs.end(), values.begin(), values.end());
//			}
//		}
//	}
//
//	// map vector to matrix
//	tPts_ = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>(pts.data(), 3, pts.size() / 3);
//	tHPts_ = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>(hpts.data(), 3, hpts.size() / 3);
//	tHVecs_ = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>(hvecs.data(), 3, hvecs.size() / 3);
//	tPts_.transposeInPlace();
//	tHPts_.transposeInPlace();
//	tHVecs_.transposeInPlace();
//	return true;
//}
//
//bool CMeshDeformation::ReadRMFTestExample(std::string fileName)
//{
//	// check open
//	std::ifstream inFile;
//	inFile.open(fileName);
//	if (!inFile) {
//		std::cerr << "Unable to open file " << fileName << std::endl;
//		return false;
//	}
//	else {
//		std::cerr << "Txt Reader " << fileName << std::endl;
//	}
//
//	// read
//	std::vector<double> pts;
//	int flag = 1; // 1 for pts, 2 for handles, 3 for handle displacement
//	int cnt = 0;
//	while (!inFile.eof()) {
//		std::string line;
//		getline(inFile, line);
//		//std::cout << "line: "<< line << std::endl;
//		std::stringstream ss(line);
//		if (line.size() == 0)
//			break;
//		for (int i = 0; i < 3; i++) {
//			double d;
//			ss >> d;
//			pts.push_back(d);
//		}
//	}
//
//	// map vector to matrix
//	tPts_ = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>(pts.data(), 3, pts.size() / 3);
//	tPts_.transposeInPlace();
//	return true;
//}

//void CMeshDeformation::TestComputeRMF()
//{
//	ReadRMFTestExample("deform_tests/curvePts.txt");
//
//	std::vector<Mat3d> rmf_list;
//	computeSkeletalRMF(tPts_, rmf_list);
//
//	std::ofstream fout("rmf.m");
//	fout << "RV = [\n";
//	for (int i = 0; i < rmf_list.size(); i++)
//	{
//		fout << rmf_list[i].row(0) << " ";
//		fout << rmf_list[i].row(1) << " ";
//		fout << rmf_list[i].row(2) << "\n";
//	}
//	fout << "];\n";
//	fout.close();
//
//	std::cout << "TestComputeRMF Done" << std::endl;
//}
//
//void CMeshDeformation::TestDeformWithRBF()
//{
//	ReadDeformTestExample("deform_tests/test.txt");
//
//	////
//	//using namespace std::chrono;
//	//high_resolution_clock::time_point t1 = high_resolution_clock::now();
//	////
//
//	Mat3d M;
//	M.setIdentity();
//
//	int numHandles = tHPts_.rows();
//	DenseMatrixXd w;
//	DenseMatrixXd KerMat(numHandles, numHandles);
//	KerMat.setZero();
//	for (int i = 0; i < numHandles; i++) {
//		VectorXd v;
//		Point3d p = tHPts_.row(i);
//		computeKernelMatrixRow(p, tHPts_, M, v);
//		KerMat.row(i) = v;
//	}
//	// solve for the weights with QR decomposition
//	w = KerMat.colPivHouseholderQr().solve(tHVecs_);
//	// check if solution exists
//	double relative_error = (KerMat*w - tHVecs_).norm() / tHVecs_.norm(); // norm() is L2 norm
//	std::cout << "The relative error is " << relative_error << std::endl;
//	
//	//std::cout << w.transpose() << std::endl;
//
//	// obtain the (embedded) deformation fields
//	int numROIPts = tPts_.rows();
//	std::cout << "numROIPts = " << numROIPts << std::endl;
//	DenseMatrixXd roi_vecs;
//	roi_vecs.resize(numROIPts, 3);
//	for (int i = 0; i < numROIPts; i++) {
//		Eigen::VectorXd v;
//		computeKernelMatrixRow(tPts_.row(i), tHPts_, M, v);
//		roi_vecs.row(i) = w.transpose()*v;
//	}
//
//	//std::ofstream fout("roi_vecs.m");
//	//fout << "RV = [\n";
//	//fout << roi_vecs;
//	//fout << "];\n";
//	//fout.close();
//
//	////
//	//high_resolution_clock::time_point t2 = high_resolution_clock::now();
//	//duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
//	//std::cout << time_span.count() << " sec"<< std::endl;
//	////
//
//	std::cout << "TestDeformWithRBF Done" << std::endl;
//}