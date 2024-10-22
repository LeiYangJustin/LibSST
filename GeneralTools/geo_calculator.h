#ifndef C_GEO_CALCULATOR_H
#define C_GEO_CALCULATOR_H

#include "general_tool_prereq.h"
#include "../DataColle/mesh_object.h"
#include <Eigen/Core>

class GENERAL_TOOLS_CLASS CGeoCalculator
{
public:
	static double ComputeDistFromPoint2Plane(COpenMeshT::Point p, COpenMeshT::Point b, COpenMeshT::Point n);
	static double ComputeLinearity(std::vector<COpenMeshT::Point> p_list);
	static double ComputeLineLinearity(std::vector<COpenMeshT::Point> p_list, COpenMeshT::Point ctr, COpenMeshT::Point dir);
	
	static void FittingLineSegmentToPointSet(std::vector<COpenMeshT::Point> p_list, 
		COpenMeshT::Point &pstart, COpenMeshT::Point &pend);
	static void FittingLineSegmentToPointSet(std::vector<COpenMeshT::Point> p_list, 
		COpenMeshT::Point base, COpenMeshT::Point direction, 
		COpenMeshT::Point &pstart, COpenMeshT::Point &pend);
	static void FittingPlaneToPointSet(std::vector<COpenMeshT::Point> p_list, std::vector<double>& transform);

	static void ComputePrincipalDirection(std::vector<COpenMeshT::Point> p_list, std::vector<double>& transform);
	static void ComputePCAwithEIG(std::vector<COpenMeshT::Point> p_list,
		std::vector<std::vector<double>>& eigenvectors, std::vector<double>& eigenvalues);
	static void ComputePCAwithSVD(std::vector<COpenMeshT::Point> p_list,
		std::vector<std::vector<double>>& eigenvectors, std::vector<double>& eigenvalues);

	static void ProjectPointToPlane(COpenMeshT::Point point, std::vector<double> plane, COpenMeshT::Point & projection);
	static double ComputeDistPointToLine(COpenMeshT::Point point, COpenMeshT::Point start_point, COpenMeshT::Point end_point);

	static inline double innerProduct(COpenMeshT::Point p1, COpenMeshT::Point p2) {
		return (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]);
	}
	//static inline double twoPointSqrDistance(COpenMeshT::Point p1, COpenMeshT::Point p2)
	//{
	//	return (p1[0] - p2[0])*(p1[0] - p2[0]) + (p1[1] - p2[1])*(p1[1] - p2[1]) + (p1[2] - p2[2])*(p1[2] - p2[2]);
	//}

	// output cosA and sinA
	static void givensTransform(double x, double y, std::vector<double> &g);

	// CGAL dependent
	static void computeOBB(const std::vector<COpenMeshT::Point> & pts, std::vector<double> & obb);
	
	// Bezier curve discretization
	static std::vector<double> getBernsteinBasis(int n, double step, Eigen::MatrixXd & Bmat);
	static void getBernsteinBasis(int n, std::vector<double> t, Eigen::MatrixXd & Bmat);
	static void getSampleFromBezier(
		Eigen::MatrixXd & sample_pts,
		Eigen::MatrixXd & ctrl_pts,
		double step=0.01);
	static void getEqualArcLengthSampleFromBezier(
		std::vector<COpenMeshT::Point> & sample_pts,
		std::vector<double> & tparas,
		const std::vector<COpenMeshT::Point> & ctrl_pts,
		int numSamples);
	// nchoosek
	static double nchoosek(int n, int k);

	static void convertOMptToEigenMat(const std::vector<COpenMeshT::Point> om_pts, Eigen::MatrixXd & Mat);
	static void convertEigenMatToOMpt(const Eigen::MatrixXd & Mat, std::vector<COpenMeshT::Point> & om_pts);

	// pt sorting alg; 
	static void pts_sorting_alg(std::vector<COpenMeshT::Point> &pts);
	static void simplify_polygon(std::vector<COpenMeshT::Point> &pts, double dT = 0.5);
	static void reconstruct_curve_from_pointset(std::vector<COpenMeshT::Point> &pts, float tolOMT=1.0);
	
	// spacing = a percentage of step size over total length
	static void sample_polygon(std::vector<COpenMeshT::Point> &pts, double spacing, bool is_closed);
	
	static void sample_polygon_parameters(
		std::vector<double> &tparas, const std::vector<COpenMeshT::Point> pts, 
		double spacing, bool is_closed);

	static void sample_polygon_with_parameters(
		std::vector<double> tparas,
		std::vector<COpenMeshT::Point> cs_pts, 
		std::vector<COpenMeshT::Point> &out_pts, 
		bool is_closed);
};

#endif // !C_GEO_CALCULATOR_H



