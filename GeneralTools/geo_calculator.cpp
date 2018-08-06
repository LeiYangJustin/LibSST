#include "geo_calculator.h"
#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <assert.h>

#include "cgal_predef.h"


//// CGAL 
// Compute OBB
#include <CGAL/Polygon_2.h>
#include <CGAL/min_quadrilateral_2.h>
// convex hull
#include <CGAL/convex_hull_2.h>
// polyline simplification
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/Polyline_simplification_2/Squared_distance_cost.h>
// curve reconstruction from unorganized point set
#include <CGAL/Optimal_transportation_reconstruction_2.h>
// random point generator
#include <CGAL/point_generators_2.h>

double CGeoCalculator::ComputeDistFromPoint2Plane(COpenMeshT::Point p, COpenMeshT::Point b, COpenMeshT::Point n)
{
	// assert n is a unit vector
	assert((n.norm() - 1.0) < 10e-5);

	COpenMeshT::Point v = p - b;
	return innerProduct(v, n);
}

double CGeoCalculator::ComputeLinearity(std::vector<COpenMeshT::Point> p_list)
{
	if (p_list.size() < 3)
	{
		return 1.0;
	}
	else
	{
		std::vector<std::vector<double>> eigenvectors;
		std::vector<double> lambda;
		ComputePCAwithSVD(p_list, eigenvectors, lambda);
		//ComputePCAwithEIG(p_list, eigenvectors, lambda);
		return lambda[0] / (lambda[0] + lambda[1] + lambda[2]);
	}
}

double CGeoCalculator::ComputeLineLinearity(std::vector<COpenMeshT::Point> p_list, COpenMeshT::Point ctr, COpenMeshT::Point dir)
{
	double sum_dev = 0.0;
	for (int i = 0; i < p_list.size(); i++)
	{
		double dev = (p_list[i] - ctr - innerProduct(p_list[i] - ctr, dir)*dir).norm();
		sum_dev += dev;
	}
	std::cout << "center = " << ctr << std::endl;
	std::cout << "direction = " << dir << std::endl;
	std::cout << sum_dev <<" "<< p_list.size() << std::endl;
	sum_dev /= double(p_list.size());
	return sum_dev;
}

void CGeoCalculator::FittingLineSegmentToPointSet(std::vector<COpenMeshT::Point> p_list, COpenMeshT::Point & pstart, COpenMeshT::Point & pend)
{
	// base 
	COpenMeshT::Point base(0, 0, 0);
	for (int i = 0; i < p_list.size(); i++)
	{
		base += p_list[i];
	}
	base = base / double(p_list.size());
	for (int i = 0; i < p_list.size(); i++)
	{
		p_list[i] -= base;
	}

	// get line direction
	std::vector<double> direction;
	ComputePrincipalDirection(p_list, direction);

	std::cout << "line direction: ";
	std::copy(direction.begin(), direction.end(), std::ostream_iterator<double>(std::cout, " "));
	std::cout << std::endl;

	// project every point to direction and find the max and min parameters
	double max_t = 0, min_t = 0;
	COpenMeshT::Point v(direction[0], direction[1], direction[2]);
	for (int i = 0; i < p_list.size(); i++) {
		double t = innerProduct(v, p_list[i]);
		std::cout << "t = " << t << std::endl;
		if (t > max_t)
			max_t = t;
		if (t < min_t)
			min_t = t;
	}
	std::cout << "min_t = "<< min_t << ", max_t = " << max_t << std::endl;
	// output
	pstart = base + v*min_t;
	pend = base + v*max_t;
}

void CGeoCalculator::FittingLineSegmentToPointSet(std::vector<COpenMeshT::Point> p_list, 
	COpenMeshT::Point base, COpenMeshT::Point dir, 
	COpenMeshT::Point & pstart, COpenMeshT::Point & pend)
{
	// base 
	for (int i = 0; i < p_list.size(); i++)
	{
		p_list[i] -= base;
	}
	// project every point to direction and find the max and min parameters
	double max_t = 0, min_t = 0;
	for (int i = 0; i < p_list.size(); i++) {
		double t = innerProduct(dir, p_list[i]);
		std::cout << "t = " << t << std::endl;
		if (t > max_t)
			max_t = t;
		if (t < min_t)
			min_t = t;
	}
	std::cout << "min_t = " << min_t << ", max_t = " << max_t << std::endl;
	// output
	pstart = base + dir*min_t;
	pend = base + dir*max_t;
}

void CGeoCalculator::FittingPlaneToPointSet(std::vector<COpenMeshT::Point> p_list, std::vector<double>& transform)
{
	transform.clear();
	int Dim = 3;
	Eigen::MatrixXd P(p_list.size(), Dim);
	Eigen::MatrixXd V(Dim, Dim);

	
	COpenMeshT::Point center_p(0, 0, 0);
	for (int i = 0; i < p_list.size(); i++)
	{
		center_p += p_list[i];
	}
	center_p = center_p / double(p_list.size());
	for (int i = 0; i < p_list.size(); i++)
	{
		P(i, 0) = p_list[i][0] - center_p[0];
		P(i, 1) = p_list[i][1] - center_p[1];
		P(i, 2) = p_list[i][2] - center_p[2];
	}
	

	Eigen::JacobiSVD<Eigen::MatrixXd> svd(P, Eigen::ComputeFullU | Eigen::ComputeFullV);
	V = svd.matrixV();

	Eigen::Vector3d n = V.col(2);
	n.normalized();
	transform.push_back(n[0]);
	transform.push_back(n[1]);
	transform.push_back(n[2]);
}

void CGeoCalculator::ComputePrincipalDirection(std::vector<COpenMeshT::Point> p_list, std::vector<double>& transform)
{
	transform.clear();
	int Dim = 3;
	Eigen::MatrixXd P(p_list.size(), Dim);
	Eigen::MatrixXd V(Dim, Dim);

	COpenMeshT::Point center_p(0, 0, 0);
	for (int i = 0; i < p_list.size(); i++)
	{
		center_p += p_list[i];
	}
	center_p = center_p / double(p_list.size());
	for (int i = 0; i < p_list.size(); i++)
	{
		P(i, 0) = p_list[i][0] - center_p[0];
		P(i, 1) = p_list[i][1] - center_p[1];
		P(i, 2) = p_list[i][2] - center_p[2];
	}

	Eigen::JacobiSVD<Eigen::MatrixXd> svd(P, Eigen::ComputeThinU | Eigen::ComputeThinV);
	V = svd.matrixV();

	//std::cout << "P = [" << P << "]" << std::endl;
	//std::cout << "V = [" << V << "]" << std::endl;

	Eigen::Vector3d n = V.col(0);
	n.normalized();
	transform.push_back(n[0]);
	transform.push_back(n[1]);
	transform.push_back(n[2]);
}

void CGeoCalculator::ComputePCAwithEIG(std::vector<COpenMeshT::Point> p_list, 
	std::vector<std::vector<double>>& eigenvectors, std::vector<double>& eigenvalues)
{
	eigenvectors.clear();
	int Dim = 3;
	Eigen::MatrixXd P(p_list.size(), Dim);
	Eigen::MatrixXd V(Dim, Dim);
	Eigen::MatrixXd S(Dim, Dim);

	// construct
	COpenMeshT::Point center_p(0, 0, 0);
	for (int i = 0; i < p_list.size(); i++)
	{
		center_p += p_list[i];
	}
	center_p = center_p / double(p_list.size());
	for (int i = 0; i < p_list.size(); i++)
	{
		P(i, 0) = p_list[i][0] - center_p[0];
		P(i, 1) = p_list[i][1] - center_p[1];
		P(i, 2) = p_list[i][2] - center_p[2];
	}

	Eigen::MatrixXd Pt = P.transpose();
	Eigen::MatrixXd PtP = P.transpose() * P;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(PtP);
	S = eigen_solver.eigenvalues().asDiagonal();
	V = eigen_solver.eigenvectors();

	//std::cout << "using ComputePCAwithEIG function" << std::endl;
	//std::cout << "S = [" << S << "]" << std::endl;
	//std::cout << "V = [" << V << "]" << std::endl;

	// make it in the descending order
	for (int i = Dim; i > 0; i--)
	{
		std::vector<double> tmp_vector;
		Eigen::Vector3d n = V.col(i-1);
		n.normalized();
		tmp_vector.push_back(n[0]);
		tmp_vector.push_back(n[1]);
		tmp_vector.push_back(n[2]);
		eigenvectors.push_back(tmp_vector);
		eigenvalues.push_back(S(i-1, i-1));
	}
}

void CGeoCalculator::ComputePCAwithSVD(std::vector<COpenMeshT::Point> p_list, 
	std::vector<std::vector<double>>& eigenvectors, std::vector<double>& eigenvalues)
{
	eigenvectors.clear();
	eigenvalues.clear();

	int Dim = 3;
	Eigen::MatrixXd P(p_list.size(), Dim);
	Eigen::MatrixXd V(Dim, Dim);
	Eigen::VectorXd S(Dim);

	COpenMeshT::Point center_p(0, 0, 0);
	for (int i = 0; i < p_list.size(); i++)
	{
		center_p += p_list[i];
	}
	center_p = center_p / double(p_list.size());
	for (int i = 0; i < p_list.size(); i++)
	{
		P(i, 0) = p_list[i][0] - center_p[0];
		P(i, 1) = p_list[i][1] - center_p[1];
		P(i, 2) = p_list[i][2] - center_p[2];
	}
	Eigen::MatrixXd CovP = P.transpose() * P;
	CovP = CovP / double(p_list.size());
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(CovP, Eigen::ComputeThinU | Eigen::ComputeThinV);
	V = svd.matrixV();
	S = svd.singularValues();

	//std::cout << "using ComputePCAwithSVD function" << std::endl;
	//std::cout << "S = [" << S << "]" << std::endl;
	//std::cout << "V = [" << V << "]" << std::endl;

	for (int i = 0; i < 3; i++)
	{
		std::vector<double> tmp_vector;
		Eigen::Vector3d n = V.col(i);
		n.normalized();
		tmp_vector.push_back(n[0]);
		tmp_vector.push_back(n[1]);
		tmp_vector.push_back(n[2]);
		eigenvectors.push_back(tmp_vector);
		double s = sqrt(S(i));
		eigenvalues.push_back(s);
	}
}

void CGeoCalculator::ProjectPointToPlane(COpenMeshT::Point point, std::vector<double> plane, COpenMeshT::Point & projection)
{
	double d = point[0] * plane[0] + point[1] * plane[1] + point[2] * plane[2] + plane[3];
	COpenMeshT::Point v(plane[0] * d, plane[1] * d, plane[2] * d);
	projection = point - v;
}

double CGeoCalculator::ComputeDistPointToLine(COpenMeshT::Point q, COpenMeshT::Point sp, COpenMeshT::Point tp)
{
	COpenMeshT::Point dir = tp - sp;
	dir /= dir.norm();
	return (q - sp - innerProduct(q - sp, dir)*dir).norm();
}

void CGeoCalculator::givensTransform(double x, double y, std::vector<double>& g)
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
	g.push_back(c);
	g.push_back(s);
}

void CGeoCalculator::computeOBB(const std::vector<COpenMeshT::Point> & pts, std::vector<double> & out_obb)
{
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef K::Point_2								Point_2;
	typedef CGAL::Polygon_2<K>						Polygon_2;

	std::vector<Point_2> pointset;
	for (int i = 0; i < pts.size(); i++)
	{
		pointset.push_back(Point_2(pts[i][1], pts[i][2]));
	}

	Polygon_2 ch, obb;
	CGAL::convex_hull_2(pointset.begin(), pointset.end(), std::back_inserter(ch));
	CGAL::min_rectangle_2(ch.vertices_begin(), ch.vertices_end(), std::back_inserter(obb));
	out_obb.clear();
	out_obb.push_back(obb.bbox().xmin());
	out_obb.push_back(obb.bbox().xmax());
	out_obb.push_back(obb.bbox().ymin());
	out_obb.push_back(obb.bbox().xmax());
}

std::vector<double>  CGeoCalculator::getBernsteinBasis(int n, double step, Eigen::MatrixXd & Bmat)
{
	// parameters
	std::vector<double> t;
	double it = 0.0;
	t.push_back(it);
	while (it < 1)
	{
		it += step;
		t.push_back(it);
	}
	if (fabs(t.back() - 1.0) > 0.00001)
		t.push_back(1.0);

	//for (int i = 0; i < t.size(); i++)
	//	std::cout << t[i] << std::endl;

	Bmat.resize(t.size(), n);
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < t.size(); i++) {
			double c = CGeoCalculator::nchoosek(n - 1, j);
			Bmat(i, j) = c*pow(t[i], j)*pow((1.0-t[i]), (n - 1 - j));
		}
	}

	return t;
}

void CGeoCalculator::getBernsteinBasis(int n, std::vector<double> t, Eigen::MatrixXd & Bmat)
{
	Bmat.resize(t.size(), n);
	Bmat.setZero();
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < t.size(); i++) {
			double c = CGeoCalculator::nchoosek(n - 1, j);
			Bmat(i, j) = c*pow(t[i], j)*pow((1.0 - t[i]), (n - 1 - j));
		}
	}
}

void CGeoCalculator::getSampleFromBezier(
	Eigen::MatrixXd & sample_pts, 
	Eigen::MatrixXd & ctrl_pts,
	double step)
{
	int numCPs = ctrl_pts.rows();
	Eigen::MatrixXd Bmat;
	getBernsteinBasis(numCPs, step, Bmat);
	sample_pts = Bmat*ctrl_pts;
}

void CGeoCalculator::getEqualArcLengthSampleFromBezier(
	std::vector<COpenMeshT::Point>& om_sample_pts, 
	std::vector<double>& tparas,
	const std::vector<COpenMeshT::Point>& om_ctrl_pts, int Nsamples)
{
	om_sample_pts.clear();
	tparas.clear();

	Eigen::MatrixXd ctrl_pts;
	convertOMptToEigenMat(om_ctrl_pts, ctrl_pts);

	int numCPs = ctrl_pts.rows();
	Eigen::MatrixXd Bmat;
	double dt = 0.0001;
	std::vector<double> ts;
	ts = getBernsteinBasis(numCPs, dt, Bmat);
	// get total arc length
	Eigen::MatrixXd dense_pts;
	dense_pts = Bmat*ctrl_pts;
	double total_length = 0.0;
	for (int i = 1; i < dense_pts.rows(); i++)
	{
		total_length += (dense_pts.row(i) - dense_pts.row(i - 1)).norm();
	}
	double dD = total_length / double(Nsamples);
	double accum_length = 0.0;
	std::vector<double> t_samples;
	t_samples.push_back(0);
	for (int i = 1; i < dense_pts.rows(); i++)
	{
		if (accum_length > dD) {
			//std::cout << accum_length << ", " << dD << std::endl;
			t_samples.push_back(ts[i]);
			accum_length = 0.0;
		}
		else {
			accum_length += (dense_pts.row(i) - dense_pts.row(i - 1)).norm();
		}
	}
	t_samples.push_back(1.0);
	getBernsteinBasis(numCPs, t_samples, Bmat);
	Eigen::MatrixXd sample_pts;
	sample_pts = Bmat*ctrl_pts;

	// convert sample_pts to om_pts;
	convertEigenMatToOMpt(sample_pts, om_sample_pts);
	tparas = t_samples;
}

double CGeoCalculator::nchoosek(int n, int k)
{
	if (k < 0 || k > n) {
		std::cerr << "Error: k should be less than n" << std::endl;
		return 0;
	}
	if (k == 0) {
		return 1;
	}
	return (n * nchoosek(n - 1, k - 1)) / k;
}

void CGeoCalculator::convertOMptToEigenMat(const std::vector<COpenMeshT::Point> om_pts, Eigen::MatrixXd & Mat)
{
	std::vector<double> om_pt_vec;
	for (int i = 0; i < om_pts.size(); i++)
	{
		om_pt_vec.push_back(om_pts[i][0]);
		om_pt_vec.push_back(om_pts[i][1]);
		om_pt_vec.push_back(om_pts[i][2]);
	}
	Mat = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>
		(om_pt_vec.data(), 3, om_pt_vec.size() / 3);
	Mat.transposeInPlace();
}

void CGeoCalculator::convertEigenMatToOMpt(const Eigen::MatrixXd & Mat, std::vector<COpenMeshT::Point>& om_pts)
{
	om_pts.clear();
	for (int i = 0; i < Mat.rows(); i++)
	{
		COpenMeshT::Point pt(Mat(i, 0), Mat(i, 1), Mat(i, 2));
		om_pts.push_back(pt);
	}
}

void CGeoCalculator::pts_sorting_alg(std::vector<COpenMeshT::Point>& pts)
{
	// use delaunay to compute the visibility graph
	double xCoord = pts[0][0];
	DelaunayT dt;
	for (int i = 0; i < pts.size(); i++)
	{
		KPoint2 kp(pts[i][1], pts[i][2]);
		dt.insert(kp);
	}

	COpenMeshT triMesh;
	std::map<DelaunayT::Vertex_handle, COpenMeshT::VertexHandle> cgal_om_vh_map;
	std::vector<std::vector<DelaunayT::Vertex_handle>> f_vhandles_cgal;
	for (DelaunayT::Finite_faces_iterator fiter = dt.finite_faces_begin();
		fiter != dt.finite_faces_end(); ++fiter)
	{
		std::vector<COpenMeshT::VertexHandle> vhandles;
		for (int i = 0; i < 3; i++) {
			DelaunayT::Vertex_handle vh_cgal = fiter->vertex(i);
			if (cgal_om_vh_map.find(vh_cgal) == cgal_om_vh_map.end()) {
				COpenMeshT::VertexHandle vh_om =
					triMesh.add_vertex(COpenMeshT::Point(0, vh_cgal->point().x(), vh_cgal->point().y()));
				cgal_om_vh_map[vh_cgal] = vh_om;
				vhandles.push_back(vh_om);
			}
			else
			{
				COpenMeshT::VertexHandle vh_om = cgal_om_vh_map[vh_cgal];
				vhandles.push_back(vh_om);
			}
		}
		triMesh.add_face(vhandles);
	}

	COpenMeshT::VertexHandle seed_vh;
	for (COpenMeshT::VertexIter viter = triMesh.vertices_begin(); viter != triMesh.vertices_end(); ++viter)
	{
		if (triMesh.is_boundary(*viter))
		{
			seed_vh = *viter;
			break;
		}
	}

	std::vector<COpenMeshT::VertexHandle> b_vertices;
	auto next_vh = seed_vh;
	do {
		b_vertices.push_back(next_vh);
		for (COpenMeshT::VertexOHalfedgeCCWIter voh_ccwiter = triMesh.voh_ccwbegin(next_vh);
			voh_ccwiter != triMesh.voh_ccwend(next_vh); ++voh_ccwiter)
		{
			if (triMesh.is_boundary(*voh_ccwiter))
			{
				next_vh = triMesh.to_vertex_handle(*voh_ccwiter);
				break;
			}
		}
	} while (next_vh != seed_vh);

	// output
	pts.clear();
	//std::cout << b_vertices.size() << std::endl;
	for (int i = 0; i < b_vertices.size(); i++)
	{
		COpenMeshT::Point p = triMesh.point(b_vertices[i]);
		p[0] = xCoord;
		pts.push_back(p);
		//std::cout << p[1] << " " << p[2] << std::endl;
	}

}

void CGeoCalculator::simplify_polygon(std::vector<COpenMeshT::Point>& pts, double dT)
{
	namespace PS = CGAL::Polyline_simplification_2;
	typedef PS::Stop_above_cost_threshold		Stop;
	typedef PS::Squared_distance_cost			Cost;

	double xCoord = pts[0][0];
	KPolygon2 polygon;
	for (int i = 0; i < pts.size(); i++)
		polygon.push_back(KPoint2(pts[i][1], pts[i][2]));

	Cost cost;
	polygon = PS::simplify(polygon, cost, Stop(dT));

	pts.clear();
	for (int i = 0; i < polygon.size(); i++) {
		pts.push_back(COpenMeshT::Point(xCoord, polygon[i].x(), polygon[i].y()));
	}
}

void CGeoCalculator::sample_polygon(std::vector<COpenMeshT::Point>& pts, double spacing, bool is_closed)
{
	if (is_closed)
		pts.push_back(pts.front());
	std::vector<COpenMeshT::Point> samples;

	// step size = ssz
	double alength = 0.0;
	for (int i = 0; i < pts.size() - 1; i++)
	{
		alength += (pts[i] - pts[i + 1]).norm();
	}
	double ssz = alength*spacing;

	// sample
	for (int i = 0; i < pts.size() - 1; i++)
	{
		double t_alength = (pts[i] - pts[i + 1]).norm();
		int cnt = 0;
		double t = 0.0;
		while (t < 1.0)
		{
			COpenMeshT::Point pt = pts[i] + t * (pts[i + 1] - pts[i]);
			if ((pt - pts[i + 1]).norm() > ssz) {
				samples.push_back(pt);
			}
			cnt++;
			t = double(cnt)*ssz / t_alength;
		}
	}

	//if (is_closed)
	//	samples.push_back(samples.front());
	pts = samples;
}

void CGeoCalculator::reconstruct_curve_from_pointset(std::vector<COpenMeshT::Point> &pts, float tolOMT)
{
	typedef CGAL::Optimal_transportation_reconstruction_2<K>    Otr_2;
	std::vector<KPoint2> points;
	for (int i = 0; i < pts.size(); i++)
		points.push_back(KPoint2(pts[i][1], pts[i][2]));

	std::cout << "init pts: " << points.size() << std::endl;
	if (points.size() == 0)
		return;

	std::vector<double> obb;
	computeOBB(pts, obb);
	double dl = (std::max)((obb[1] - obb[0]) / 2.,
		(obb[3] - obb[2]) / 2.);
	tolOMT = dl*tolOMT;
	std::cout << "tolOMT: " << tolOMT << std::endl;
	Otr_2 otr2(points);
	//otr2.set_random_sample_size(15);
	//otr2.set_relevance(0.3);
	std::vector<KPoint2> ptss;
	std::vector<std::size_t> isolated_vertices;
	std::vector<std::pair<std::size_t, std::size_t> > edges;
	int cntIter = 0;
	//float tolOMT = 0.6;
	otr2.run_under_wasserstein_tolerance(tolOMT);

	otr2.indexed_output(
		std::back_inserter(ptss),
		std::back_inserter(isolated_vertices),
		std::back_inserter(edges));
	std::cout << "(-------------List output---------- )" << tolOMT << " " << cntIter << std::endl;
	// points
	std::vector<KPoint2>::iterator pit;
	for (pit = ptss.begin(); pit != ptss.end(); pit++)
		std::cout << *pit << std::endl;
	// isolated vertices
	std::vector<std::size_t>::iterator vit;
	for (vit = isolated_vertices.begin(); vit != isolated_vertices.end(); vit++)
		std::cout << "1 " << *vit << std::endl;
	// edges
	std::vector<std::pair<std::size_t, std::size_t> >::iterator eit;
	for (eit = edges.begin(); eit != edges.end(); eit++)
		std::cout << "2 " << eit->first << " " << eit->second << std::endl;

	//KPoint2 orig(0.0, 0.0);
	//for (eit = edges.begin(); eit != edges.end(); eit++)
	//{
	//	KVector2 v12 = ptss[eit->second] - ptss[eit->first];
	//	KVector2 v01 = ptss[eit->first] - orig;
	//	double a = v01.x*v12*
	//}


	

}