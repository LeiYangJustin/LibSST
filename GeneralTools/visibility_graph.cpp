#include "visibility_graph.h"

#include <queue>

#include <CGAL/convex_hull_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/Polyline_simplification_2/Squared_distance_cost.h>


CVisibilityGraph::CVisibilityGraph()
{
}


CVisibilityGraph::~CVisibilityGraph()
{
}

void CVisibilityGraph::ghpr_operater(KPoint2 c, std::vector<KPoint2> in_pts, std::vector<int> &out_tags, double gamma)
{
	// step 1: point transformation
	std::vector<KPoint2> q_list, out_pts;
	for (int i = 0; i < in_pts.size(); i++)
	{
		KVector2 t = in_pts[i] - c;
		double d = sqrt(t.squared_length());
		// inversion kernel
		q_list.push_back(KPoint2(0, 0) + t/d* pow(d, gamma));
	}
	// step 2: compute convex hull
	CGAL::convex_hull_2(q_list.begin(), q_list.end(), std::back_inserter(out_pts));

	//pull back
	for (int i = 0; i < out_pts.size(); i++)
	{
		double dmin = 10e-5;
		int idmin = -1;
		for (int j = 0; j < q_list.size(); j++)
		{
			double dtmp = (q_list[j] - out_pts[i]).squared_length();
			if (dtmp < dmin)
			{
				idmin = j;
				dmin = dtmp;
			}
		}
		if (idmin != -1)
			out_tags.push_back(idmin);
	}
}

void CVisibilityGraph::construct_vgraph(std::vector<KPoint2> pt_list, Eigen::MatrixXi &A, double gamma)
{
	typedef std::pair<int, double> Id_Dist_pair;
	typedef DelaunayT::Vertex_handle dtvh;
	// compute dt of pt_list
	// use delaunay to compute the visibility graph
	DelaunayT dt;
	std::map<DelaunayT::Vertex_handle, int> vh_vid_map;
	std::map<int, DelaunayT::Vertex_handle> vid_vh_map;
	for (int i = 0; i < pt_list.size(); i++) {
		DelaunayT::Vertex_handle vh = dt.push_back(pt_list[i]);
		vh_vid_map[vh] = i;
		vid_vh_map[i] = vh;
	}
	
	// construct bi-directional map between vh and vid
	A.resize(pt_list.size(), pt_list.size());
	A.setZero();
	std::map<dtvh, std::vector<dtvh>> vh_vcir_map;
	for (DelaunayT::Finite_vertices_iterator vit = dt.finite_vertices_begin();
		vit != dt.finite_vertices_end(); ++vit)
	{
		//int ir = vh_vid_map[vit];
		DelaunayT::Vertex_circulator vh_cir = dt.incident_vertices(vit), vcend = vh_cir;
		std::vector<dtvh> vcir_list;
		do {
			if (!dt.is_infinite(vh_cir)) {
				//int ic = vh_vid_map[vh_cir];
				//A(ir, ic) = 1;
				vcir_list.push_back(vh_cir);
			}
		} while (++vh_cir != vcend);
		//
		vh_vcir_map[vit] = vcir_list;
	}

	// removing edges from A if they are not visibile via GHPR
	for (auto itmap = vh_vcir_map.begin();
		itmap!= vh_vcir_map.end(); ++itmap) 
	{
		std::vector<int> in_vids;
		std::vector<KPoint2> in_pts;
		
		KPoint2 cp = itmap->first->point();
		int ir = vh_vid_map[itmap->first];

		std::vector<dtvh> vcir_list = itmap->second;
		for (int i = 0; i < vcir_list.size(); ++i) 
		{
			int ic = vh_vid_map[vcir_list[i]];
			in_vids.push_back(ic);
			in_pts.push_back(vcir_list[i]->point());
		}


		std::vector<int> out_tags;
		ghpr_operater(cp, in_pts, out_tags, gamma);

		// recording the edges
		std::priority_queue<Id_Dist_pair, std::vector<Id_Dist_pair>, compare_id_dist_pair> pq_edge;
		for (int i = 0; i < out_tags.size(); i++) {
			int ic = in_vids[out_tags[i]];
			double dpr = (cp - in_pts[out_tags[i]]).squared_length();
			pq_edge.push(std::make_pair(ic, dpr));
			//std::cout << ic << " - " << dpr << std::endl;
		}
		// keep only the shortest two edges and remove the others
		int cnt = 0;
		do {
			if (cnt++ < 2) {
				std::cout << ir << ", " << pq_edge.top().first << ": " << pq_edge.top().second << std::endl;
				A(ir, pq_edge.top().first) = 1;
				//A(pq_edge.top().first, ir) = 1;
			}
			pq_edge.pop();
		} while (pq_edge.size() > 0);
	}

	// output A
	for (int ir = 0; ir < A.rows()-1; ir++)
	{
		for (int ic = ir+1; ic < A.cols(); ic++)
		{
			if (A(ir, ic) != A(ic, ir))
			{
				A(ir, ic) = 0;
				A(ic, ir) = 0;
			}
		}
	}
}

void CVisibilityGraph::construct_vgraph(std::vector<COpenMeshT::Point> omesh_pts, Eigen::MatrixXi &A, double gamma)
{
	std::vector<KPoint2> pt_list;
	for (int i = 0; i < omesh_pts.size(); i++)
		pt_list.push_back(KPoint2(omesh_pts[i][1], omesh_pts[i][2]));

	construct_vgraph(pt_list, A, gamma);
}

void CVisibilityGraph::construct_all(std::vector<COpenMeshT::Point>& omesh_pts, double gamma)
{
	std::vector<KPoint2> pt_list;
	for (int i = 0; i < omesh_pts.size(); i++)
		pt_list.push_back(KPoint2(omesh_pts[i][1], omesh_pts[i][2]));

	KVector2 cpvec(0, 0);
	for (int i = 0; i < pt_list.size(); i++) {
		cpvec = cpvec + (pt_list[i] - KPoint2(0, 0));
	}
	cpvec = cpvec / double(pt_list.size());

	KPoint2 cp = KPoint2(0, 0) + cpvec;
	std::vector<int> out_tags;
	ghpr_operater(cp, pt_list, out_tags, gamma);

	
	std::vector<KPoint2> qt_list;
	for (int i = 0; i < out_tags.size(); i++) {
		qt_list.push_back(pt_list[out_tags[i]]);
		//omesh_pts.push_back(COpenMeshT::Point(0, pt_list[out_tags[i]].x(), pt_list[out_tags[i]].y()));
		//std::cout << pt_list[out_tags[i]] << std::endl;
	}

	namespace PS = CGAL::Polyline_simplification_2;
	typedef PS::Stop_above_cost_threshold		Stop;
	typedef PS::Squared_distance_cost			Cost;

	KPolygon2 polygon;
	polygon.insert(polygon.vertices_end(), qt_list.begin(), qt_list.end());
	Cost cost;
	polygon = PS::simplify(polygon, cost, Stop(0.5));

	omesh_pts.clear();
	for (int i = 0; i < polygon.size(); i++) {
		omesh_pts.push_back(COpenMeshT::Point(0, polygon[i].x(), polygon[i].y()));
		//std::cout << pt_list[out_tags[i]] << std::endl;
	}
}
