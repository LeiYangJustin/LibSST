#include "sorting_2d_pts.h"
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/Polyline_simplification_2/Squared_distance_cost.h>

CSorting2DPts::CSorting2DPts(std::vector<COpenMeshT::Point>& pt_list, int verbose_Lvl, double dT)
	: verbose_Lvl_(verbose_Lvl), dThreshold_(dT)
{
	construct_all(pt_list);
}

CSorting2DPts::~CSorting2DPts()
{
}

void CSorting2DPts::construct_all(std::vector<COpenMeshT::Point>& pt_list)
{
	// CONVERT
	double xCoord = pt_list[0][0];
	std::vector<KPoint2> in_pts, cv_pts;
	for (int i = 0; i < pt_list.size(); i++) {
		in_pts.push_back(KPoint2(pt_list[i][1], pt_list[i][2]));
	}

	// COMPUTE CONVEX HULL
	CGAL::convex_hull_2(in_pts.begin(), in_pts.end(), std::back_inserter(cv_pts));
	std::vector<KSegment2> segment_list;
	for (int j = 0; j < cv_pts.size(); j++) {
		if (j == cv_pts.size() - 1)
			segment_list.push_back(KSegment2(cv_pts[j], cv_pts[0]));
		else
			segment_list.push_back(KSegment2(cv_pts[j], cv_pts[j+1]));
	}

	// ITERATIVELY REFINE THE SHAPE
	int cntIter = 0, maxIter = 5;
	while (cntIter++ < maxIter) {
		std::cout << "Iter: " << cntIter << std::endl;
		double maxDist = 0;
		// find inserted position
		std::map<int, int, std::greater<int>> segment_inpid_map;
		std::map<int, double, std::greater<int>> segment_dist_map;
		for (int i = 0; i < in_pts.size(); i++) 
		{
			double dmin = std::numeric_limits<double>::max();
			int idmin = -1;
			for (int j = 0; j < segment_list.size(); j++)
			{
				KPoint2 po(0.0, 0.0);
				double d = compute_point_to_segment(in_pts[i], segment_list[j], po);
				if (d < dmin) {
					dmin = d;
					idmin = j;
				}
			}
			//std::cout << i << "-" << idmin << ": " << dmin << std::endl;
			if (maxDist < dmin)
				maxDist = dmin;

			if (dmin > dThreshold_)
			{
				if (segment_inpid_map.find(idmin) == segment_inpid_map.end()) {
					segment_inpid_map[idmin] = i;
					segment_dist_map[idmin] = dmin;
				}
				else if (segment_dist_map[idmin] < dmin) {
					segment_inpid_map[idmin] = i;
					segment_dist_map[idmin] = dmin;
				}
			}
		}
		if (segment_inpid_map.size() == 0) {
			std::cout << "maxDist = " << maxDist << " is less than dThreshold" << std::endl;
			break;
		}
		// insert segments
		if (verbose_Lvl_ == 1)
			std::cout << "#segment = " << segment_list.size() << std::endl;

		for (auto itmap = segment_inpid_map.begin();
			itmap != segment_inpid_map.end(); ++itmap)
		{
			auto liter = segment_list.begin();
			for (int j = 0; j < segment_list.size(); j++, ++liter)
			{
				if (j == itmap->first)
				{
					KPoint2 src = segment_list[itmap->first].source();
					KPoint2 tgt = segment_list[itmap->first].target();
					KPoint2 inserted = in_pts[itmap->second];
					KSegment2 seg_prev(src, inserted);
					// replace
					segment_list[itmap->first] = seg_prev;
					// then insert
					KSegment2 seg_next(inserted, tgt);
					segment_list.insert(++liter, seg_next);
					break;
				}
			}
			if (verbose_Lvl_==1) {
				std::cout << "seg_id = " << itmap->first
					<< "(" << segment_list[itmap->first].source()
					<< ", " << segment_list[itmap->first].target() << ")"
					<< ", pid = " << itmap->second << std::endl;
			}
		}
	}
	std::vector<KPoint2> out_pts;
	for (int i = 0; i < segment_list.size(); i++) {
		out_pts.push_back(segment_list[i].source());
	}


	// correct the added points
	simplify_polygon(out_pts, dThreshold_);
	//pt_list.clear();
	for (int i = 0; i < out_pts.size(); i++)
		pt_list.push_back(COpenMeshT::Point(xCoord, out_pts[i].x(), out_pts[i].y()));

	//// OUTPUT INTERMEDIATE INFO
	//if (verbose_Lvl_ == 1) 
	//{
	//	std::ofstream fout_secpts;
	//	fout_secpts.open("fDataSkeleton/section_pts.txt");
	//	for (int i = 0; i < in_pts.size(); i++)
	//		fout_secpts << in_pts[i] << std::endl;
	//	fout_secpts.close();

	//	//
	//	std::ofstream fout_ch;
	//	fout_ch.open("fDataSkeleton/convex_hull.txt");
	//	for (int j = 0; j < cv_pts.size(); j++)
	//		fout_ch << cv_pts[j] << std::endl;
	//	fout_ch.close();

	//	//
	//	std::ofstream fout_ch_after;
	//	fout_ch_after.open("fDataSkeleton/convex_hull_after.txt");
	//	for (int i = 0; i < segment_list.size(); i++) {
	//		fout_ch_after << segment_list[i].source() << std::endl;
	//	}
	//	fout_ch_after.close();
	//}

}

double CSorting2DPts::compute_point_to_segment(const KPoint2 p, const KSegment2 s, KPoint2 & po)
{
	double dist = 0.0;

	// projection
	KPoint2 p0 = s.source();
	KPoint2 p1 = s.target();

	KVector2 v = p - p0;
	KVector2 sv = s.to_vector();
	
	po = p0 + sv*(v*sv) / sv.squared_length();
	double t = 0.0;
	if (s.to_vector().x() != 0)
		t = (po - p0).x() / s.to_vector().x();
	else
		t = (po - p0).y() / s.to_vector().y();

	if (t >= 0 && t <= 1)
		dist = sqrt((p - po).squared_length());
	else if (t < 0)
		dist = sqrt((p - p0).squared_length());
	else if (t > 1)
		dist = sqrt((p - p1).squared_length());

	return dist;
}

void CSorting2DPts::simplify_polygon(std::vector<KPoint2> &pts, double dT)
{
	namespace PS = CGAL::Polyline_simplification_2;
	typedef PS::Stop_above_cost_threshold		Stop;
	typedef PS::Squared_distance_cost			Cost;

	KPolygon2 polygon;
	polygon.insert(polygon.vertices_end(), pts.begin(), pts.end());
	Cost cost;
	polygon = PS::simplify(polygon, cost, Stop(dT));

	pts.clear();
	for (int i = 0; i < polygon.size(); i++) {
		pts.push_back(KPoint2(0, polygon[i].x(), polygon[i].y()));
	}

	//
	if (verbose_Lvl_ == 1) {
		std::ofstream fout_ch_after;
		fout_ch_after.open("fDataSkeleton/polygon.txt");
		for (int i = 0; i < polygon.size(); i++) {
			fout_ch_after << polygon[i] << std::endl;
		}
		fout_ch_after.close();
	}

}
