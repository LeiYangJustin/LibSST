#ifndef C_SORTING_2D_PTS
#define C_SORTING_2D_PTS

#include "general_tool_prereq.h"
#include "../DataColle/mesh_object.h"
//#include <Eigen/Core>
#include "cgal_predef.h"

class GENERAL_TOOLS_CLASS CSorting2DPts
{
public:
	CSorting2DPts(std::vector<COpenMeshT::Point> &pt_list, int verbose_Lvl, double dT = 0.5);
	~CSorting2DPts();
	
	//void construct_all(std::vector<COpenMeshT::Point> &pt_list, int verbose_Lvl, double dT = 0.5);

private:
	int verbose_Lvl_;
	double dThreshold_;

private:
	void construct_all(std::vector<COpenMeshT::Point> &pt_list);
	double compute_point_to_segment(const KPoint2 p, const KSegment2 s, KPoint2 &po);
	void simplify_polygon(std::vector<KPoint2> &pts, double dT = 0.5);
};

#endif //C_SORTING_2D_PTS