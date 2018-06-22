#ifndef C_VISIBILITY_GRAPH
#define C_VISIBILITY_GRAPH

#include "general_tool_prereq.h"
#include "../DataColle/mesh_object.h"
//#include <Eigen/Core>
#include "cgal_predef.h"

// the class is developed mainly based on the paper http://webee.technion.ac.il/~ayellet/Ps/17-KT.pdf

class GENERAL_TOOLS_CLASS CVisibilityGraph
{
public:
	CVisibilityGraph();
	~CVisibilityGraph();
	void construct_vgraph(std::vector<COpenMeshT::Point> pt_list, Eigen::MatrixXi &A, double gamma = -1.0);
	void construct_all(std::vector<COpenMeshT::Point> &pt_list, double gamma = -1.0);


private:
	void construct_vgraph(std::vector<KPoint2> pt_list, Eigen::MatrixXi &A, double gamma);
	void ghpr_operater(KPoint2 c, std::vector<KPoint2> in_pts, std::vector<int> &out_tags, double gamma = -1.0);	

	struct compare_id_dist_pair {
		bool operator()(const std::pair<int, double> &lhs,
			const std::pair<int, double> &rhs) const {
			return lhs.second > rhs.second;
		}
	};

};


#endif // ! C_VISIBILITY_GRAPH


