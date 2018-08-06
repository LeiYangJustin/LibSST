#include "euclidean_MST.h"
#include "cgal_predef.h"

void CfindChainLoopUsingEMST::delaunay(std::vector<COpenMeshT::Point> pts, std::vector<Edge> & edges)
{
	std::map<DelaunayT::Vertex_handle, int> vh_id_map;
	edges.clear();
	//
	double xCoord = pts[0][0];
	DelaunayT dt;
	for (int i = 0; i < pts.size(); i++)
	{
		KPoint2 kp(pts[i][1], pts[i][2]);
		DelaunayT::Vertex_handle vh = dt.insert(kp);
		vh_id_map.insert(std::make_pair(vh, i));
	}
	//
	for (DelaunayT::Finite_edges_iterator eiter = dt.finite_edges_begin();
		eiter != dt.finite_edges_end(); ++eiter)
	{
		// 
		DelaunayT::Face_handle fh = eiter->first;
		int vid = eiter->second;
		DelaunayT::Vertex_handle vsrc = fh->vertex(fh->cw(vid));
		DelaunayT::Vertex_handle vtgt = fh->vertex(fh->ccw(vid));
		//
		Node v1(vh_id_map[vsrc], COpenMeshT::Point(0, vsrc->point()[0], vsrc->point()[1]));
		Node v2(vh_id_map[vtgt], COpenMeshT::Point(0, vtgt->point()[0], vtgt->point()[1]));
		Edge e(v1, v2);
		edges.push_back(e);
	}
}

void CfindChainLoopUsingEMST::Compute()
{
	// delaunay point set
	std::vector<Edge> edges;
	delaunay(pts_, edges);
	// sort edges in the increasing order
	std::priority_queue<std::pair<double, int>, 
		std::vector<std::pair<double, int>>, 
		std::greater<std::pair<double, int>>> p_queue;
	for (int i = 0; i < edges.size(); i++)
		p_queue.push(std::make_pair(edges[i].GetWeightedLength(), i));
	// incrementally add the shortest edges to s_
	while (!p_queue.empty())
	{
		//std::cout << p_queue.top().first << std::endl;
		Edge e = edges[p_queue.top().second];
		p_queue.pop();
		s_.AddEdge(e);
	}
	//std::cout << "# edges in triangulation: " << edges.size() << std::endl;
	//std::cout << "# edges in subgraph: " << s_.edges_.size() << std::endl;

	s_.Pruning();
	//std::cout << "# edges after pruning: " << s_.edges_.size() << std::endl;

	s_.SortNodes(sorted_pid_);

	//std::ofstream pts_file("pts.txt");
	//if (pts_file.is_open())
	//{
	//	for (int i = 0; i < pts_.size(); i++) {
	//		pts_file << pts_[i] << std::endl;
	//	}
	//}
	//pts_file.close();

	//std::ofstream edges_file("edges_file.txt");
	//if (edges_file.is_open())
	//{
	//	for (int i = 0; i < edges.size(); i++) {
	//		edges_file << edges[i].GetSrcPt() << ", " << edges[i].GetTgtPt() << std::endl;
	//	}
	//}
	//edges_file.close();

	//std::ofstream mst_file("mst_file.txt");
	//if (mst_file.is_open())
	//{
	//	for (int i = 0; i < s_.edges_.size(); i++) {
	//		mst_file << s_.edges_[i].GetSrcPt() << ", " << s_.edges_[i].GetTgtPt() << std::endl;
	//	}
	//}
	//mst_file.close();




	
}
