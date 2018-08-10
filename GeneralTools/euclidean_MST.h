#ifndef C_EUCLIDEAN_MINIMAL_SPANNING_TREE_H
#define C_EUCLIDEAN_MINIMAL_SPANNING_TREE_H

#include "general_tool_prereq.h"
#include "../DataColle/mesh_object.h"
#include "geo_calculator.h"
#include <queue>

class Node
{
public:
	Node(int id, COpenMeshT::Point pt, 
		COpenMeshT::Point vec = COpenMeshT::Point(0.0, 0.0, 0.0)) 
	{
		point_ = pt;
		normal_ = vec;
		id_ = id;
	}
	~Node() {};

	COpenMeshT::Point	point_;
	COpenMeshT::Point	normal_;
	int					id_;
};

class Edge
{
public:
	Edge(Node a, Node b)
	{
		length_ = (a.point_ - b.point_).norm();
		weight_ = 1.0; // normal weight?
		src_id_ = a.id_;
		p_src_ = a.point_;
		tgt_id_ = b.id_;
		p_tgt_ = b.point_;
	};
	~Edge() {} ;
	double GetEdgeLength() { return length_; };
	double GetWeightedLength() { return weight_*length_; };
	int GetSrcId() { return src_id_; };
	int GetTgtId() { return tgt_id_; };
	COpenMeshT::Point GetSrcPt() { return p_src_; };
	COpenMeshT::Point GetTgtPt() { return p_tgt_; };
private:
	double length_;
	double weight_;
	COpenMeshT::Point p_src_;
	COpenMeshT::Point p_tgt_;
	int src_id_;
	int tgt_id_;
};

class CTree
{
public:
	CTree() {};
	~CTree() {};

	bool AddEdge(Edge e)
	{
		int src_id = e.GetSrcId();
		int tgt_id = e.GetTgtId();
		int s_set_id = find_set(src_id);
		int t_set_id = find_set(tgt_id);
	
		if (s_set_id == -1 && t_set_id == -1)
		{
			// add as a new set
			edges_.push_back(e);
			std::set<int> tmp_set;
			tmp_set.insert(src_id);
			tmp_set.insert(tgt_id);
			set_list_.push_back(tmp_set);
			return true;
		}
		else if (s_set_id == -1) {
			edges_.push_back(e);
			set_list_[t_set_id].insert(src_id);
			return true;
		}
		else if (t_set_id == -1) {
			edges_.push_back(e);
			set_list_[s_set_id].insert(tgt_id);
			return true;
		}
		else if (t_set_id != s_set_id)
		{
			// merge set
			edges_.push_back(e);
			set_list_[s_set_id].insert(set_list_[t_set_id].begin(), set_list_[t_set_id].end());
			set_list_[t_set_id].clear();
			return true;
		}
		else {
			return false;
		}
	};
	void Pruning()
	{
		int cntIter = 0;
		int cntRemoval;
		do {
			// init
			std::map<int, int> map_id_ndeg;
			compute_node_degs(map_id_ndeg);
			cntRemoval = 0;
			// processing
			std::vector<Edge> new_edges;
			for (int i = 0; i < edges_.size(); i++)
			{
				Edge e = edges_[i];
				if (map_id_ndeg[e.GetSrcId()] == 1 && map_id_ndeg[e.GetTgtId()] > 2) {
					cntRemoval++;
				}
				else if (map_id_ndeg[e.GetSrcId()] > 2 && map_id_ndeg[e.GetTgtId()] == 1) {
					cntRemoval++;
				}
				else {
					new_edges.push_back(e);
				}
			}
			edges_.clear();
			edges_ = new_edges;
			//std::cout << "# removal: " << cntRemoval << "; # edges: " << edges_.size() << std::endl;
		} while (cntRemoval > 0 && cntIter++ < 10);
	}
	void SortNodes(std::vector<int> & sorted_nodes)
	{
		// find free end
		std::vector<int> endPids;
		std::map<int, int> map_id_ndeg;
		compute_node_degs(map_id_ndeg);
		int numJoints = 0;
		for (auto mapIter = map_id_ndeg.begin();
			mapIter != map_id_ndeg.end(); ++mapIter)
		{
			if (mapIter->second == 1)
				endPids.push_back(mapIter->first);
			else if (mapIter->second > 2)
				numJoints++;
		}

		// should be endPids.size() - 2 = numJoints
		//std::cerr << "#End points: " << endPids.size() << ", #Joints: " << numJoints << std::endl;
		
		// traverse to the other end
		int start_id = endPids.front();
		//std::vector<int> sorted_nodes;
		std::set<int> visitor_set;
		int cntIter = 0;
		while (visitor_set.size() < edges_.size() && cntIter++ < 5000)
		{
			for (int ei = 0; ei < edges_.size(); ei++)
			{
				if (visitor_set.find(ei) == visitor_set.end()) {
					// 
					int sid = edges_[ei].GetSrcId();
					int tid = edges_[ei].GetTgtId();
					if (start_id == sid) {
						visitor_set.insert(ei);
						sorted_nodes.push_back(sid);
						start_id = tid;
						break;
					}
					else if (start_id == tid) {
						visitor_set.insert(ei);
						sorted_nodes.push_back(tid);
						start_id = sid;
						break;
					}
				}
			}
			//std::cout << "# visited edges: " << visitor_set.size() << std::endl;
		}
	};

public:
	std::vector<Edge> edges_;
	std::vector<std::set<int>> set_list_;


private:
	void compute_node_degs(std::map<int, int> & map_id_ndeg)
	{
		map_id_ndeg.clear();
		for (int i = 0; i < edges_.size(); i++)
		{
			if (map_id_ndeg.find(edges_[i].GetSrcId()) == map_id_ndeg.end()) {
				map_id_ndeg[edges_[i].GetSrcId()] = 1;
			}
			else {
				map_id_ndeg[edges_[i].GetSrcId()]++;
			}
			if (map_id_ndeg.find(edges_[i].GetTgtId()) == map_id_ndeg.end()) {
				map_id_ndeg[edges_[i].GetTgtId()] = 1;
			}
			else {
				map_id_ndeg[edges_[i].GetTgtId()]++;
			}
		}
	};
	int find_set(int s)
	{
		int s_set = -1, t_set = -1;
		for (int i = 0; i < set_list_.size(); i++)
		{
			if (set_list_[i].find(s) != set_list_[i].end()) {
				s_set = i;
				break;
			}
		}
		return s_set;
	}
};

class GENERAL_TOOLS_CLASS CfindChainLoopUsingEMST
{
public:
	CfindChainLoopUsingEMST() {};
	~CfindChainLoopUsingEMST() {};
	void SetPoints(std::vector<COpenMeshT::Point> pts)
	{
		pts_ = pts;
	};
	void GetSortedPoints(std::vector<COpenMeshT::Point> & spts)
	{
		spts.clear();
		for (int i = 0; i < sorted_pid_.size(); i++)
		{
			assert(sorted_pid_[i] < pts_.size());
			spts.push_back(pts_[sorted_pid_[i]]);
		}
	};
	void GetSortedPids(std::vector<int> & sorted_pid)
	{
		sorted_pid = sorted_pid_;
	};
	void Compute();
private:
	CTree s_;
	std::vector<COpenMeshT::Point> pts_;
	std::vector<int> sorted_pid_;

	void delaunay(std::vector<COpenMeshT::Point> pts, std::vector<Edge> & edges);
};
#endif //C_EUCLIDEAN_MINIMAL_SPANNING_TREE_H