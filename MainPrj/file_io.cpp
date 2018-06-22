#include "file_io.h"
#include "../DataColle/data_io.h"

SwitchType CFileIO::read_config_file(std::string fname)
{
	//Switch = 1 -> global
	//Switch = 2 -> local
	//Switch = 0 -> pending
	//Switch = -1 -> exit

	std::ifstream in_file(fname);
	if (!in_file.is_open())
	{
		std::cerr << "Failed to open the config_switch file" << std::endl;
		return SwitchType::sstExit;
	}
	std::string line;
	std::getline(in_file, line);
	//while (std::getline(in_file, line))
	{
		std::cout << line << std::endl;
		if (strcmp(line.c_str(), "GLOBAL") == 0) {
			return SwitchType::sstGlobal; // global
		}
		else if (strcmp(line.c_str(), "LOCAL") == 0) {
			return SwitchType::sstLocal; // local 
		}
		else if (strcmp(line.c_str(), "PENDING") == 0) {
			return  SwitchType::sstPending; // pending
		}
		else {
			return SwitchType::sstExit; // exit
		}
	}
	in_file.close();
	return true;
}

bool CFileIO::read_local_def_cs_id_from_config_file(std::string fname, std::vector<int> &cs_ids)
{
	std::ifstream in_file(fname);
	if (!in_file.is_open())
	{
		std::cerr << "Failed to open the config_switch file" << std::endl;
		return false;
	}
	std::string line;
	std::getline(in_file, line);
	std::getline(in_file, line);
	std::istringstream ss(line);
	// coord
	while(ss.good())
	{
		int a;
		ss >> a;
		cs_ids.push_back(a);
	}
	in_file.close();
	return true;
}

bool CFileIO::read_cross_sections_file(std::string fname, std::vector<CCrossSection>& cs_list)
{
	cs_list.clear();

	std::ifstream in_file(fname);
	if (!in_file.is_open())
	{
		std::cerr << "Failed to open the CS file" << std::endl;
		return false;
	}
	std::string line;
	std::getline(in_file, line);
	if (strcmp(line.c_str(), "##CrossSection") != 0)
	{
		std::cerr << "Failed to open the CS file" << std::endl;
		return false;
	}

	int flag = 0; // 1- points; 2-rmf;
	//std::vector<double> tmpCoord, tmpArcLength, tmpRMF;
	std::map<int, COpenMeshT::Point> id_pt_map;
	
	std::vector<COpenMeshT::Point> pt_list;
	while (std::getline(in_file, line))
	{
		if (strcmp(line.c_str(), "#Points") == 0) {
			flag = 1;
			continue;
		}
		else if (strcmp(line.c_str(), "#CS") == 0) {
			flag = 2;
			CCrossSection cs;
			if (pt_list.size() > 0)
			{
				cs.SetEmbProfPts(pt_list); // set cs
				cs_list.push_back(cs); // add cs
				pt_list.clear(); // clear pt for another run
			}
			continue;
		}
		else if (strcmp(line.c_str(), "#END") == 0)
		{
			std::cout << "skeleton reading finished" << std::endl;
			break;
		}
		else if (flag == 1) {
			std::istringstream ss(line);
			// id
			int id;
			ss >> id;
			// coord
			COpenMeshT::Point pt(0.0, 0.0, 0.0);
			for (int i = 0; i < 3; i++)
			{
				double a;
				ss >> a;
				pt[i] = a;
			}
			id_pt_map[id] = pt;
		}
		else if (flag == 2) {
			std::istringstream ss(line);
			int a;
			do {
				ss >> a;
				pt_list.push_back(id_pt_map[a]);
			} while (!ss.eof());
		}
	}
	in_file.close();
	return true;
}

bool CFileIO::read_cross_section_pts(std::string fname, std::vector<COpenMeshT::Point>& cs_pts)
{
	std::ifstream in_file(fname);
	if (!in_file.is_open())
	{
		std::cerr << "Failed to open the CS file" << std::endl;
		return false;
	}
	std::string line;
	while (std::getline(in_file, line)) {
		std::istringstream ss(line);
		// coord
		COpenMeshT::Point pt(0.0, 0.0, 0.0);
		for (int i = 0; i < 3; i++)
		{
			double a;
			ss >> a;
			pt[i] = a;
		}
		cs_pts.push_back(pt);
	};
	return true;
}

bool CFileIO::read_skeleton(std::string fname, CSkeleton &s)
{
	CCurveObject* p_def_skeleton = new CCurveObject;
	CDataIO::ReadCurve(fname, *p_def_skeleton);
	CSkeleton s_tmp(p_def_skeleton->GetCurve());
	s.CopyFrom(s_tmp);

	return true;
}
