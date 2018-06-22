#include "fe_meshIO.h"

bool CFEMeshIO::readKeyFile(std::string filename, 
	std::map<int, std::vector<double>> &Vmap, 
	std::vector<std::vector<int>>& F)
{
	// 
	std::ifstream inFile;
	inFile.open(filename);
	if (!inFile) {
		std::cerr << "Unable to open file " << filename << std::endl;
		return false;
	}
	else {
		std::cerr << "Opened file " << filename << std::endl;
	}

	std::vector<std::string> keyword_list;
	keyword_list.push_back("*END");
	keyword_list.push_back("*NODE");
	keyword_list.push_back("*ELEMENT_SHELL");
	
	std::vector<int> vid_list;
	std::vector<int> fid_list;

	int flag = -1;
	while (!inFile.eof()) {
		std::string line;
		getline(inFile, line);
		if (inFile.fail()) break;

		// $$ comment lines
		if (line[0] == '$')
			continue;
		// * keyword lines
		else if (line[0] == '*') {
			flag = -1; //reset
			for (int i = 0; i < keyword_list.size(); i++) {
				if (line == keyword_list[i]) {
					flag = i; break;
				}
			}
		}
		// data line
		else {
			// exit loop
			if (flag == 0) {
				break;
			}
			// read vertices
			else if (flag == 1) {
				std::vector<double> data;
				std::istringstream iss(line);
				int id;
				iss >> id;
				double d;
				for(int i = 0; i < 3; i++)
				{
					iss >> d;
					data.push_back(d);
				}
				Vmap.insert(std::make_pair(id, data));
			}
			// read shell elements
			else if (flag == 2) {
				std::vector<int> data;
				std::istringstream iss(line);
				int d;
				iss >> d; // element_id
				iss >> d; // part_id
				for (int i = 0; i < 4; i++)
				{
					iss >> d;
					data.push_back(d);
				}
				F.push_back(data);
			}
		}
	}
	inFile.close();
	return true;
}

bool CFEMeshIO::readKeyFile(std::string filename, COpenMeshT & mesh)
{
	std::map<int, std::vector<double>> nodes;
	std::vector<std::vector<int>> elements;
	readKeyFile(filename, nodes, elements);
	
	std::cout << "#nodes = " << nodes.size() << "; #elements = " << elements.size() << std::endl;

	// make mesh
	mesh.clear();
	std::map<int, COpenMeshT::VertexHandle> pid_vh_map;
	for (int i = 0; i < elements.size(); i++)
	{
		std::vector<COpenMeshT::VertexHandle> handles;
		for (int j = 0; j < 4; j++) {
			int pid = elements[i][j];
			if (pid_vh_map.find(pid) == pid_vh_map.end()) {
				auto p = COpenMeshT::Point(nodes[pid][0], nodes[pid][1], nodes[pid][2]);
				COpenMeshT::VertexHandle vh = mesh.add_vertex(p);
				mesh.data(vh).set_vlabel(pid);
				pid_vh_map.insert(std::make_pair(pid, vh));
			}
			handles.push_back(pid_vh_map[pid]);
		}
		std::vector<COpenMeshT::VertexHandle> fvhandles;
		fvhandles.push_back(handles[0]);
		fvhandles.push_back(handles[1]);
		fvhandles.push_back(handles[2]);
		mesh.add_face(fvhandles);
		if (handles[2] != handles[3]) {
			fvhandles.clear();
			fvhandles.push_back(handles[0]);
			fvhandles.push_back(handles[2]);
			fvhandles.push_back(handles[3]);
			mesh.add_face(fvhandles);
		}
	}

	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter)
	{
		int cnt = 0;
		for (auto vviter = mesh.vv_ccwiter(*viter); vviter.is_valid(); ++vviter)
			cnt++;
		if (cnt == 0)
			mesh.delete_vertex(*viter);
	}
	mesh.garbage_collection();

	return true;
}

bool CFEMeshIO::updateNodesInKeyFile(
	std::string filename_src, std::string filename_tgt, 
	std::map<int, std::vector<double>> &Vmap)
{
	// 
	std::ifstream srcFile;
	srcFile.open(filename_src);
	if (!srcFile) {
		std::cerr << "Unable to open file " << filename_src << std::endl;
		return false;
	}

	std::ofstream outFile(filename_tgt, std::ios::out);
	if (!outFile.is_open()) {
		std::cerr << "Unable to open file " << filename_tgt << std::endl;
		return false;
	}

	int cntNodes = 0;
	int flag = -1;
	while (!srcFile.eof()) {
		std::string line;
		getline(srcFile, line);
		
		// check
		if (srcFile.fail()) break;

		// $$ comment lines
		if (line[0] == '$') {
			outFile << line << std::endl;
		}
		// * keyword lines
		else if (line[0] == '*') {
			outFile << line << std::endl;
			flag = -1; //reset
			if ("*NODE" == line)
				flag = 1;
		}
		// data line
		else {
			// copy and update vertices
			if (flag == 1) {
				std::istringstream iss(line);
				int id;
				iss >> id;
				outFile << std::setw(8) << id
					<< std::setw(16) << Vmap[id][0]
					<< std::setw(16) << Vmap[id][1]
					<< std::setw(16) << Vmap[id][2]
					<< std::endl;
				cntNodes++;
			}
			// COPY other lines
			else {
				outFile << line << std::endl;
			}
		}
	}
	srcFile.close();
	outFile.close();
	return true;
}

bool CFEMeshIO::readFEResultSECFORC(std::string fpath, double & data)
{
	std::string filename = fpath;
	filename.append("\\secforc");

	std::ifstream inFile;
	inFile.open(filename);
	if (!inFile) {
		std::cerr << "Unable to open file " << filename << std::endl;
		return false;
	}

	double d_max = 0;
	int flag = -1;
	int cnt = 0; 
	int cntMax = 100;
	while (!inFile.eof()) {
		std::string line;
		getline(inFile, line);
		// skip to data
		if (line == "{END LEGEND}") {
			flag = 0;
			// skip six lines
			for (int i = 0; i < 6; i++)
				getline(inFile, line);
			continue;
		}
		// read line
		if (flag == 0) {
			std::istringstream iss(line);
			double d;
			for (int i = 0; i < 6; i++) {
				iss >> d;
			}
			if (d > d_max) {
				cnt = 0;
				d_max = d;
			}
			if (cnt++ > cntMax) {
				data = d_max;
				return true;
			}

			// skip three lines
			for (int i = 0; i < 3; i++)
				getline(inFile, line);
		}
	}
	return true;
}

bool CFEMeshIO::readFEResultSECFORC(std::string fpath, std::vector<double>& data_list, std::vector<double> &timestep)
{
	std::string filename = fpath;
	filename.append("\\secforc");

	std::ifstream inFile;
	inFile.open(filename);
	if (!inFile) {
		std::cerr << "Unable to open file " << filename << std::endl;
		return false;
	}

	//double d_max = 0;
	int flag = -1;
	//int cnt = 0;
	//int cntMax = 5;
	while (!inFile.eof()) {
		std::string line;
		getline(inFile, line);
		if (inFile.fail()) break;
		// skip to data
		if (line == "{END LEGEND}") {
			flag = 0;
			// skip six lines
			for (int i = 0; i < 6; i++)
				getline(inFile, line);
			continue;
		}
		// read line
		if (flag == 0) {
			std::istringstream iss(line);
			double d;
			for (int i = 0; i < 6; i++) {
				iss >> d;
				if (i == 1)
				{
					timestep.push_back(d);
				}
			}
			data_list.push_back(d);
			// skip three lines
			for (int i = 0; i < 3; i++)
				getline(inFile, line);
		}
	}
	return true;
}

bool CFEMeshIO::readFEResultGLSTATIE(std::string fpath, double & data)
{
	std::string filename = fpath;
	filename.append("\\glstat");

	std::ifstream inFile;
	inFile.open(filename);
	if (!inFile) {
		std::cerr << "Unable to open file " << filename << std::endl;
		return false;
	}

	// goto the 8th line
	std::string line;
	for (int i = 0; i < 9; i++)
	{
		getline(inFile, line);
		//std::cout << line << std::endl;
	}

	//double d_max = 0;
	while (!inFile.eof()) {
		std::string line;
		getline(inFile, line);
		if (inFile.fail()) break;
		// get double value
		std::istringstream iss(line);
		std::string temp;
		double d = 0.0;
		while (!iss.eof()) {
			iss >> temp;
			if (std::stringstream(temp) >> d) {
				data = d;
			}
			temp = "";
		}

		// skip 22 lines
		for (int i = 0; i < 22; i++) {
			if (getline(inFile, line).eof()) {
				break;
			}
		}
	}
	return true;
}

bool CFEMeshIO::readFEResultGLSTATIE(std::string fpath, std::vector<double>& data_list)
{
	std::string filename = fpath;
	filename.append("\\glstat");

	std::ifstream inFile;
	inFile.open(filename);
	if (!inFile) {
		std::cerr << "Unable to open file " << filename << std::endl;
		return false;
	}

	// goto the 8th line
	std::string line;
	for (int i = 0; i < 9; i++)
	{
		getline(inFile, line);
		//std::cout << line << std::endl;
	}

	//double d_max = 0;
	while (!inFile.eof()) {
		std::string line;
		getline(inFile, line);

		// get double value
		std::istringstream iss(line);
		std::string temp;
		double d = 0.0;
		while (!iss.eof()) {
			iss >> temp;
			if (std::stringstream(temp) >> d) {
				data_list.push_back(d);
			}
			temp = "";
		}

		// skip 22 lines
		for (int i = 0; i < 22; i++) {
			if (getline(inFile, line).eof()) {
				break;
			}
		}
	}
	return true;
}
