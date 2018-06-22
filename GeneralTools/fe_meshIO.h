#ifndef C_FE_MESH_IO_H
#define C_FE_MESH_IO_H

#include "general_tool_prereq.h"
#include "../DataColle/mesh_object.h"

class GENERAL_TOOLS_CLASS CFEMeshIO
{
public:

	// read key file and store data in V and F
	static bool readKeyFile(std::string filename, 
		std::map<int, std::vector<double>> &Vmap, 
		std::vector<std::vector<int>> &F);

	// read key file and store data in COpenMeshT class
	static bool readKeyFile(std::string filename, COpenMeshT & m);

	//// read key file and store data in COpenMeshT class
	//static bool readKeyFile(std::string filename, MyOpenMeshP & pm);
	
	//static bool writeKeyFile(std::string filename, std::vector<std::vector<double>> V, std::vector<std::vector<int>> F);

	static bool updateNodesInKeyFile(std::string filename, 
		std::string filename_src, 
		std::map<int, std::vector<double>> &V);



	//
	static bool readFEResultSECFORC(std::string fpath, double &data);
	static bool readFEResultSECFORC(std::string fpath, std::vector<double> &data_list, std::vector<double> &timestep);
	
	static bool readFEResultGLSTATIE(std::string fpath, double &data);
	static bool readFEResultGLSTATIE(std::string fpath, std::vector<double> &data_list);

};


#endif // !C_FE_MESH_IO_H



