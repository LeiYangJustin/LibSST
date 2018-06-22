#ifndef C_FILE_READER_H
#define C_FILE_READER_H
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "../DataColle/mesh_object.h"
#include "../DeformationAlg/skeleton_object.h"

enum SwitchType { sstExit, sstGlobal, sstLocal, sstPending};

class CFileIO
{
public:
	static SwitchType read_config_file(std::string fname);
	static bool read_local_def_cs_id_from_config_file(std::string fname, std::vector<int> &cs_ids);

	//static bool write_cross_sections_file(std::string fname);
	static bool read_cross_sections_file(std::string fname, std::vector<CCrossSection> &cs_list);
	static bool read_cross_section_pts(std::string fname, std::vector<COpenMeshT::Point> &cs_pts);

	//static bool write_skeleton(std::string fname, CSkeleton s);
	static bool read_skeleton(std::string fname, CSkeleton &s);
};

#endif