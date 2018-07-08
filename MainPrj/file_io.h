#ifndef C_FILE_READER_H
#define C_FILE_READER_H
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "../DataColle/mesh_object.h"
#include "../DeformationAlg/skeleton_object.h"
#include "../DeformationAlg/sst_object.h"

enum SwitchType { sstInit, sstUpdate, sstGlobal, sstLocal, sstPending, sstExit};

class CFileIO
{
public:
	struct InputPaths
	{
		std::string path_to_config_in;
		std::string path_to_mesh_in;
		std::string path_to_src;
		std::string path_to_dst;
	};

	struct OutputPaths
	{
		std::string path_to_folder;
		std::string path_to_mesh_in;
		std::string path_to_mesh_emb;
		std::string path_to_mesh_out;
		std::string path_to_src_skel;
		std::string path_to_def_skel;
		std::string path_to_src_cs;
		std::string path_to_def_cs;
	};

public:

	static SwitchType read_config_file(std::string fname);

	static bool read_local_def_cs_id_from_config_file(std::string fname, std::vector<int> &cs_ids);
	static bool read_cross_section_pts(std::string fname, std::vector<COpenMeshT::Point> &cs_pts);
	static bool read_skeleton(std::string fname, CSkeleton &s);


	static bool xml_read_skeleton(std::string fname, CSkeleton &s);
	static bool xml_write_skeleton(std::string fname, const CSkeleton &s);

	static bool xml_read_a_section(std::string fname, const int csid, CCrossSection &cs);
	static bool xml_write_a_section(std::string fname, const int csid, const CCrossSection &cs);
	static bool xml_write_sections(std::string fname, const std::vector<CCrossSection> &cs_list);

	static SwitchType xml_read_config_file(
		const std::string fname, CFileIO::InputPaths in_paths);

	static bool write_mesh(std::string fname, COpenMeshT &mesh, bool is_emb);
	static bool print_SST(CSstObject* sst);
	static void set_output_path(std::string path, CFileIO::OutputPaths o_paths);

};

#endif