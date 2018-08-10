#ifndef C_FILE_READER_H
#define C_FILE_READER_H
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "../DataColle/mesh_object.h"
#include "../DeformationAlg/skeleton_object.h"
#include "../DeformationAlg/sst_object.h"

enum SwitchType { sstInit, sstUpdateCS, sstGlobal, sstLocal, sstPending, sstExit};

class CFileIO
{
public:
	struct InputPaths
	{
		std::string path_to_config_in;
		std::string path_to_mesh_in;
		std::string path_to_src;
		std::string path_to_dst;
		std::string path_to_outfolder;
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

	static bool xml_read_skeleton(std::string fname, CSkeleton &s);
	static bool xml_write_skeleton_polyline(std::string fname, const CSkeleton &s);
	static bool xml_write_skeleton_bezier(std::string fname, const CSkeleton &s);

	static bool xml_read_a_section(std::string fname, const int csid, CCrossSection &cs);
	//static bool xml_write_a_section(std::string fname, const int csid, const CCrossSection &cs);
	
	static bool xml_read_sections(std::string fname, std::vector<CCrossSection> &cs_list);
	static bool xml_write_sections(std::string fname, const std::vector<CCrossSection> &cs_list);

	//static SwitchType xml_read_config_file(
	//	const std::string fname, CFileIO::InputPaths &in_paths);

	static SwitchType ini_read_config_file(
		const std::string fname, CFileIO::InputPaths &in_paths);
	static bool ini_write_config_file_to_pending(const std::string fname);

	static bool write_mesh(std::string fname, COpenMeshT &mesh, bool is_emb = false);
	static bool write_mesh_to_stl(std::string fname, COpenMeshT &mesh, bool is_emb = false);
	static bool write_mesh_to_obj(std::string fname, COpenMeshT &mesh, bool is_emb = false);

	static bool output_SST(CSstObject* sst, CFileIO::OutputPaths out_paths);
	static bool output_CS(CSstObject* sst, std::string fname);
	static void set_output_path(std::string path, CFileIO::OutputPaths &output_paths);

};

#endif