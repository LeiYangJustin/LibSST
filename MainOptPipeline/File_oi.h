#ifndef FILE_OI_H
#define FILE_OI_H
#include<fstream>
#include<string>
#include<vector>
#include<sstream>
#include "pugixml\pugixml.hpp"
#include "pugixml\pugiconfig.hpp"
#include "ini_parser\ini.hpp"
#include"General_function.h"
#include "../DeformationAlg/skeleton_object.h"
#include "../DeformationAlg/sst_object.h"

enum SwitchType { sstInit, sstUpdateCS, sstGlobal, sstLocal, sstPending, sstExit };

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

class file_io//
{
private:
	
	std::vector<Basic_data_of_DesSection> def_section_list_;//截面的基本信息
	int times_=0;//有多少个模型变体
	bool xml_is_empty_ = true;

	void Read_section_data_from_xml(std::string Dessection_fname);//读取并设置截面基本信息
	void Write_xml(Basic_data_of_DesSection def_section, std::vector<Point2> def_pt_list, std::string outfname);
public:
	file_io() {};
	~file_io() {};
	void set_basic_data(std::string Dessection_fname);
	void get_outxml(Basic_data_of_DesSection def_section, std::vector<Point2> def_pt_list, std::string outfname);
	std::vector<Basic_data_of_DesSection> get_def_section_data_list();
	int get_times()
	{
		return times_;
	}
	
};


#endif // !FILE_OI_H

