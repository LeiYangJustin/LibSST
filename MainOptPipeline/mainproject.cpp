#include"CrossSection.h"
#include"File_OI.h"
#include"General_function.h"
#include <stdio.h>

#include "../DataColle/data_io.h"
#include "../DataColle/curve_object.h"
#include "../DataColle/mesh_object.h"
#include "../DeformationAlg/skeleton_object.h"
#include "../DeformationAlg/sst_object.h"
#include "../GeneralTools/geo_calculator.h"

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

int main()
{
	int cs_spacing=5;//取样截面之间的间隔
	CFileIO::InputPaths in_paths;
	CSstObject* sst_obj = new CSstObject;
	char cCurrentPath[FILENAME_MAX];

	if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
	{
		std::cerr << "No available current dir" << std::endl;
		return 0;
	}

	std::string path_to_gen_mesh = cCurrentPath;
	path_to_gen_mesh.append("\\gen_mesh.stl");//网格模型文件的路径
	

	//读取用户给与的数据文件
	file_io file;
	file.set_basic_data("prac.xml");
	std::vector<Basic_data_of_DesSection> def_cs_list = file.get_def_section_data_list();


	std::cout << "Temporary init..." << std::endl;
	// read mesh
	CMeshObject * p_mesh_src_obj = new CMeshObject;
	if (!CDataIO::ReadMesh("Data\\sf\\Sframe.stl", *p_mesh_src_obj)) {
		std::cerr << "Path to the mesh model cannot be found!\n";
		return 0;
	}
	// read skeleton
	CSkeleton skeleton;
	std::string skel_fname = "Data\\sf\\src_ctrl_skeleton.xml";
	if (!CFileIO::xml_read_skeleton(skel_fname, skeleton)) {
		std::cerr << "Path to the input skeleton cannot be found!\n";
		return 0;
	}
	// visualize skeleton
	// Initialization SST
	sst_obj->SetSkeleton(skeleton);
	sst_obj->SetMesh(p_mesh_src_obj->GetMesh());
	sst_obj->Encode();
	// Get cross-sections
	std::vector<int> sid_list;
	int cnt = 5;
	while (cnt < skeleton.GetAccumArcLength().size()) {
		sid_list.push_back(cnt);
		cnt += cs_spacing;
	}
	sst_obj->InitCrossSections(sid_list);
	

	std::cout << "Local deformation starts..." << std::endl;
	// reset cross-sections
	sst_obj->ResetCrossSectionsAfterDeformation();
	//将控制点的变形转换为所有点的变形
	std::vector<int> def_csid_list;
	for (int cs_id = 0; cs_id < def_cs_list.size(); cs_id++)
	{
		CrossSection cs;
		Basic_data_of_DesSection sec = def_cs_list[cs_id];
		cs.set_the_number_of_edgepts(sec.pts_num_on_single_edge);
		cs.set_ctrl_pid_list_and_def_list(sec.ctrl_pid_list_, sec.ctrl_pt_def_list_);
		cs.set_rectangular_cross_section(sec.length_, sec.width_);
		cs.deform();
		std::vector<Point2> dst_cs_pt = cs.get_def_pt_list();
		std::vector<Point2> src_cs_pt = cs.get_src_pt_list();
		std::vector<COpenMeshT::Point> dst_cs_pts;
		std::vector<COpenMeshT::Point> src_cs_pts;

		for (int i = 0; i < dst_cs_pt.size()-1; i++)
		{
			double x = def_cs_list[cs_id].x_coordinate_;
			dst_cs_pts.push_back(COpenMeshT::Point(x, dst_cs_pt[i].y_, dst_cs_pt[i].z_));
			src_cs_pts.push_back(COpenMeshT::Point(x, src_cs_pt[i].y_, src_cs_pt[i].z_));
		}

		int sid = def_cs_list[cs_id].spid_;
		def_csid_list.push_back(sid);//变形截面的spid
		if (!sst_obj->RenewCSInfo(sid, src_cs_pts, dst_cs_pts)) {
			std::cerr << "please make sure CS" << sid << "'s info is correct" << std::endl;
		}
	}


	// deform
	sst_obj->SetDefCSList(def_csid_list);
	sst_obj->LocalDeformSetup();
	sst_obj->LocalDeformSolve();

	COpenMeshT def_mesh, src_mesh;
	if (sst_obj->OutputDefMesh(def_mesh))
	{
		/*src_mesh = sst_obj->GetMesh();
		CFileIO::write_mesh_to_stl("Data\\sf\\Sframe.stl", src_mesh, false);*/
		CFileIO::write_mesh_to_stl(path_to_gen_mesh, def_mesh, false);
		std::cerr << "Gen_Mesh output to path: " << path_to_gen_mesh << std::endl;
		std::cerr << "Local deformation done" << std::endl;


		//CFileIO::write_mesh(path_to_src_mesh_txt, src_mesh);
		//CFileIO::write_mesh(path_to_gen_mesh_txt, def_mesh);
	}
	else {
		std::cerr << "No deformation result is available" << std::endl;
		return 0;
	}
	
	delete p_mesh_src_obj;
	delete sst_obj;
	system("pause");
	return 0;

}