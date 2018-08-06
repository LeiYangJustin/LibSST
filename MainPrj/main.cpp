#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>  /* defines FILENAME_MAX */

#include "../DataColle/data_io.h"
#include "../DataColle/curve_object.h"
#include "../DataColle/mesh_object.h"
#include "../DeformationAlg/skeleton_object.h"
#include "../DeformationAlg/sst_object.h"
#include "../GeneralTools/geo_calculator.h"
#include "file_io.h"

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

int main(int argc ,char** argv)
{
	
	std::cout << "argc: " << argc << std::endl;

	char cCurrentPath[FILENAME_MAX];

	if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
	{
		std::cerr << "No available current dir" << std::endl;
		return 0;
	}

	// parameters
	// argv[0]: the exe file
	// argv[1]: configure file (.ini)
	// argv[2]: max iteration to stop
	// argv[3]: spacing
	int cs_spacing, maxIter;
	std::string path_to_cfg;
	bool is_testing = false;
	SwitchType solver_state = SwitchType::sstPending;
	if (argc == 1)
	{
		std::cerr << "This is C++ testing mode with argc == 1" << std::endl;
		path_to_cfg = "Data\\SolverInit.ini";
		maxIter = 1;
		cs_spacing = 5;
		is_testing = true;
	}
	else if (argc == 2)
	{
		path_to_cfg = argv[1];
		maxIter = 10000;
		cs_spacing = 5;
	}
	else if (argc == 3) {
		path_to_cfg = argv[1];
		maxIter = atoi(argv[2]);
		cs_spacing = 5;
	}
	else {
		return 0;
	}
	CSstObject* sst_obj = new CSstObject;
	// LOOP TO CATCH ANY CALL
	while (solver_state != SwitchType::sstExit)
	{
		std::cout << "\n... " << std::endl;
		// step 1: read the config file and update the two global switchs
		CFileIO::InputPaths in_paths;
		solver_state = CFileIO::ini_read_config_file(path_to_cfg, in_paths);
		std::string path_to_gen_mesh = cCurrentPath;
		path_to_gen_mesh.append("\\gen_mesh.stl");

		// exit when all switches are off
		switch (solver_state)
		{
		// in this part, we shall allow users to do any modification to the selected cross-sections.
		case SwitchType::sstUpdateCS:
		{
			std::cout << "Update is under construction" << std::endl;
			if (is_testing)
			{
				// read mesh
				CMeshObject * p_mesh_src_obj = new CMeshObject;
				if (!CDataIO::ReadMesh(in_paths.path_to_mesh_in, *p_mesh_src_obj)) {
					std::cerr << "Path to the mesh model cannot be found!\n";
					return 0;
				}
				// read skeleton
				CSkeleton skeleton;
				if (!CFileIO::xml_read_skeleton(in_paths.path_to_src, skeleton)) {
					std::cerr << "Path to the input skeleton cannot be found!\n";
					return 0;
				}
				// visualize skeleton
				std::cout << skeleton.GetCtrlPts().size() << std::endl;
				//CFileIO::xml_write_skeleton_bezier(in_paths.path_to_src, skeleton);
				// Initialization SST
				sst_obj->SetSkeleton(skeleton);
				sst_obj->SetMesh(p_mesh_src_obj->GetMesh());
				sst_obj->Encode();
			}
			std::vector<int> add_spids;
			add_spids.push_back(14);
			// csid: 8, 10
			sst_obj->InsertCrossSections(add_spids);
			//// get sids
			//std::vector<int> delete_spids;
			//sst_obj->DeleteCrossSections(delete_spids);
			
			sst_obj->Update();
			return 1;
			break;
		}
		case SwitchType::sstInit:
		{
			// read mesh
			CMeshObject * p_mesh_src_obj = new CMeshObject;
			if (!CDataIO::ReadMesh(in_paths.path_to_mesh_in, *p_mesh_src_obj)) {
				std::cerr << "Path to the mesh model cannot be found!\n";
				return 0;
			}
			// read skeleton
			CSkeleton skeleton;
			if (!CFileIO::xml_read_skeleton(in_paths.path_to_src, skeleton)) {
				std::cerr << "Path to the input skeleton cannot be found!\n";
				return 0;
			}
			// visualize skeleton
			std::cout << skeleton.GetCtrlPts().size() << std::endl;
			if (skeleton.GetCtrlPts().size()) {
				CFileIO::xml_write_skeleton_bezier(in_paths.path_to_src, skeleton);
			}
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
			// output
			std::cout << "Initialization done" << std::endl;
			std::string path_to_gen_cs = cCurrentPath;
			path_to_gen_cs.append("\\src_sections.xml");
			CFileIO::output_CS(sst_obj, path_to_gen_cs);

			CFileIO::ini_write_config_file_to_pending(path_to_cfg);
			
			break;
		}
		case SwitchType::sstGlobal:
		{
			if (is_testing) {
				std::cerr << "Temporary init..." << std::endl;
				// read mesh
				CMeshObject * p_mesh_src_obj = new CMeshObject;
				if (!CDataIO::ReadMesh(in_paths.path_to_mesh_in, *p_mesh_src_obj)) {
					std::cerr << "Path to the mesh model cannot be found!\n";
					return 0;
				}
				// read skeleton
				CSkeleton skeleton;
				if (!CFileIO::xml_read_skeleton(in_paths.path_to_src, skeleton)) {
					std::cerr << "Path to the input skeleton cannot be found!\n";
					return 0;
				}
				// visualize skeleton
				//std::cout << skeleton.GetCtrlPts().size() << std::endl;
				CFileIO::xml_write_skeleton_bezier(in_paths.path_to_src, skeleton);
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
			}
			if (!sst_obj->IsEncode())
			{
				std::cerr << "SST object is yet to initialize!" << std::endl;
				std::this_thread::sleep_for(std::chrono::seconds(10));
				return 0;
			}

			//
			std::cout << "Global deformation starts..." << std::endl;
			// reset cross-sections
			sst_obj->ResetCrossSectionsAfterDeformation();
			// input
			CSkeleton def_skeleton;
			if (!CFileIO::xml_read_skeleton(in_paths.path_to_dst, def_skeleton)) {
				return 0;
			}
			// visualize def_skeleton
			 CFileIO::xml_write_skeleton_bezier(in_paths.path_to_dst, def_skeleton);
			// set deformed skeleton
			sst_obj->SetDefSkeleton(def_skeleton);
			// deform
			sst_obj->GlobalDeformSetup();
			sst_obj->GlobalDeformSolve();
			// output
			COpenMeshT def_mesh;
			if (sst_obj->OutputDefMesh(def_mesh))
			{
				CFileIO::write_mesh_to_stl(path_to_gen_mesh, def_mesh);
				std::cout << "Global deformation done" << std::endl;
				std::cout << "Gen_Mesh output to path: " << path_to_gen_mesh << std::endl;
				if (is_testing)
					return 0;
				else
					CFileIO::ini_write_config_file_to_pending(path_to_cfg);
			}
			else {
				std::cerr << "No deformation result is available" << std::endl;
				return 0;
			}
			break;
		}
		case SwitchType::sstLocal:
		{
			if (is_testing) {
				std::cout << "Temporary init..." << std::endl;
				// read mesh
				CMeshObject * p_mesh_src_obj = new CMeshObject;
				if (!CDataIO::ReadMesh(in_paths.path_to_mesh_in, *p_mesh_src_obj)) {
					std::cerr << "Path to the mesh model cannot be found!\n";
					return 0;
				}
				// read skeleton
				CSkeleton skeleton;
				std::string skel_fname = "Data\\src_ctrl_skeleton.xml";
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
			}
			if (!sst_obj->IsEncode())
			{
				std::cerr << "SST object is yet to initialize!" << std::endl;
				std::this_thread::sleep_for(std::chrono::seconds(10));
				return 0;
			}

			std::cout << "Local deformation starts..." << std::endl;
			// reset cross-sections
			sst_obj->ResetCrossSectionsAfterDeformation();
			// read def cross-sections
			std::vector<CCrossSection> def_cs_list;
			std::vector<int> def_csid_list;
			std::cout << "load dst sections" << std::endl;
			CFileIO::xml_read_sections(in_paths.path_to_dst, def_cs_list);
			for (int i = 0; i < def_cs_list.size(); i++)
			{
				int sid = def_cs_list[i].GetSid();
				def_csid_list.push_back(sid);
				std::vector<COpenMeshT::Point> dst_cs_pts = def_cs_list[i].GetEmbProfPts();
				CCrossSection cs;
				std::cout << "load src sections" << std::endl;
				CFileIO::xml_read_a_section(in_paths.path_to_src, sid, cs);
				std::vector<COpenMeshT::Point> src_cs_pts = cs.GetEmbProfPts();
				if (!sst_obj->RenewCSInfo(sid, src_cs_pts, dst_cs_pts)) {
					std::cerr << "please make sure CS" << sid << "'s info is correct" << std::endl;
				}
			}
			// deform
			sst_obj->SetDefCSList(def_csid_list);
			sst_obj->LocalDeformSetup();
			sst_obj->LocalDeformSolve();

			COpenMeshT def_mesh;
			if (sst_obj->OutputDefMesh(def_mesh))
			{
				CFileIO::write_mesh_to_stl(path_to_gen_mesh, def_mesh, false);
				std::cerr << "Local deformation done" << std::endl;
				std::cerr << "Gen_Mesh output to path: " << path_to_gen_mesh << std::endl;
				if (is_testing)
					return 0;
				else
					CFileIO::ini_write_config_file_to_pending(path_to_cfg);
			}
			else {
				std::cerr << "No deformation result is available" << std::endl;
				return 0;
			}

			break;
		}
		case SwitchType::sstPending:
		{
			//std::cerr << "Pending..." << std::endl;
			std::this_thread::sleep_for(std::chrono::seconds(2));
			break;
		}
		default : {
			//std::cerr << "Pending..." << std::endl;
			std::this_thread::sleep_for(std::chrono::seconds(2));
			break;
		}

		} // end switch
	}

	std::cerr << "EXIT" << std::endl;
	return 1;

}