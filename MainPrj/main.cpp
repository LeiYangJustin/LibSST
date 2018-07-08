#include <iostream>
#include <fstream>
#include <string>

#include "../DataColle/data_io.h"
#include "../DataColle/curve_object.h"
#include "../DataColle/mesh_object.h"
#include "../DeformationAlg/skeleton_object.h"
#include "../DeformationAlg/sst_object.h"
#include "../GeneralTools/geo_calculator.h"
#include "file_io.h"

int main(int argc ,char* argv)
{
	argv = "E:\\Research\\SSTsystem\\LibSST\\MainPrj\\Data\\SolveController_global_example.xml";

	// init
	CSstObject* sst_obj = new CSstObject;

	

	// parameters
	// 1 - cs spacing
	int cs_spacing = 5;
	SwitchType solver_state(SwitchType::sstPending); 

	//CSkeleton s;
	//CFileIO::xml_read_skeleton("E:\\Research\\SSTsystem\\LibSST\\MainPrj\\xml\\Skeleton.xml", s);
	//CFileIO::xml_write_skeleton("E:\\Research\\SSTsystem\\LibSST\\MainPrj\\xml\\defSkeleton.xml", s);

	//CCrossSection cs;
	//CFileIO::xml_read_a_section("E:\\Research\\SSTsystem\\LibSST\\MainPrj\\xml\\Section.xml", 2, cs);
	//CFileIO::xml_write_a_section("E:\\Research\\SSTsystem\\LibSST\\MainPrj\\xml\\defSection.xml", 2, cs);


	//////////////////////////////////////////////////////////
	//std::vector<COpenMeshT::Point> pts;
	//std::string fname = "E:\\Research\\SSTsystem\\LibSST\\MainPrj\\fDataRepo\\CSData\\src_cross_sections18.txt";
	//CFileIO::read_cross_section_pts(fname, pts);
	//CGeoCalculator::reconstruct_curve_from_pointset(pts, 0.05);
	//////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////
	//// read mesh
	//CMeshObject * p_mesh_src_obj = new CMeshObject;
	//if (!CDataIO::ReadMesh("Data/Sframe.stl", *p_mesh_src_obj)) {
	//	std::cerr << "Path to the mesh model cannot be found!\n";
	//	return 0;
	//}
	//CCurveObject * p_skel_src_obj = new CCurveObject;
	//if (!CDataIO::ReadCurve("Data/in_skeleton.txt", *p_skel_src_obj)) {
	//	std::cerr << "Path to the input skeleton cannot be found!\n";
	//	return 0;
	//}
	//// Initialization SST
	//CSkeleton skeleton(p_skel_src_obj->GetCurve());
	//sst_obj->SetSkeleton(skeleton);
	//sst_obj->SetMesh(p_mesh_src_obj->GetMesh());
	//sst_obj->Encode();
	//// Get cross-sections
	//std::vector<int> sid_list;
	//int cnt = 5;
	//while (cnt < skeleton.GetAccumArcLength().size()) {
	//	sid_list.push_back(cnt);
	//	cnt += cs_spacing;
	//}
	//sst_obj->InitCrossSections(sid_list);


	//// Output
	////sst_obj->PrintSSTandMesh();
	//CFileIO::print_SST(sst_obj);
	//////////////////////////////////////////////////////////

	//return 0;

	// LOOP TO CATCH ANY CALL
	int cntIter = 0, maxIter = 10000;
	while (++cntIter < 4)
	{
		std::cout << cntIter << std::endl;

		// step 1: read the config file and update the two global switchs
		CFileIO::InputPaths in_paths;
		solver_state = CFileIO::xml_read_config_file(argv, in_paths);

		// exit when all switches are off
		switch (solver_state)
		{
		case SwitchType::sstExit: 
		{
			std::cerr << "EXIT" << std::endl;
			return 1;
		}
		case SwitchType::sstUpdate:
		{
			std::cout << "Update is under construction" << std::endl;
			//std::cout << "SST has been updated" << std::endl;

			// in this part, we shall allow users to do any modification to the selected cross-sections.

			//// get sids
			//std::vector<int> sid_list;
			//sst_obj->InsertCrossSections(sid_list);
			sst_obj->Update();
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
			sst_obj->PrintSSTandMesh();

			break;
		}
		case SwitchType::sstGlobal:
		{
			std::cout << "Global deformation..." << std::endl;
			
			// reset cross-sections
			sst_obj->ResetCrossSectionsAfterDeformation();
			
			// input
			CSkeleton def_skeleton;
			if (!CFileIO::xml_read_skeleton(in_paths.path_to_dst, def_skeleton))
			{
				return 0;
			}
			//CFileIO::read_skeleton("fDataRepo/def_skeleton.txt", def_skeleton);
			sst_obj->SetDefSkeleton(def_skeleton);

			// deform
			sst_obj->GlobalDeformSetup();
			sst_obj->GlobalDeformSolve();
			bool is_deformed = true;
			//sst_obj->PrintMesh("fDataRepo/dst_gmesh.txt", "fDataRepo/dst_emb_gmesh.txt", is_deformed);
			////sst_obj->Update(); 
			////
			std::cout << "Global deformation done" << std::endl;
			std::this_thread::sleep_for(std::chrono::seconds(10));
			break;
		}
		case SwitchType::sstLocal:
		{
			std::cout << "Local deformation..." << std::endl;

			// get path

			// reset cross-sections
			sst_obj->ResetCrossSectionsAfterDeformation();

			// read def cs

			// read src cs

			// get def_csid_list

			std::vector<int> def_csid_list;
			//CFileIO::read_local_def_cs_id_from_config_file(CONFIG_FILE, def_csid_list);
			//std::string fname_dst = "fDataRepo/CSdata/dst_cross_sections";
			//std::string fname_src = "fDataRepo/CSdata/src_cross_sections";

			//// set up
			//sst_obj->SetDefCSList(def_csid_list);
			//for (int i = 0; i < def_csid_list.size(); i++)
			//{
			//	std::ostringstream id_str;
			//	id_str << def_csid_list[i];
			//	// get src_cs
			//	std::string fname_cs_src = fname_src;
			//	fname_cs_src.append(id_str.str());
			//	fname_cs_src.append(".txt");
			//	std::vector<COpenMeshT::Point> src_cs_pts;
			//	CFileIO::read_cross_section_pts(fname_cs_src, src_cs_pts);
			//	// get dst_cs
			//	std::string fname_cs_dst = fname_dst;
			//	fname_cs_dst.append(id_str.str());
			//	fname_cs_dst.append(".txt");
			//	std::vector<COpenMeshT::Point> dst_cs_pts;
			//	CFileIO::read_cross_section_pts(fname_cs_dst, dst_cs_pts);
			//	if (!sst_obj->RenewCSInfo(def_csid_list[i], src_cs_pts, dst_cs_pts)) {
			//		std::cerr << "please make sure CS" << def_csid_list[i] << "'s info is correct" << std::endl;
			//	}
			//}

			// deform
			sst_obj->LocalDeformSetup();
			sst_obj->LocalDeformSolve();
			bool is_deformed = true;
			sst_obj->PrintMesh("fDataRepo/dst_lmesh.txt", "fDataRepo/dst_emb_lmesh.txt", is_deformed);

			// 
			std::cout << "Local deformation done" << std::endl;
			std::this_thread::sleep_for(std::chrono::seconds(10));
			break;
		}
		case SwitchType::sstPending:
		{
			std::cout << "Pending for 10 seconds..." << std::endl;
			std::this_thread::sleep_for(std::chrono::seconds(10));
			break;
		}
		default : {
			std::cout << "Pending for 10 seconds..." << std::endl;
			std::this_thread::sleep_for(std::chrono::seconds(10));
			break;
		}

		} // end switch
	}

	return 1;
}