#include <iostream>
#include <fstream>
#include <string>

#include "../DataColle/data_io.h"
#include "../DataColle/curve_object.h"
#include "../DataColle/mesh_object.h"
#include "../DeformationAlg/skeleton_object.h"
#include "../DeformationAlg/sst_object.h"
#include "file_io.h"

//using namespace std;

static std::string CONFIG_FILE{ "fDataInput/config_switch.txt" };   //输入文件名称
static std::string DEF_CROSS_SECTIONS_FILE{ "fDataInput/def_cross_sections.txt" };   //输入文件名称

//enum SwitchType {sstExit, sstGlobal, sstLocal, sstPending};

static SwitchType SWITCH(SwitchType::sstExit); // 全局开关

int main(int argc ,char* argv)
{
	//////////////////////////////////////////////////////////////////////////

	// read mesh
	CMeshObject * p_mesh_src_obj = new CMeshObject;
	if (!CDataIO::ReadMesh("fDataInput/in_mesh_dense.stl", *p_mesh_src_obj))
		return 0;

	// read skeleton
	CCurveObject * p_skel_src_obj = new CCurveObject;
	if (!CDataIO::ReadCurve("fDataInput/in_skeleton.txt", *p_skel_src_obj))
		return 0;
	//////////////////////////////////////////////////////////////////////////	

	CSkeleton skeleton(p_skel_src_obj->GetCurve());
	CSstObject* sst_obj = new CSstObject(skeleton, p_mesh_src_obj->GetMesh());

	//sst_obj->SetParameters(
	//	SRC_SKELETON_FILE,
	//	DEF_SKELETON_FILE,
	//	SRC_CROSS_SECTIONS_FILE,
	//	DEF_CROSS_SECTIONS_FILE,
	//	IN_MESH_FILE,
	//	IN_EMB_MESH_FILE,
	//	OUT_MESH_FILE,
	//	OUT_EMB_MESH_FILE
	//);

	// ENCODING
	sst_obj->Encode();

	// SET SID_LIST
	std::cout << "Sid list: " << "(# skeletal pts " << skeleton.GetAccumArcLength().size() << ")" << std::endl;
	std::vector<int> sid_list;
	int cnt = 5;
	while (cnt < skeleton.GetAccumArcLength().size())
	{
		std::cout << cnt << " "; // std::endl;
		sid_list.push_back(cnt);
		cnt += 5;
	}
	std::cout << std::endl;

	sst_obj->InitCrossSections(sid_list);
	sst_obj->PrintSSTandMesh();

	// LOOP TO CATCH ANY CALL
	int cntIter = 0, maxIter = 10000;
	while (++cntIter < 3)
	{
		std::cout << cntIter << std::endl;

		// step 1: read the config file and update the two global switchs
		SWITCH = CFileIO::read_config_file(CONFIG_FILE);

		// exit when all switches are off
		if (SWITCH == SwitchType::sstExit)
		{
			std::cerr << "EXIT" << std::endl;
			return 1;
		}
		// GO ON GLOBAL, USING THE PRE-FACTORIZATION
		else if (SWITCH == SwitchType::sstPending) {
			// pending
			//std::cout << "pending" << std::endl;
		}
		else if (SWITCH == SwitchType::sstGlobal)
		{
			std::cout << "Continue the Global deformation...(UNDER CONSTRUCTION)" << std::endl;
			CSkeleton def_skeleton;
			CFileIO::read_skeleton("fDataRepo/def_skeleton.txt", def_skeleton);
			sst_obj->SetDefSkeleton(def_skeleton);
			sst_obj->GlobalDeformSetup();
			sst_obj->GlobalDeformSolve();
			//
			std::cout << "Global deformation done" << std::endl;
			bool is_deformed = true;
			sst_obj->PrintMesh("fDataRepo/dst_gmesh.txt", "fDataRepo/dst_emb_gmesh.txt", is_deformed);
		}
		// GO ON LOCAL, USING THE PRE-FACTORIZATION
		else if (SWITCH == SwitchType::sstLocal)
		{
			std::cout << "Continue the Local deformation...(UNDER CONSTRUCTION)" << std::endl;
			std::vector<int> def_csid_list;
			CFileIO::read_local_def_cs_id_from_config_file(CONFIG_FILE, def_csid_list);
			std::string fname_dst = "fDataRepo/CSdata/dst_cross_sections";
			std::string fname_src = "fDataRepo/CSdata/src_cross_sections";
			for (int i = 0; i < def_csid_list.size(); i++)
			{
				std::ostringstream id_str;
				id_str << def_csid_list[i];
				// get src_cs
				std::string fname_cs_src = fname_src;
				fname_cs_src.append(id_str.str());
				fname_cs_src.append(".txt");
				std::vector<COpenMeshT::Point> src_cs_pts;
				CFileIO::read_cross_section_pts(fname_cs_src, src_cs_pts);
				// get dst_cs
				std::string fname_cs_dst = fname_dst;
				fname_cs_dst.append(id_str.str());
				fname_cs_dst.append(".txt");
				std::vector<COpenMeshT::Point> dst_cs_pts;
				CFileIO::read_cross_section_pts(fname_cs_dst, dst_cs_pts);
				// if not true, we need to update the original cs as well
				if (src_cs_pts.size() == dst_cs_pts.size()) {
					sst_obj->RenewCSInfo(def_csid_list[i], src_cs_pts);
					sst_obj->SetDefCSInfo(def_csid_list[i], dst_cs_pts);
				}
				else {
					std::cerr << "please make sure CS" << def_csid_list[i] << "'s info is correct" << std::endl;
				}
			}
			sst_obj->LocalDeformSetup();
			sst_obj->LocalDeformSolve();
			//
			std::cout << "Local deformation done" << std::endl;
			bool is_deformed = true;
			sst_obj->PrintMesh("fDataRepo/dst_lmesh.txt", "fDataRepo/dst_emb_lmesh.txt", is_deformed);
		}
	}

	return 1;
}