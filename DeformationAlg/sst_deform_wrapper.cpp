#include "sst_deform_wrapper.h"
#include <direct.h> 
#include <stdlib.h>  

CSSTDeformWrapper::CSSTDeformWrapper()
{
}


CSSTDeformWrapper::~CSSTDeformWrapper()
{
	delete sst_obj_;
	for (auto liter = p_deformer_list_.begin(); liter != p_deformer_list_.end(); ++liter)
	{
		delete *liter;
	}
}

void CSSTDeformWrapper::SetUp(COpenMeshT triMesh, std::vector<COpenMeshT::Point> skelPts)
{
	skelPts_ = skelPts; // read-in
	triMesh_ = triMesh; // read-in
}

void CSSTDeformWrapper::InitSST()
{
	sst_obj_ = new CSSTparameterization();

	sst_obj_->SetUp(triMesh_, skelPts_);

	sst_obj_->GetEmbeddedMesh(embMesh_);
	sst_obj_->GetPidSkeidMap(pid_skeid_map_);
	sst_obj_->GetPidVidMap(pid_vid_map_);
}

void CSSTDeformWrapper::DeformWithSST(COpenMeshT & out_mesh, 
	std::vector<std::vector<COpenMeshT::Point>> & ori_cs_pts_list,
	std::vector<int> & cs_type_list)
{
	// init
	out_mesh = triMesh_;

	/*--------- CROSS-SECTIONAL CONSTRAINTS ---------------*/
	// Initialize the cross-sections and roi data
	// cross-sections
	const int numCS = 8;
	double id_array[numCS] = { 20, 25, 35, 40, 60, 65, 75, 80 };
	int fix_mask[numCS] = { 0, 1, 1, 0, 0, 2, 2, 0 };
	std::vector<double> id_cs;
	id_cs.insert(id_cs.end(), id_array, id_array + numCS);
	// get constraints
	std::vector<std::vector<COpenMeshT::Point>> cs_pts_list;
	std::vector<std::vector<COpenMeshT::Point>> cs_vecs_list;
	for (int i = 0; i < numCS; i++) 
	{
		// this is a temporary treatment for getting the cross-sectional constraints
		// we shall have a mechanism, either UI-based or other methods, for assigning constraints
		std::vector<COpenMeshT::Point> cs_pts, cs_vecs;
		//readCSdisplacement("deform_tests/Constraint.txt", cs_pts, cs_vecs);
		
		// points
		readCSData("deform_tests/pts_handle.txt", cs_pts);
		double arclength = 0;
		for (int j = 0; j < id_array[i]; j++)
			arclength += (skelPts_[j + 1] - skelPts_[j]).norm();
		for (int j = 0; j < cs_pts.size(); j++)
		{
			cs_pts[j][0] += arclength;
		}
		cs_pts_list.push_back(cs_pts);

		// displacements
		if (fix_mask[i] == 0) {
			for (int j = 0; j < cs_pts.size(); j++)
			{
				COpenMeshT::Point p(0, 0, 0);
				cs_vecs.push_back(p);
			}
		}
		else if (fix_mask[i] == 1)
		{
			readCSData("deform_tests/vecs_handle11.txt", cs_vecs);
		}
		else if (fix_mask[i] == 2) {
			readCSData("deform_tests/vecs_handle12.txt", cs_vecs);
		}
		cs_vecs_list.push_back(cs_vecs);

		//
		std::cout << "# pts: " << cs_pts.size() << ", # vecs: " << cs_vecs.size() << std::endl;
	}
	/*--------- CROSS-SECTIONAL CONSTRAINTS END------------*/


	/*------------- OUTPUT CROSS-SECTIONS -----------------*/
	for (int i = 0; i < numCS; i++) {
		std::vector<COpenMeshT::Point> cs_pts;
		for (int j = 0; j < cs_pts_list[i].size(); j++)
		{
			Vector3d pt;
			Vector3d emb_pt(cs_pts_list[i][j][0], cs_pts_list[i][j][1], cs_pts_list[i][j][2]);
			sst_obj_->DecodePosition(pt, emb_pt, id_cs[i]);
			COpenMeshT::Point p(pt(0), pt(1), pt(2));
			cs_pts.push_back(p);
		}
		ori_cs_pts_list.push_back(cs_pts);
		cs_type_list.push_back(fix_mask[i]);
	}
	/*------------- OUTPUT CROSS-SECTIONS END -------------*/

	/*------------------- GET PTS IN ROI ---------------------*/
	// define the pairs of cross-sections that bound roi's
	int numROI = 6;
	DenseMatrixXd cspair_map(numROI, 2);
	cspair_map << 0, 1,
		1, 2,
		2, 3,
		4, 5,
		5, 6,
		6, 7;
	std::vector<std::vector<int>> rid_keeper_list(numROI);
	for (auto itmap = pid_skeid_map_.begin(); itmap != pid_skeid_map_.end(); ++itmap) 
	{
		for (int i = 0; i < numROI; i++)
		{
			int id_first = cspair_map(i, 0);
			int id_second = cspair_map(i, 1);

			double al_first = 0, al_second = 0;
			for (int j = 0; j < id_array[id_first]; j++)
				al_first += (skelPts_[j + 1] - skelPts_[j]).norm();
			for (int j = 0; j < id_array[id_second]; j++)
				al_second += (skelPts_[j + 1] - skelPts_[j]).norm();

			// 
			if (itmap->second > id_array[id_first] && itmap->second < id_array[id_second])
			{
				rid_keeper_list[i].push_back(itmap->first);
			}
			else if (itmap->second == id_array[id_first])
			{
				auto p = embMesh_.point(embMesh_.vertex_handle(itmap->first));
				if (p[0] - al_first > 0)
					rid_keeper_list[i].push_back(itmap->first);
			}
			else if (itmap->second == id_array[id_second]) 
			{
				auto p = embMesh_.point(embMesh_.vertex_handle(itmap->first));
				if (p[0] - al_second < 0)
					rid_keeper_list[i].push_back(itmap->first);
			}
		}
	}
	/*------------------- GET PTS IN ROI END -----------------*/

	/*------------------- RBF DEFORMATION --------------------*/
	std::vector<DenseMatrixXd> ori_roi_def_field_list(numROI);
	for (int i = 0; i < numROI; i++) 
	{
		DenseMatrixXd roi_pts;
		roi_pts.resize(rid_keeper_list[i].size(), 3);
		for (int j = 0; j < rid_keeper_list[i].size(); j++) 
		{
			int rid = rid_keeper_list[i][j];
			auto p = embMesh_.point(embMesh_.vertex_handle(pid_vid_map_[rid]));
			roi_pts.row(j) = Vector3d(p[0], p[1], p[2]);
		}
	
		// RBF DEFORM 
		int id_first = cspair_map(i, 0);
		int id_second = cspair_map(i, 1);
		std::vector<double> handlepts, handlevecs;
		for (int j = 0; j < cs_pts_list[id_first].size(); j++)
		{
			handlepts.push_back(cs_pts_list[id_first][j][0]);
			handlepts.push_back(cs_pts_list[id_first][j][1]);
			handlepts.push_back(cs_pts_list[id_first][j][2]);
	
			// handlevecs
			handlevecs.push_back(cs_vecs_list[id_first][j][0]);
			handlevecs.push_back(cs_vecs_list[id_first][j][1]);
			handlevecs.push_back(cs_vecs_list[id_first][j][2]);
		}
		for (int j = 0; j < cs_pts_list[i+1].size(); j++)
		{
			handlepts.push_back(cs_pts_list[id_second][j][0]);
			handlepts.push_back(cs_pts_list[id_second][j][1]);
			handlepts.push_back(cs_pts_list[id_second][j][2]);
	
			// handlevecs
			handlevecs.push_back(cs_vecs_list[id_second][j][0]);
			handlevecs.push_back(cs_vecs_list[id_second][j][1]);
			handlevecs.push_back(cs_vecs_list[id_second][j][2]);
		}
		DenseMatrixXd pts_handles, vecs_handles;
		pts_handles = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>(handlepts.data(), 3, handlepts.size() / 3);
		vecs_handles = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>(handlevecs.data(), 3, handlevecs.size() / 3);
		pts_handles.transposeInPlace();
		vecs_handles.transposeInPlace();
	
		DenseMatrixXd roi_vecs, ori_roi_vecs;
		CMeshDeformation deformer;
		// this precompute the qr decomposition and store it
		deformer.SetUp(pts_handles);
		deformer.Solve(vecs_handles, roi_pts, roi_vecs);
		// decode the field
		sst_obj_->DecodeVectorField(ori_roi_vecs, roi_vecs, rid_keeper_list[i]);
		ori_roi_def_field_list[i] = ori_roi_vecs;
	}
	/*------------------- RBF DEFORMATION END ----------------*/

	/*------------------- PRODUCE OUTPUT MESH ----------------*/
	// add decoded deformation vector fields to original mesh
	std::map<COpenMeshT::VertexHandle, std::pair<COpenMeshT::Point, int>> vh_disp_cnt_map;
	for (int i = 0; i < numROI; i++)
	{
		for (int j = 0; j < rid_keeper_list[i].size(); j++) 
		{
			int rid = rid_keeper_list[i][j];
			Vector3d v = ori_roi_def_field_list[i].row(j);
			auto vh = out_mesh.vertex_handle(pid_vid_map_[rid]);
			if (vh_disp_cnt_map.find(vh) == vh_disp_cnt_map.end()) {
				auto pair = std::make_pair(COpenMeshT::Point(v(0), v(1), v(2)), 1);
				vh_disp_cnt_map.insert(std::make_pair(vh, pair));
			}
			else {
				auto pair = vh_disp_cnt_map[vh];
				pair.first += COpenMeshT::Point(v(0), v(1), v(2));
				pair.second++;
			}
		}
	}
	// average deformation
	for (auto itmap = vh_disp_cnt_map.begin(); itmap != vh_disp_cnt_map.end(); ++itmap)
	{
		COpenMeshT::Point p = out_mesh.point(itmap->first) + itmap->second.first / double(itmap->second.second);
		out_mesh.set_point(itmap->first, p);
	}
	/*------------------- PRODUCE OUTPUT MESH END ------------*/
}

void CSSTDeformWrapper::DeformationLoop(COpenMeshT & out_mesh,
	std::vector<std::vector<COpenMeshT::Point>> & ori_cs_pts_list,
	std::vector<int> & cs_type_list,
	const std::string &fpath,
	double scale)
{
	// Initialize the cross-sections and roi data
	const int numCS = 8;
	double id_array[numCS] = { 20, 25, 35, 40, 60, 65, 75, 80 };
	int fix_mask[numCS] = { 0, 1, 1, 0, 0, 2, 2, 0 };
	std::vector<double> id_cs;
	id_cs.insert(id_cs.end(), id_array, id_array + numCS);
	std::vector<int> fix_mask_list;
	fix_mask_list.insert(fix_mask_list.end(), fix_mask, fix_mask + numCS);

	/*--------- CROSS-SECTIONAL CONSTRAINTS ---------------*/
	std::vector<std::vector<COpenMeshT::Point>> cs_pts_list;
	for (int i = 0; i < numCS; i++)
	{
		// points
		std::vector<COpenMeshT::Point> cs_pts;
		std::string file_pts_handle = fpath;
		file_pts_handle.append("/pts_handle.txt");
		readCSData(file_pts_handle, cs_pts);
		double arclength = 0;
		for (int j = 0; j < id_array[i]; j++)
			arclength += (skelPts_[j + 1] - skelPts_[j]).norm();
		for (int j = 0; j < cs_pts.size(); j++)	{
			cs_pts[j][0] += arclength;
		}
		cs_pts_list.push_back(cs_pts);
	}
	/*--------- CROSS-SECTIONAL CONSTRAINTS END------------*/

	/*------------- OUTPUT CROSS-SECTIONS -----------------*/
	for (int i = 0; i < numCS; i++) {
		std::vector<COpenMeshT::Point> cs_pts;
		for (int j = 0; j < cs_pts_list[i].size(); j++)
		{
			Vector3d pt;
			Vector3d emb_pt(cs_pts_list[i][j][0], cs_pts_list[i][j][1], cs_pts_list[i][j][2]);
			sst_obj_->DecodePosition(pt, emb_pt, id_cs[i]);
			COpenMeshT::Point p(pt(0), pt(1), pt(2));
			cs_pts.push_back(p);
		}
		ori_cs_pts_list.push_back(cs_pts);
		cs_type_list.push_back(fix_mask[i]);
	}
	/*------------- OUTPUT CROSS-SECTIONS END -------------*/

	/*------------------- GET PTS IN ROI ---------------------*/
	// define the pairs of cross-sections that bound roi's
	int numROI = 6;
	DenseMatrixXd cspair_map(numROI, 2);
	cspair_map << 0, 1,
		1, 2,
		2, 3,
		4, 5,
		5, 6,
		6, 7;
	std::vector<std::vector<int>> rid_keeper_list(numROI);
	for (auto itmap = pid_skeid_map_.begin(); itmap != pid_skeid_map_.end(); ++itmap)
	{
		for (int i = 0; i < numROI; i++)
		{
			int id_first = cspair_map(i, 0);
			int id_second = cspair_map(i, 1);

			double al_first = 0, al_second = 0;
			for (int j = 0; j < id_array[id_first]; j++)
				al_first += (skelPts_[j + 1] - skelPts_[j]).norm();
			for (int j = 0; j < id_array[id_second]; j++)
				al_second += (skelPts_[j + 1] - skelPts_[j]).norm();

			// 
			if (itmap->second > id_array[id_first] && itmap->second < id_array[id_second]) {
				rid_keeper_list[i].push_back(itmap->first);
			}
			else if (itmap->second == id_array[id_first]) {
				auto p = embMesh_.point(embMesh_.vertex_handle(itmap->first));
				if (p[0] - al_first > 0)
					rid_keeper_list[i].push_back(itmap->first);
			}
			else if (itmap->second == id_array[id_second]) {
				auto p = embMesh_.point(embMesh_.vertex_handle(itmap->first));
				if (p[0] - al_second < 0)
					rid_keeper_list[i].push_back(itmap->first);
			}
		}
	}
	/*------------------- GET PTS IN ROI END -----------------*/

	/*------------------- DEFORMATION LOOP AND OUTPUT --------*/
	DeformationSetup(cspair_map, cs_pts_list);
	std::cout << "setup success" << std::endl;

	CMeshObject *p_mesh_obj = new CMeshObject();

	int numSamples = 1;
	for (int sample_id = 0; sample_id < numSamples; ) {
		//sample_id++;

		sample_id = 2;
		std::cout << "sample_id = " << sample_id << std::endl;
		
		DeformationSolveAndOutput(sample_id, p_mesh_obj->GetMesh(), cspair_map, rid_keeper_list, fix_mask_list,  fpath);

		// visualization in the window
		out_mesh = p_mesh_obj->GetMesh();
		//
		//p_mesh_obj->Rescale(1.0 / scale);
		////p_mesh_obj->swapXYCoords();
		//std::map<int, std::vector<double>> Vmap;
		//for (auto viter = p_mesh_obj->GetMesh().vertices_begin(); viter != p_mesh_obj->GetMesh().vertices_end(); ++viter) {
		//	std::vector<double> data;
		//	auto p = p_mesh_obj->GetMesh().point(*viter);
		//	data.push_back(p[0]);
		//	data.push_back(p[1]);
		//	data.push_back(p[2]);
		//	int id = p_mesh_obj->GetMesh().data(*viter).get_vlabel();
		//	Vmap[id] = data;
		//}

		////std::string srcFile = "data_deform/test/MeshSF.k";
		//std::string srcFile = fpath;
		//srcFile.append("/MeshSF_rw.k");

		////std::string tgtFile = "data_deform/test/MeshSF_var";
		//std::string sampleDir = "E:\\Research\\SOSST\\sample_folder\\sample";
		//std::ostringstream str_sample_id;
		//str_sample_id << sample_id;
		//sampleDir.append(str_sample_id.str());
		//// create a new folder
		//_mkdir(sampleDir.c_str());
		//std::string tgtFile = sampleDir;
		//tgtFile.append("/MeshSF_rw_v");
		//tgtFile.append(".k");
		//CFEMeshIO::updateNodesInKeyFile(srcFile, tgtFile, Vmap);
	}
	/*------------------- DEFORMATION LOOP AND OUTPUT END -----*/
}

void CSSTDeformWrapper::DeformationSetup(
	const DenseMatrixXd & cspair_map, 
	const std::vector<std::vector<COpenMeshT::Point>> & cs_pts_list)
{
	int numROI = cspair_map.rows();
	for (int i = 0; i < numROI; i++)
	{
		// RBF DEFORM 
		int id_first = cspair_map(i, 0);
		int id_second = cspair_map(i, 1);
		std::vector<double> handlepts, handlevecs;
		for (int j = 0; j < cs_pts_list[id_first].size(); j++)
		{
			handlepts.push_back(cs_pts_list[id_first][j][0]);
			handlepts.push_back(cs_pts_list[id_first][j][1]);
			handlepts.push_back(cs_pts_list[id_first][j][2]);
		}
		for (int j = 0; j < cs_pts_list[id_second].size(); j++)
		{
			handlepts.push_back(cs_pts_list[id_second][j][0]);
			handlepts.push_back(cs_pts_list[id_second][j][1]);
			handlepts.push_back(cs_pts_list[id_second][j][2]);
		}
		DenseMatrixXd pts_handles;
		pts_handles = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>(handlepts.data(), 3, handlepts.size() / 3);
		pts_handles.transposeInPlace();

		CMeshDeformation * deformer = new CMeshDeformation();
		deformer->SetUp(pts_handles);
		p_deformer_list_.push_back(deformer);
	}
	/*------------------- RBF DEFORMATION END ----------------*/
}

void CSSTDeformWrapper::DeformationSolveAndOutput(int sample_id,
	COpenMeshT & out_mesh,
	const DenseMatrixXd & cspair_map,
	const std::vector<std::vector<int>> & rid_keeper_list,
	const std::vector<int> &fix_mask,
	const std::string &fpath)
{
	// init
	std::map<COpenMeshT::VertexHandle, std::pair<COpenMeshT::Point, int>> vh_disp_cnt_map;
	out_mesh = triMesh_;

	int numROI = cspair_map.rows();
	/*------------------- RBF DEFORMATION --------------------*/
	for (int i = 0; i < numROI; i++)
	{
		// ROI
		DenseMatrixXd roi_pts;
		roi_pts.resize(rid_keeper_list[i].size(), 3);
		for (int j = 0; j < rid_keeper_list[i].size(); j++)
		{
			int rid = rid_keeper_list[i][j];
			auto p = embMesh_.point(embMesh_.vertex_handle(pid_vid_map_[rid]));
			roi_pts.row(j) = Vector3d(p[0], p[1], p[2]);
		}

		// RBF DEFORM 
		int id_first = cspair_map(i, 0);
		int id_second = cspair_map(i, 1);
		DenseMatrixXd vecs_h_left, vecs_h_right;
		std::string vecfielname = "/vecs_handle";
		// displacements
		if (fix_mask[id_first] == 1)
		{
			std::ostringstream str_sample_id;
			str_sample_id << sample_id;
			//
			//std::string filename = "data_deform/test/vecs_handle";
			std::string filename = fpath;
			filename.append(vecfielname);
			filename.append(str_sample_id.str());
			filename.append("1.txt");
			//
			std::vector<double> hvecs;
			readCSData(filename, hvecs);
			//
			vecs_h_left = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>(hvecs.data(), 3, hvecs.size() / 3);
			vecs_h_left.transposeInPlace();
		}
		else if (fix_mask[id_first] == 2) {
			std::ostringstream str_sample_id;
			str_sample_id << sample_id;
			//
			std::string filename = fpath;
			filename.append(vecfielname);
			filename.append(str_sample_id.str());
			filename.append("2.txt");
			//
			std::vector<double> hvecs;
			readCSData(filename, hvecs);
			//
			vecs_h_left = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>(hvecs.data(), 3, hvecs.size() / 3);
			vecs_h_left.transposeInPlace();
		}
		else {
			vecs_h_left.resize(0, 0);
		}

		if (fix_mask[id_second] == 1)
		{
			std::ostringstream str_sample_id;
			str_sample_id << sample_id;
			//
			std::string filename = fpath;
			filename.append(vecfielname);
			filename.append(str_sample_id.str());
			filename.append("1.txt");
			//
			std::vector<double> hvecs;
			readCSData(filename, hvecs);
			//
			vecs_h_right = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>(hvecs.data(), 3, hvecs.size() / 3);
			vecs_h_right.transposeInPlace();
		}
		else if (fix_mask[id_second] == 2) {
			std::ostringstream str_sample_id;
			str_sample_id << sample_id;
			//
			std::string filename = fpath;
			filename.append(vecfielname);
			filename.append(str_sample_id.str());
			filename.append("2.txt");
			//
			std::vector<double> hvecs;
			readCSData(filename, hvecs);
			vecs_h_right = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>(hvecs.data(), 3, hvecs.size() / 3);
			vecs_h_right.transposeInPlace();
		}
		else {
			vecs_h_right.resize(0, 0);
		}

		// solve
		DenseMatrixXd roi_vecs, ori_roi_vecs;
		CMeshDeformation * deformer = p_deformer_list_[i];
		deformer->Solve(vecs_h_left, vecs_h_right, roi_pts, roi_vecs);
		
		// decode the field
		sst_obj_->DecodeVectorField(ori_roi_vecs, roi_vecs, rid_keeper_list[i]);
		for (int j = 0; j < rid_keeper_list[i].size(); j++)
		{
			int rid = rid_keeper_list[i][j];
			Vector3d v = ori_roi_vecs.row(j);
			auto vh = out_mesh.vertex_handle(pid_vid_map_[rid]);
			if (vh_disp_cnt_map.find(vh) == vh_disp_cnt_map.end()) {
				auto pair = std::make_pair(COpenMeshT::Point(v(0), v(1), v(2)), 1);
				vh_disp_cnt_map.insert(std::make_pair(vh, pair));
			}
			else {
				auto pair = vh_disp_cnt_map[vh];
				pair.first += COpenMeshT::Point(v(0), v(1), v(2));
				pair.second++;
			}
		}
	}
	/*------------------- RBF DEFORMATION END ----------------*/

	/*------------------- PRODUCE OUTPUT MESH ----------------*/
	// average deformation
	for (auto itmap = vh_disp_cnt_map.begin(); itmap != vh_disp_cnt_map.end(); ++itmap)
	{
		COpenMeshT::Point p = out_mesh.point(itmap->first) + itmap->second.first / double(itmap->second.second);
		out_mesh.set_point(itmap->first, p);
	}
	/*------------------- PRODUCE OUTPUT MESH END ------------*/
}

void CSSTDeformWrapper::GlobalDeformationLoop(
	COpenMeshT & out_mesh, 
	std::vector<std::vector<COpenMeshT::Point>> & out_cs_pts_list,
	const std::string &fpath,
	double scale)
{
	// Initialize the cross-sections and roi data
	/*--------- CROSS-SECTIONAL CONSTRAINTS ---------------*/
	int numCS;
	std::vector<std::vector<COpenMeshT::Point>> cs_pts_list, emb_cs_pts_list;
	std::vector<int> cs_id_list;
	getCrossSections(emb_cs_pts_list, cs_pts_list, cs_id_list);
	numCS = cs_pts_list.size();
	//out_cs_pts_list = cs_pts_list;
	/*--------- CROSS-SECTIONAL CONSTRAINTS END------------*/
	
	/*------------------- GET PTS IN ROI ---------------------*/
	// define the pairs of cross-sections that bound roi's
	int numROI = numCS - 1;
	DenseMatrixXd cspair_map;
	std::vector<std::vector<int>> rid_keeper_list;
	getROIPts(numROI, cs_id_list, cspair_map, rid_keeper_list);
	/*------------------- GET PTS IN ROI END -----------------*/

	/*------------------- DEFORMATION LOOP AND OUTPUT --------*/
	GlobalDeformationSetup(cs_pts_list);
	std::cout << "setup success" << std::endl;

	CMeshObject *p_mesh_obj = new CMeshObject();
	Eigen::MatrixXd Bmat;
	CGeoCalculator::getBernsteinBasis(4, 0.01, Bmat);

	int numSamples = 10;
	for (int sample_id = 0; sample_id < numSamples; ) {
		sample_id++;
		std::cout << "sample_id = " << sample_id << std::endl;

		std::ostringstream str_sample_id;
		str_sample_id << sample_id;

		// compute displacement
		////CP = [0 0 0; 300 0 0; 200 0 -150; 500 0 -150];
		std::string skel_file = fpath;
		skel_file.append("/skel_ctrp_oo");
		skel_file.append(str_sample_id.str());
		skel_file.append(".txt");
		std::vector<double> data;
		readSkelCtrlPts(skel_file, data);

		DenseMatrixXd def_skel_pts, ctrl_pts;
		ctrl_pts = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>(data.data(), 3, 4);
		ctrl_pts.transposeInPlace();
		ctrl_pts *= scale;
		def_skel_pts = Bmat*ctrl_pts;

		//
		std::vector<DenseMatrixXd> vecs_array;
		std::vector<std::vector<COpenMeshT::Point>> def_cs_pts_list;
		sst_obj_->GetDeformedCrossSectionViaSkelDeformation(
			vecs_array, def_cs_pts_list, emb_cs_pts_list, def_skel_pts, cs_id_list);
		int numHandles = 0;
		for (int i = 0; i < vecs_array.size(); i++)
			numHandles += vecs_array[i].rows();
		DenseMatrixXd vecs_handles(numHandles, 3);
		int cnt = 0;
		for (int i = 0; i < vecs_array.size(); i++) {
			for (int j = 0; j < vecs_array[i].rows(); j++) {
				vecs_handles.row(cnt) = vecs_array[i].row(j);
				cnt++;
			}
		}
		GlobalDeformationSolveAndOutput(sample_id, p_mesh_obj->GetMesh(), vecs_handles, fpath);

		// visualization in the window
		out_cs_pts_list = def_cs_pts_list;
		out_mesh = p_mesh_obj->GetMesh();

		// write file
		//p_mesh_obj->Rescale(1.0 / scale);
		std::map<int, std::vector<double>> Vmap;
		for (auto viter = p_mesh_obj->GetMesh().vertices_begin(); viter != p_mesh_obj->GetMesh().vertices_end(); ++viter) {
			std::vector<double> data;
			auto p = p_mesh_obj->GetMesh().point(*viter);
			data.push_back(p[0]);
			data.push_back(p[1]);
			data.push_back(p[2]);
			int id = p_mesh_obj->GetMesh().data(*viter).get_vlabel();
			Vmap[id] = data;
		}

		//std::string srcFile = "data_deform/test/MeshSF.k";
		//std::string tgtFile = "data_deform/test/MeshSF_skel.K";
		//CFEMeshIO::updateNodesInKeyFile(srcFile, tgtFile, Vmap);

		std::string srcFile = fpath;
		srcFile.append("/MeshSF_rw.k");
		std::string sampleDir = "E:\\Research\\SOSST\\opt_folder_cs1\\sample";
		sampleDir.append(str_sample_id.str());
		// create a new folder
		_mkdir(sampleDir.c_str());
		std::string tgtFile = sampleDir;
		tgtFile.append("/MeshSF_rw_skel.k");
		CFEMeshIO::updateNodesInKeyFile(srcFile, tgtFile, Vmap);
	}
}

void CSSTDeformWrapper::GlobalDeformationSetup(
	const std::vector<std::vector<COpenMeshT::Point>> & cs_pts_list)
{
	std::vector<double> handlepts;
	for (int i = 0; i < cs_pts_list.size(); i++) {
		for (int j = 0; j < cs_pts_list[i].size(); j++)
		{
			handlepts.push_back(cs_pts_list[i][j][0]);
			handlepts.push_back(cs_pts_list[i][j][1]);
			handlepts.push_back(cs_pts_list[i][j][2]);
		}
	}

	DenseMatrixXd pts_handles;
	pts_handles = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>(handlepts.data(), 3, handlepts.size() / 3);
	pts_handles.transposeInPlace();

	CMeshDeformation * deformer = new CMeshDeformation();
	deformer->SetUp(pts_handles);
	p_deformer_list_.push_back(deformer);
}

void CSSTDeformWrapper::GlobalDeformationSolveAndOutput(int sample_id,
	COpenMeshT & out_mesh,
	DenseMatrixXd & vecs_handles,
	const std::string &fpath)
{
	// init
	std::map<COpenMeshT::VertexHandle, std::pair<COpenMeshT::Point, int>> vh_disp_cnt_map;
	out_mesh = triMesh_;

	DenseMatrixXd roi_vecs, roi_pts;
	roi_pts.resize(out_mesh.n_vertices(), 3);
	roi_vecs.resize(out_mesh.n_vertices(), 3);

	int cnt = 0;
	for (auto viter = out_mesh.vertices_begin(); viter != out_mesh.vertices_end(); ++viter)
	{
		roi_pts(cnt, 0) = out_mesh.point(*viter)[0];
		roi_pts(cnt, 1) = out_mesh.point(*viter)[1];
		roi_pts(cnt, 2) = out_mesh.point(*viter)[2];
		cnt++;
	}
	CMeshDeformation * deformer = p_deformer_list_[0];
	deformer->Solve(vecs_handles, roi_pts, roi_vecs);
	roi_pts += roi_vecs;

	/*------------------- PRODUCE OUTPUT MESH ----------------*/
	// average deformation
	cnt = 0;
	for (auto viter = out_mesh.vertices_begin(); viter != out_mesh.vertices_end(); ++viter)
	{
		COpenMeshT::Point p(roi_pts(cnt, 0), roi_pts(cnt, 1), roi_pts(cnt, 2));
		out_mesh.set_point(*viter, p);
		cnt++;
	}
	/*------------------- PRODUCE OUTPUT MESH END ------------*/

}

void CSSTDeformWrapper::GlobalDeformationSolveAndOutput(int sample_id,
	COpenMeshT & out_mesh,
	const DenseMatrixXd & cspair_map,
	const std::vector<DenseMatrixXd> & vecs_array,
	const std::vector<std::vector<int>> & rid_keeper_list,
	const std::string &fpath) 
{
	// init
	std::map<COpenMeshT::VertexHandle, std::pair<COpenMeshT::Point, int>> vh_disp_cnt_map;
	out_mesh = triMesh_;

	int numROI = cspair_map.rows();
	/*------------------- RBF DEFORMATION --------------------*/
	for (int i = 0; i < numROI; i++)
	{
		// ROI
		DenseMatrixXd roi_pts;
		roi_pts.resize(rid_keeper_list[i].size(), 3);
		for (int j = 0; j < rid_keeper_list[i].size(); j++)
		{
			int rid = rid_keeper_list[i][j];
			auto p = embMesh_.point(embMesh_.vertex_handle(pid_vid_map_[rid]));
			roi_pts.row(j) = Vector3d(p[0], p[1], p[2]);
		}

		// RBF DEFORM 
		int id_first = cspair_map(i, 0);
		int id_second = cspair_map(i, 1);
		DenseMatrixXd vecs_h_left, vecs_h_right;
		vecs_h_left = vecs_array[id_first];
		vecs_h_right = vecs_array[id_second];

		// solve
		DenseMatrixXd roi_vecs;
		CMeshDeformation * deformer = p_deformer_list_[i];
		deformer->Solve(vecs_h_left, vecs_h_right, roi_pts, roi_vecs);

		for (int j = 0; j < rid_keeper_list[i].size(); j++)
		{
			int rid = rid_keeper_list[i][j];
			Vector3d v = roi_vecs.row(j);
			auto vh = out_mesh.vertex_handle(pid_vid_map_[rid]);
			if (vh_disp_cnt_map.find(vh) == vh_disp_cnt_map.end()) {
				auto pair = std::make_pair(COpenMeshT::Point(v(0), v(1), v(2)), 1);
				vh_disp_cnt_map.insert(std::make_pair(vh, pair));
			}
			else {
				auto pair = vh_disp_cnt_map[vh];
				pair.first += COpenMeshT::Point(v(0), v(1), v(2));
				pair.second++;
			}
		}
	}
	/*------------------- RBF DEFORMATION END ----------------*/

	/*------------------- PRODUCE OUTPUT MESH ----------------*/
	// average deformation
	for (auto itmap = vh_disp_cnt_map.begin(); itmap != vh_disp_cnt_map.end(); ++itmap)
	{
		COpenMeshT::Point p = out_mesh.point(itmap->first) + itmap->second.first / double(itmap->second.second);
		out_mesh.set_point(itmap->first, p);
	}
	/*------------------- PRODUCE OUTPUT MESH END ------------*/
}

void CSSTDeformWrapper::getROIPts(int numROI, 
	std::vector<int> id_array, 
	DenseMatrixXd & cspair_map, 
	std::vector<std::vector<int>> &rid_keeper_list)
{
	cspair_map.resize(numROI, 2);
	for (int i = 0; i < numROI; i++)
	{
		cspair_map(i, 0) = i;
		cspair_map(i, 1) = i + 1;
	}

	rid_keeper_list.resize(numROI);
	for (auto itmap = pid_skeid_map_.begin(); itmap != pid_skeid_map_.end(); ++itmap)
	{
		for (int i = 0; i < numROI; i++)
		{
			int id_first = cspair_map(i, 0);
			int id_second = cspair_map(i, 1);

			double al_first = 0, al_second = 0;
			for (int j = 0; j < id_array[id_first]; j++)
				al_first += (skelPts_[j + 1] - skelPts_[j]).norm();
			for (int j = 0; j < id_array[id_second]; j++)
				al_second += (skelPts_[j + 1] - skelPts_[j]).norm();

			// 
			if (itmap->second > id_array[id_first] && itmap->second < id_array[id_second]) {
				rid_keeper_list[i].push_back(itmap->first);
			}
			else if (itmap->second == id_array[id_first]) {
				auto p = embMesh_.point(embMesh_.vertex_handle(itmap->first));
				if (p[0] - al_first > 0)
					rid_keeper_list[i].push_back(itmap->first);
			}
			else if (itmap->second == id_array[id_second]) {
				auto p = embMesh_.point(embMesh_.vertex_handle(itmap->first));
				if (p[0] - al_second < 0)
					rid_keeper_list[i].push_back(itmap->first);
			}
		}
	}
}

bool CSSTDeformWrapper::readCSdisplacement(const std::string filename,
	std::vector<COpenMeshT::Point> &con_pts,
	std::vector<COpenMeshT::Point> &con_vecs)
{
	// 
	std::ifstream inFile;
	inFile.open(filename);
	if (!inFile) {
		std::cerr << "Unable to open file " << filename << std::endl;
		return false;
	}
	//else {
	//	std::cerr << "Opened file " << filename << std::endl;
	//}

	int flag = 0;
	while (!inFile.eof()) {
		std::string line;
		getline(inFile, line);
		
		if (inFile.fail()) break;
		std::istringstream iss(line);
		char ch;
		iss >> ch;
		if (ch == 'P') flag = 0;
		else if (ch == 'D') flag = 1;
		else {
			COpenMeshT::Point p(0, 0, 0);
			for (int i = 0; i < 3; i++) iss >> p[i];
			if (flag == 0) con_pts.push_back(p);
			else if (flag == 1) con_vecs.push_back(p);
		}
	}
	return true;
}

bool CSSTDeformWrapper::readCSData(const std::string filename, std::vector<COpenMeshT::Point>& data)
{
	// 
	std::ifstream inFile;
	inFile.open(filename);
	if (!inFile) {
		std::cerr << "Unable to open file " << filename << std::endl;
		return false;
	}
	//else {
	//	std::cerr << "Opened file " << filename << std::endl;
	//}
	while (!inFile.eof()) {
		std::string line;
		getline(inFile, line);
		if (inFile.fail()) break;
		std::istringstream iss(line);
		COpenMeshT::Point p(0, 0, 0);
		for (int i = 0; i < 3; i++) {
			iss >> p[i];
		}
		data.push_back(p);
	}
	return true;
}
bool CSSTDeformWrapper::readCSData(const std::string filename, std::vector<double>& data)
{
	// 
	std::ifstream inFile;
	inFile.open(filename);
	if (!inFile) {
		std::cerr << "Unable to open file " << filename << std::endl;
		return false;
	}
	//else {
	//	std::cerr << "Opened file " << filename << std::endl;
	//}
	while (!inFile.eof()) {
		std::string line;
		getline(inFile, line);
		if (inFile.fail()) break;
		std::istringstream iss(line);
		for (int i = 0; i < 3; i++) { 
			double d = 0;
			iss >> d;
			data.push_back(d);
		}
	}
	return true;
}

bool CSSTDeformWrapper::readSkelCtrlPts(
	const std::string filename,
	std::vector<double> &data) 
{
	// 
	std::ifstream inFile;
	inFile.open(filename);
	if (!inFile) {
		std::cerr << "Unable to open file " << filename << std::endl;
		return false;
	}
	//else {
	//	std::cerr << "Opened file " << filename << std::endl;
	//}
	while (!inFile.eof()) {
		std::string line;
		getline(inFile, line);
		if (inFile.fail()) break;
		std::istringstream iss(line);
		for (int i = 0; i < 3; i++) {
			double d = 0;
			iss >> d;
			data.push_back(d);
		}
	}
	return true;
}

bool CSSTDeformWrapper::writeMeshVertsToTxt(COpenMeshT & out_mesh, int sample_id)
{
	std::string filename = "deform_tests/result/sample";
	std::ostringstream str_sample_id;
	str_sample_id  << sample_id;
	filename.append(str_sample_id.str());
	filename.append(".txt");
	std::ofstream fout(filename);
	for (auto viter = out_mesh.vertices_begin(); viter != out_mesh.vertices_end(); ++viter)
		fout << std::setw(5) << viter->idx() << " " << out_mesh.point(*viter) << "\n";
	fout.close();
	return true;
}


//void CSSTDeformWrapper::GetEmbeddedMesh(COpenMeshT & triMesh)
//{
//	triMesh = embMesh_;
//}

void CSSTDeformWrapper::getCrossSections(
	std::vector<std::vector<COpenMeshT::Point>> & emb_cs_pts_list,
	std::vector<std::vector<COpenMeshT::Point>> & cs_pts_list, 
	std::vector<int> & id_cs)
{
	/*--------- CROSS-SECTIONAL CONSTRAINTS ---------------*/
	
	//std::vector<int> id_cs;
	//std::cout << skelPts_.size() << std::endl;
	for (int i = 0; i < skelPts_.size(); i+=10)
	{
		id_cs.push_back(i);
		// points
		std::vector<COpenMeshT::Point> cs_pts;
		std::string file_pts_handle = "data_deform/test";
		file_pts_handle.append("/pts_handle.txt");
		readCSData(file_pts_handle, cs_pts);
		double arclength = 0;
		for (int j = 0; j < i; j++)
			arclength += (skelPts_[j + 1] - skelPts_[j]).norm();
		for (int j = 0; j < cs_pts.size(); j++) {
			cs_pts[j][0] += arclength;
		}
		emb_cs_pts_list.push_back(cs_pts);
	}
	/*--------- CROSS-SECTIONAL CONSTRAINTS END------------*/

	/*------------- OUTPUT CROSS-SECTIONS -----------------*/
	for (int i = 0; i < emb_cs_pts_list.size(); i++) {
		std::vector<COpenMeshT::Point> cs_pts;
		for (int j = 0; j < emb_cs_pts_list[i].size(); j++)
		{
			Vector3d pt;
			Vector3d emb_pt(emb_cs_pts_list[i][j][0], emb_cs_pts_list[i][j][1], emb_cs_pts_list[i][j][2]);
			sst_obj_->DecodePosition(pt, emb_pt, id_cs[i]);
			COpenMeshT::Point p(pt(0), pt(1), pt(2));
			cs_pts.push_back(p);
		}
		cs_pts_list.push_back(cs_pts);
	}
	/*------------- OUTPUT CROSS-SECTIONS END -------------*/
}

void CSSTDeformWrapper::outputEmbeddingMesh(std::string fpath, double scale)
{
	//std::string srcFile = "data_deform/test/MeshSF.k";
	std::string srcFile = fpath;
	srcFile.append("/MeshSF.k");

	// get embedded mesh
	CMeshObject *p_mesh_obj = new CMeshObject();
	p_mesh_obj->GetMesh() = embMesh_;
	//p_mesh_obj->Rescale(1.0 / scale);
	//p_mesh_obj->swapXYCoords();
	
	// construct vmap
	std::map<int, std::vector<double>> Vmap;
	for (auto viter = p_mesh_obj->GetMesh().vertices_begin(); viter != p_mesh_obj->GetMesh().vertices_end(); ++viter) {
		std::vector<double> data;
		auto p = p_mesh_obj->GetMesh().point(*viter);
		data.push_back(p[0]);
		data.push_back(p[1]);
		data.push_back(p[2]);
		int id = p_mesh_obj->GetMesh().data(*viter).get_vlabel();
		Vmap[id] = data;
	}

	// create a new folder
	std::string tgtFile = fpath;
	tgtFile.append("/MeshSF_embed.k");
	CFEMeshIO::updateNodesInKeyFile(srcFile, tgtFile, Vmap);
}