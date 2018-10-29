#include "../DataColle/data_io.h"
#include "File_oi.h"
#include "pugixml\pugixml.hpp"
#include "pugixml\pugiconfig.hpp"
#include "ini_parser\ini.hpp"

SwitchType CFileIO::read_config_file(std::string fname)
{
	std::ifstream in_file(fname);
	if (!in_file.is_open())
	{
		std::cerr << "Failed to open the config_switch file" << std::endl;
		return SwitchType::sstExit;
	}
	std::string line;
	std::getline(in_file, line);
	in_file.close();
	if (strcmp(line.c_str(), "GLOBAL") == 0) {
		return SwitchType::sstGlobal; // global
	}
	else if (strcmp(line.c_str(), "LOCAL") == 0) {
		return SwitchType::sstLocal; // local 
	}
	else if (strcmp(line.c_str(), "PENDING") == 0) {
		return  SwitchType::sstPending; // pending
	}
	else if (strcmp(line.c_str(), "UPDATECS") == 0) {
		return  SwitchType::sstUpdateCS; // pending
	}
	else if (strcmp(line.c_str(), "INIT") == 0) {
		return  SwitchType::sstInit; // pending
	}
	else {
		return SwitchType::sstExit; // exit
	}
}

bool CFileIO::xml_read_skeleton(std::string fname, CSkeleton & s)
{
	pugi::xml_document doc;
	pugi::xml_parse_result result = doc.load_file(fname.c_str());
	if (!result)
	{
		std::cerr << "Error at CFileIO::xml_read_skeleton()" << std::endl;
		return false;
	}

	std::vector<COpenMeshT::Point> ctrlPts, skelPts;
	std::vector<double>ts;
	pugi::xml_node skel_node = doc.child("DUT_AutoMorpher").child("Skeleton");
	pugi::xml_node pt_node = skel_node.first_child();
	for (; pt_node; pt_node = pt_node.next_sibling())
	{
		COpenMeshT::Point pt(0.0, 0.0, 0.0);
		double t;
		std::istringstream sx(pt_node.attribute("X").value());
		sx >> pt[0];
		std::istringstream sy(pt_node.attribute("Y").value());
		sy >> pt[1];
		std::istringstream sz(pt_node.attribute("Z").value());
		sz >> pt[2];
		std::istringstream st(pt_node.attribute("t").value());
		st >> t;
		if (strcmp(pt_node.name(), "Point") == 0) {
			skelPts.push_back(pt);
			ts.push_back(t);
		}
		else if (strcmp(pt_node.name(), "ctrlPt")==0) {
			ctrlPts.push_back(pt);
		}
	}

	// change
	bool is_bezier = false;
	if (strcmp(skel_node.attribute("IsBezier").value(), "T")==0) {
		is_bezier = true;
		std::cout << "bezier true" << std::endl;
	}
	if (is_bezier/* && ts.size() == 0*/) {
		std::cout << ctrlPts.size() << std::endl;
		int numSamples = 100;
		CSkeleton tmp_s(ctrlPts, numSamples);
		s.CopyFrom(tmp_s);
	}
	else {
		s.InitData(skelPts);
	}
	// change/
	return true;
}

bool CFileIO::xml_write_skeleton_polyline(std::string fname, const CSkeleton & s)
{
	std::vector<COpenMeshT::Point> skelPts = s.GetSkeletalPts();

	pugi::xml_document dst_doc;
	pugi::xml_node decl = dst_doc.prepend_child(pugi::node_declaration);
	decl.append_attribute("version") = "1.0";
	decl.append_attribute("encoding") = "UTF-8";

	dst_doc.append_child("DUT_AutoMorpher");
	pugi::xml_node dst_ss = dst_doc.child("DUT_AutoMorpher");
	pugi::xml_node skel_node = dst_ss.append_child("Skeleton");
	skel_node.append_attribute("IsBezier").set_value("F");
	for (int i = 0; i < skelPts.size(); i++) {
		pugi::xml_node vertex = dst_ss.child("Skeleton").append_child("Point");
		vertex.append_attribute("X");
		vertex.append_attribute("Y");
		vertex.append_attribute("Z");
		int cnt = 0;
		for (pugi::xml_attribute attr = vertex.first_attribute();
			attr; attr = attr.next_attribute(), cnt++) {
			attr.set_value(float(skelPts[i][cnt]));
		}
		vertex.prepend_attribute("ID");
		pugi::xml_attribute attr = vertex.first_attribute();
		attr.set_value(i);
	}
	std::cout << "Saving result: " << dst_doc.save_file(fname.c_str()) << std::endl;

	return true;
}

bool CFileIO::xml_write_skeleton_bezier(std::string fname, const CSkeleton & s)
{
	std::vector<COpenMeshT::Point> ctrlPts = s.GetCtrlPts();
	std::vector<COpenMeshT::Point> skelPts = s.GetSkeletalPts();
	std::vector<double> ts = s.GetTParas();
	assert(skelPts.size() == ts.size());

	pugi::xml_document dst_doc;
	pugi::xml_node decl = dst_doc.prepend_child(pugi::node_declaration);
	decl.append_attribute("version") = "1.0";
	decl.append_attribute("encoding") = "UTF-8";

	dst_doc.append_child("DUT_AutoMorpher");
	pugi::xml_node dst_ss = dst_doc.child("DUT_AutoMorpher");
	pugi::xml_node skel_node = dst_ss.append_child("Skeleton");
	skel_node.append_attribute("IsBezier").set_value("T");
	for (int i = 0; i < ctrlPts.size(); i++) {
		pugi::xml_node vertex = dst_ss.child("Skeleton").append_child("ctrlPt");
		vertex.append_attribute("X");
		vertex.append_attribute("Y");
		vertex.append_attribute("Z");
		int cnt = 0;
		for (pugi::xml_attribute attr = vertex.first_attribute();
			attr; attr = attr.next_attribute(), cnt++) {
			attr.set_value(float(ctrlPts[i][cnt]));
		}
		vertex.prepend_attribute("ID");
		pugi::xml_attribute attr = vertex.first_attribute();
		attr.set_value(i);
	}
	for (int i = 0; i < skelPts.size(); i++) {
		pugi::xml_node vertex = dst_ss.child("Skeleton").append_child("Point");
		vertex.append_attribute("X");
		vertex.append_attribute("Y");
		vertex.append_attribute("Z");
		int cnt = 0;
		for (pugi::xml_attribute attr = vertex.first_attribute();
			attr; attr = attr.next_attribute(), cnt++) {
			attr.set_value(float(skelPts[i][cnt]));
		}
		vertex.append_attribute("t").set_value(float(ts[i]));
		vertex.prepend_attribute("ID");
		pugi::xml_attribute attr = vertex.first_attribute();
		attr.set_value(i);
	}
	std::cout << "Saving result: " << dst_doc.save_file(fname.c_str()) << std::endl;

	return true;
}

bool CFileIO::xml_read_a_section(std::string fname, const int sid, CCrossSection & cs)
{
	pugi::xml_document doc;
	pugi::xml_parse_result result = doc.load_file(fname.c_str());
	if (!result) {
		std::cerr << "Error at CFileIO::xml_read_skeleton()" << std::endl;
		return false;
	}

	std::vector<COpenMeshT::Point> cs_pts, emb_cs_pts;
	pugi::xml_node cs_node = doc.child("DUT_AutoMorpher").child("Section");	
	for (cs_node; cs_node; cs_node = cs_node.next_sibling())
	{
		std::string spid_value = cs_node.attribute("SPID").value();
		int tmp_sid = atoi(spid_value.c_str());
		std::string id_value = cs_node.attribute("ID").value();
		int tmp_csid = atoi(id_value.c_str());
		if (sid == tmp_sid) {
			//std::cout << "ID" << tmp_csid << ", SPID " << tmp_sid << std::endl;
			pugi::xml_node pt_node = cs_node.first_child();
			for (; pt_node; pt_node = pt_node.next_sibling())
			{
				COpenMeshT::Point pt(0.0, 0.0, 0.0);
				std::istringstream sx(pt_node.attribute("X").value());
				sx >> pt[0];
				std::istringstream sy(pt_node.attribute("Y").value());
				sy >> pt[1];
				std::istringstream sz(pt_node.attribute("Z").value());
				sz >> pt[2];
				if (strcmp(pt_node.name(), "Point") == 0)
					cs_pts.push_back(pt);
				else if (strcmp(pt_node.name(), "EmbPoint") == 0)
					emb_cs_pts.push_back(pt);
			}
			cs.SetSid(tmp_sid);
			cs.SetEmbProfPts(emb_cs_pts);
			cs.SetProfPts(cs_pts);
			if (strcmp(cs_node.attribute("IsClosed").value(), "T") == 0)
				cs.SetClosed();
			break;
		}
	}
	return true;
}

bool CFileIO::xml_read_sections(std::string fname, std::vector<CCrossSection>& cs_list)
{
	cs_list.clear();

	pugi::xml_document doc;
	pugi::xml_parse_result result = doc.load_file(fname.c_str());
	if (!result) {
		std::cerr << "Error at CFileIO::xml_read_sections()" << std::endl;
		return false;
	}

	pugi::xml_node cs_node = doc.child("DUT_AutoMorpher").child("Section");
	for (cs_node; cs_node; cs_node = cs_node.next_sibling())
	{
		//
		std::vector<COpenMeshT::Point> cs_pts, emb_cs_pts;
		int csid, spid;
		bool is_closed, is_deformed;
		//
		std::string str_id_value = cs_node.attribute("ID").value();
		csid = atoi(str_id_value.c_str());
		std::string str_rid_value = cs_node.attribute("SPID").value();
		spid = atoi(str_rid_value.c_str());
		if (strcmp(cs_node.attribute("IsClosed").value(),"T") == 0)
			is_closed = true;
		else 
			is_closed = false;
		if (strcmp(cs_node.attribute("IsDeformed").value(), "T") == 0)
			is_deformed = true;
		else
			is_deformed = false;

		std::cout << "deformed CS ID: " << csid << ", SPID: " << spid << ", IsClosed: " << is_closed << std::endl;

		//std::cout << node.attribute("ID").value() << "   Y" << std::endl;
		pugi::xml_node pt_node = cs_node.first_child();
		for (; pt_node; pt_node = pt_node.next_sibling())
		{
			COpenMeshT::Point pt(0.0, 0.0, 0.0);
			std::istringstream sx(pt_node.attribute("X").value());
			sx >> pt[0];
			std::istringstream sy(pt_node.attribute("Y").value());
			sy >> pt[1];
			std::istringstream sz(pt_node.attribute("Z").value());
			sz >> pt[2];
			if (strcmp(pt_node.name(), "Point") == 0)
				cs_pts.push_back(pt);
			else if (strcmp(pt_node.name(), "EmbPoint") == 0)
				emb_cs_pts.push_back(pt);
		}

		//
		CCrossSection cs;
		cs.SetSid(spid);
		cs.SetEmbProfPts(emb_cs_pts);
		cs.SetProfPts(cs_pts);
		if (is_closed) cs.SetClosed();
		// 
		cs_list.push_back(cs);
	}
	return true;
}

bool CFileIO::xml_write_sections(std::string fname, const std::vector<CCrossSection> & cs_list)
{
	pugi::xml_document dst_doc;
	pugi::xml_node decl = dst_doc.prepend_child(pugi::node_declaration);
	// heading
	decl.append_attribute("version") = "1.0";
	decl.append_attribute("encoding") = "UTF-8";
	// root node
	dst_doc.append_child("DUT_AutoMorpher");
	pugi::xml_node dst_ss = dst_doc.child("DUT_AutoMorpher");


	// write cross_section nodes
	int csid = 0;
	std::vector<CCrossSection>::const_iterator cs_iter = cs_list.begin();
	for (; cs_iter != cs_list.end(); ++cs_iter, ++csid) {
		std::vector<COpenMeshT::Point> cs_pts, emb_cs_pts;
		cs_pts = cs_iter->GetProfPts();
		emb_cs_pts = cs_iter->GetEmbProfPts();
		int sid = cs_iter->GetSid();

		// cross section node and their attributes
		pugi::xml_node tmp_cs = dst_ss.append_child("Section");
		tmp_cs.append_attribute("ID").set_value(csid);
		tmp_cs.append_attribute("SPID").set_value(sid);
		if (cs_iter->IsClosed())
			tmp_cs.append_attribute("IsClosed").set_value("T");
		else
			tmp_cs.append_attribute("IsClosed").set_value("F");
		if (cs_iter->IsDeformed())
			tmp_cs.append_attribute("IsDeformed").set_value("T");
		else
			tmp_cs.append_attribute("IsDeformed").set_value("F");

		// add point nodes and their attributes
		for (int i = 0; i < cs_pts.size(); i++) {
			pugi::xml_node vertex = tmp_cs.append_child("Point");
			vertex.append_attribute("X").set_value(float(cs_pts[i][0]));
			vertex.append_attribute("Y").set_value(float(cs_pts[i][1]));
			vertex.append_attribute("Z").set_value(float(cs_pts[i][2]));
			vertex.prepend_attribute("ID").set_value(i);
		}
		for (int i = 0; i < emb_cs_pts.size(); i++) {
			pugi::xml_node vertex = tmp_cs.append_child("EmbPoint");
			vertex.append_attribute("X").set_value(float(emb_cs_pts[i][0]));
			vertex.append_attribute("Y").set_value(float(emb_cs_pts[i][1]));
			vertex.append_attribute("Z").set_value(float(emb_cs_pts[i][2]));
			vertex.prepend_attribute("ID").set_value(i);
		}
	}
	std::cout << "Saving Cross-sections: " << dst_doc.save_file(fname.c_str()) << std::endl;

	return true;
}

bool CFileIO::write_mesh(std::string fname, COpenMeshT & m, bool is_emb)
{
	// open file
	std::ofstream  fout_file;
	fout_file.open(fname);
	if (!fout_file.is_open()) {
		std::cerr << "PrintMesh(): cannot open the file" << std::endl;
		return false;
	}

	// write
	//fout_file << "##Mesh" << std::endl;
	//fout_file << "#Vertices" << std::endl;
	if (!is_emb) {
		for (auto viter = m.vertices_begin();
			viter != m.vertices_end(); ++viter)
		{
			fout_file << viter->idx() << " "
				<< m.point(*viter)
				<< std::endl;
		}
	}
	else {
		for (auto viter = m.vertices_begin();
			viter != m.vertices_end(); ++viter)
		{
			fout_file << viter->idx() << " "
				<< m.data(*viter).get_emb_coord()
				<< std::endl;
		}
	}
	//fout_file << "#Faces" << std::endl;
	//for (auto fiter = m.faces_begin();
	//	fiter != m.faces_end(); ++fiter)
	//{
	//	for (auto ccwfv_iter = m.cfv_ccwbegin(*fiter);
	//		ccwfv_iter != m.cfv_ccwend(*fiter); ++ccwfv_iter)
	//	{
	//		fout_file << ccwfv_iter->idx() << " ";
	//	}
	//	fout_file << std::endl;
	//}
	fout_file.close();
	return true;
}

bool CFileIO::write_mesh_to_stl(std::string fname, COpenMeshT & mesh, bool is_emb)
{
	std::ofstream out_file(fname);
	if (out_file.is_open())
	{
		mesh.update_face_normals();
		char* obj_name = "SF";
		out_file << "solid " << obj_name << std::endl;
		for (COpenMeshT::FaceIter fiter = mesh.faces_begin();
			fiter != mesh.faces_end(); ++fiter)
		{
			out_file << "    facet normal " << mesh.normal(*fiter) << std::endl;
			out_file << "        outer loop " << std::endl;
			for (COpenMeshT::FaceVertexCCWIter fviter = mesh.fv_ccwbegin(*fiter);
				fviter != mesh.fv_ccwend(*fiter);  ++fviter)
			{
				if (!is_emb) 
					out_file << "            vertex " << mesh.point(*fviter) << std::endl;
				else
					out_file << "            vertex " << mesh.data(*fviter).get_emb_coord() << std::endl;
			}
			out_file << "        endloop" << std::endl;
			out_file << "    endfacet" << std::endl;
		}
		out_file << "endsolid" << std::endl;
		out_file.close();
		return true;
	}
	return false;
}

bool CFileIO::write_mesh_to_obj(std::string fname, COpenMeshT & m, bool is_emb)
{
	// open file
	std::ofstream  fout_file;
	fout_file.open(fname);
	if (!fout_file.is_open()) {
		std::cerr << "PrintMesh(): cannot open the file" << std::endl;
		return false;
	}

	// write
	//fout_file << "##Mesh" << std::endl;
	//fout_file << "#Vertices" << std::endl;
	if (!is_emb) {
		for (auto viter = m.vertices_begin();
			viter != m.vertices_end(); ++viter)
		{
			fout_file<< "v "
				<< m.point(*viter)
				<< std::endl;
		}
	}
	else {
		for (auto viter = m.vertices_begin();
			viter != m.vertices_end(); ++viter)
		{
			fout_file << viter->idx() << " "
				<< m.data(*viter).get_emb_coord()
				<< std::endl;
		}
	}
	for (auto fiter = m.faces_begin();
		fiter != m.faces_end(); ++fiter)
	{
		fout_file << "f ";
		for (COpenMeshT::FaceVertexCCWIter fviter = m.fv_ccwbegin(*fiter);
			fviter != m.fv_ccwend(*fiter);  ++fviter)
			fout_file << fviter->idx()+1 << " ";
		fout_file << std::endl;
	}
	fout_file.close();
	return true;
}

//bool CFileIO::read_config_file_ini(std::string fname)
//{
//	// open file
//	std::ifstream  in_file;
//	in_file.open(fname);
//	if (!in_file.is_open()) {
//		std::cerr << "read_config_file_ini(): cannot open the file" << std::endl;
//		return false;
//	}
//
//	std::string line;
//	std::getline(in_file, line);
//
//	size_t pos = line.find_first_of("=");
//	char line_name[125];
//	char line_path[125];
//	for (int i = 0; i < pos; i++) {
//		line_name[i] = line[i];
//	}
//	for (int i = pos+1; i < line.size(); i++) {
//		line_path[i] = line[i];
//	}
//	std::cout << line_name << "\t" << line_path << std::endl;
//	
//
//	in_file.close();
//	return false;
//}

//SwitchType CFileIO::xml_read_config_file(const std::string fname, CFileIO::InputPaths &in_paths)
//{
//	pugi::xml_document doc;
//	pugi::xml_parse_result result = doc.load_file(fname.c_str());
//	if (!result)
//	{
//		std::cout << fname << std::endl;
//		std::cerr << "Error at CFileIO::xml_read_config_file()" << std::endl;
//		return SwitchType::sstExit;
//	}
//
//	// data folder
//	pugi::xml_node folder_node = doc.child("DUT_AutoMorpher").child("DataFolder");
//	std::string data_path = folder_node.attribute("Path").value();
//
//	// mesh 
//	pugi::xml_node mesh_node = doc.child("DUT_AutoMorpher").child("Mesh");
//	std::string mesh_path = mesh_node.attribute("Path").value();
//
//	// solver
//	pugi::xml_node solver_node = doc.child("DUT_AutoMorpher").child("Solver");
//	std::string solver_type = solver_node.attribute("Type").value();
//
//	// output 
//	pugi::xml_node output_node = doc.child("DUT_AutoMorpher").child("Mesh");
//	std::string output_path = output_node.attribute("Path").value();
//
//	// set paths and output
//	if (strcmp(solver_type.c_str(), "Init") == 0)
//	{
//		std::string src_skel_path = solver_node.child("Src").attribute("Path").value();
//
//		std::cout << "SolverType: " << "Initilization" << std::endl;
//		std::cout << "MeshInPath: " << mesh_path << std::endl;
//		std::cout << "SrcPath: " << src_skel_path << std::endl;
//		
//		// set path
//		in_paths.path_to_config_in = fname;
//		in_paths.path_to_mesh_in = mesh_path;
//		in_paths.path_to_src = src_skel_path;
//
//		return SwitchType::sstInit;
//	}
//	else if (strcmp(solver_type.c_str(), "Pending") == 0)
//	{
//		return SwitchType::sstPending;
//	}
//	else if (strcmp(solver_type.c_str(), "UpdateCS") == 0)
//	{
//		std::cerr << "sstUpdate is not ready to use" << std::endl;
//		return SwitchType::sstUpdateCS;
//	}
//	else if (strcmp(solver_type.c_str(), "Global") == 0)
//	{
//		std::string src_skel_path, dst_skel_path;
//		src_skel_path = solver_node.child("Src").attribute("Path").value();
//		dst_skel_path = solver_node.child("Dst").attribute("Path").value();
//		std::cout << "SolverType: " << "Global" << std::endl;
//		std::cout << "MeshInPath: " << mesh_path << std::endl;
//		std::cout << "SrcPath: " << src_skel_path << std::endl;
//		std::cout << "DstPath: " << dst_skel_path << std::endl;
//
//		// set path
//		in_paths.path_to_config_in = fname;
//		in_paths.path_to_mesh_in = mesh_path;
//		in_paths.path_to_src = src_skel_path;
//		in_paths.path_to_dst = dst_skel_path;
//
//		return SwitchType::sstGlobal;
//	}
//	else if (strcmp(solver_type.c_str(), "Local") == 0)
//	{
//		std::string src_cs_path, dst_cs_path;
//		src_cs_path = solver_node.child("Src").attribute("Path").value();
//		dst_cs_path = solver_node.child("Dst").attribute("Path").value();
//		std::cout << "SolverType: " << "Local" << std::endl;
//		std::cout << "MeshInPath: " << mesh_path << std::endl;
//		std::cout << "SrcPath: " << src_cs_path << std::endl;
//		std::cout << "DstPath: " << dst_cs_path << std::endl;
//
//		// set path
//		in_paths.path_to_config_in = fname;
//		in_paths.path_to_mesh_in = mesh_path;
//		in_paths.path_to_src = src_cs_path;
//		in_paths.path_to_dst = dst_cs_path;
//
//		return SwitchType::sstLocal;
//	}
//	else
//	{
//		return SwitchType::sstExit;
//	}
//}

SwitchType CFileIO::ini_read_config_file(const std::string fname, CFileIO::InputPaths & in_paths)
{
	INI::Parser ini_parser(fname.c_str());
	INI::Level topLevel = ini_parser.top();

	if (topLevel("Exit").depth)
	{
		return SwitchType::sstExit;
	}
	else if (topLevel("Init").depth)
	{
		std::cout << "State: " << topLevel("Init")["State"] << std::endl;

		in_paths.path_to_config_in = fname;
		in_paths.path_to_mesh_in = topLevel("Init")["Mesh"];
		in_paths.path_to_outfolder = topLevel("Init")["OutputFolder"];
		in_paths.path_to_src = topLevel("Init")["SrcPath"];
		
		return SwitchType::sstInit;
	}
	else if (topLevel("UpdateCS").depth)
	{
		std::cout << "State: " << topLevel("UpdateCS")["State"] << std::endl;

		in_paths.path_to_config_in = fname;
		in_paths.path_to_mesh_in = topLevel("UpdateCS")["Mesh"];
		in_paths.path_to_outfolder = topLevel("UpdateCS")["OutputFolder"];
		in_paths.path_to_src = topLevel("UpdateCS")["SrcPath"];

		return SwitchType::sstUpdateCS;
	}
	else if (topLevel("Global").depth)
	{
		std::cout << "State: " << topLevel("Global")["State"] << std::endl;

		in_paths.path_to_config_in = fname;
		in_paths.path_to_mesh_in = topLevel("Global")["Mesh"];
		//in_paths.path_to_outfolder = topLevel("Global")["OutputFolder"];
		in_paths.path_to_src = topLevel("Global")["SrcPath"];
		in_paths.path_to_dst = topLevel("Global")["DstPath"];

		return SwitchType::sstGlobal;
	}
	else if (topLevel("Local").depth)
	{
		std::cout << "State: " << topLevel("Local")["State"] << std::endl;

		in_paths.path_to_config_in = fname;
		in_paths.path_to_mesh_in = topLevel("Local")["Mesh"];
		//in_paths.path_to_outfolder = topLevel("Local")["OutputFolder"];
		in_paths.path_to_src = topLevel("Local")["SrcPath"];
		in_paths.path_to_dst = topLevel("Local")["DstPath"];

		return SwitchType::sstLocal;
	}
	else{
		std::cout << "pending" << std::endl;
		return SwitchType::sstPending;
	}
}

bool CFileIO::ini_write_config_file_to_pending(const std::string fname)
{
	std::ofstream out_file;
	out_file.open(fname);
	if (out_file.is_open()) {
		out_file << "[Pending]\n";
		out_file << "State = Pending\n";
		out_file.close();
		return true;
	}

	return false;
}

bool CFileIO::output_SST(CSstObject* sst, CFileIO::OutputPaths out_paths)
{
	// get path
	std::string fname_src_skeleton_ = out_paths.path_to_src_skel;
	std::string fname_dst_skeleton_ = out_paths.path_to_def_skel;
	std::string fname_src_cross_sections_ = out_paths.path_to_src_cs;
	std::string fname_dst_cross_sections_ = out_paths.path_to_def_skel;
	std::string fname_src_mesh_ = out_paths.path_to_mesh_in;
	std::string fname_dst_mesh_ = out_paths.path_to_mesh_out;
	std::string fname_src_emb_mesh_ = out_paths.path_to_mesh_emb;

	CSkeleton s = sst->GetSkeleton();
	std::vector<CCrossSection> cs_list;
	sst->GetCrossSections(cs_list);
	COpenMeshT m = sst->GetMesh();
	

	// print src
	if (!CFileIO::xml_write_skeleton_polyline(fname_src_skeleton_, s))
		return false;
	if (!CFileIO::xml_write_sections(fname_src_cross_sections_, cs_list))
		return false;
	if (!write_mesh_to_stl(fname_src_mesh_, m, false))
		return false;
	if (!write_mesh_to_stl(fname_src_emb_mesh_, m, true) && sst->IsEncode())
		return false;

	if (sst->IsDeformed())
	{
		// def data
		CSkeleton ds = sst->GetDefSkeleton();
		COpenMeshT dm; sst->OutputDefMesh(dm);
		std::vector<CCrossSection> def_cross_sections;
		for (int i = 0; i < cs_list.size(); i++)
		{
			if (cs_list[i].IsDeformed()) {
				def_cross_sections.push_back(cs_list[i]);
			}
		}
		// print def data
		if (!CFileIO::xml_write_skeleton_polyline(fname_dst_skeleton_, ds))
			return false;
		if (!CFileIO::xml_write_sections(fname_dst_cross_sections_, def_cross_sections))
			return false;
		if (!write_mesh_to_stl(fname_dst_mesh_, dm, false))
			return false;
		//if (!write_mesh(fname_dst_emb_mesh_, dm, true) && sst->IsEncode())
		//	return false;
	}

	return true;
}
bool CFileIO::output_CS(CSstObject * sst, std::string fname)
{
	std::vector<CCrossSection> cs_list;
	sst->GetCrossSections(cs_list);
	if (!CFileIO::xml_write_sections(fname, cs_list))
		return false;

	return true;
}
void CFileIO::set_output_path(std::string path, CFileIO::OutputPaths &output_paths)
{
	// set paths

	//std::string folder_path;
	std::string folder_path;
	size_t a = 0;
	if (folder_path.length() != 0) {
		a = path.find_last_of("\\");
		std::cout << "a = path.find_last_of() --> " << a << std::endl;
		for (int i = 0; i < a; i++) {
			//std::cout << path[i];
			folder_path.push_back(path[i]);
		}
	}
	folder_path.append("\\OutputData");

	output_paths.path_to_folder = folder_path;
	
	output_paths.path_to_mesh_in = folder_path;
	output_paths.path_to_mesh_in.append("\\Src_Mesh.stl");
	
	output_paths.path_to_mesh_emb = folder_path;
	output_paths.path_to_mesh_emb.append("\\Emb_Mesh.stl");
	
	output_paths.path_to_mesh_out = folder_path;
	output_paths.path_to_mesh_out.append("\\Gen_Mesh.stl");
	
	output_paths.path_to_src_skel = folder_path;
	output_paths.path_to_src_skel.append("\\src_skeleton.xml");
	
	output_paths.path_to_def_skel = folder_path;
	output_paths.path_to_def_skel.append("\\dst_skeleton.xml");
	
	output_paths.path_to_src_cs = folder_path;
	output_paths.path_to_src_cs.append("\\src_cross_sections.xml");
	
	output_paths.path_to_def_cs = folder_path;
	output_paths.path_to_def_cs.append("\\dst_cross_sections.xml");
}
//
