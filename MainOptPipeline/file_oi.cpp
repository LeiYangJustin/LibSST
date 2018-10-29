#include"File_oi.h"

void file_io::Read_section_data_from_xml(std::string Dessection_fname)
{

	pugi::xml_document doc;
    doc.load_file(Dessection_fname.c_str());
	pugi::xml_node random_num = doc.child("DUT_AutoMorpher").child("Time");

	for (random_num; random_num; random_num = random_num.next_sibling())
	{
		times_ += 1;
		pugi::xml_node cs_node = random_num.child("Section");
		for (cs_node; cs_node; cs_node = cs_node.next_sibling())
		{
			Basic_data_of_DesSection data;
			std::string id;//id
			id = cs_node.attribute("ID").value();
			int i = atoi(id.c_str());
			data.id_ = i;

			std::string spid;//spid
			spid = cs_node.attribute("SPID").value();
			int	s = atoi(spid.c_str());
			data.spid_ = s;

			//if (strcmp(cs_node.attribute("IsClosed").value(), "T") == 0)//closed
			//	data.isclosed_ = true;
			//else
			//	data.isclosed_ = false;

			//if (strcmp(cs_node.attribute("IsDeformed").value(), "T") == 0)//deformed
			//	data.isdeformed_ = true;
			//else
			//	data.isdeformed_ = false;

			std::string pts_num;//pts_number_on_single_edge
			pts_num = cs_node.attribute("EndPoint").value();
			data.pts_num_on_single_edge = atoi(pts_num.c_str());

			std::string X;//x_coordinate
			X = cs_node.attribute("X").value();
			data.x_coordinate_ = atof(X.c_str());
			
			std::string length;
			length = cs_node.attribute("Length").value();
			data.length_ = atof(length.c_str());

			std::string width;
			width = cs_node.attribute("Width").value();
			data.width_ = atof(width.c_str());

			pugi::xml_node ctrlpt_node = cs_node.first_child();
			for (;ctrlpt_node; ctrlpt_node=ctrlpt_node.next_sibling())
			{
				std::string pt_spid;//ctrl_id
				pt_spid = ctrlpt_node.attribute("ID").value();
				int i = atoi(pt_spid.c_str());
				data.ctrl_pid_list_.push_back(i);
				std::string ctrlpt_def;
				ctrlpt_def = ctrlpt_node.attribute("Deformation").value();
				double def = atof(ctrlpt_def.c_str());
				data.ctrl_pt_def_list_.push_back(def);
			}
			def_section_list_.push_back(data);
		}
	}
}
std::vector<Basic_data_of_DesSection> file_io::get_def_section_data_list()
{
	return def_section_list_;
}
void file_io::set_basic_data(std::string Dessection_fname)
{
	Read_section_data_from_xml(Dessection_fname);
}
void file_io::Write_xml(Basic_data_of_DesSection def_section, std::vector<Point2> def_pt_list, std::string outfname)
{
	if (xml_is_empty_)
	{
		pugi::xml_document def_doc;
		pugi::xml_node decl = def_doc.prepend_child(pugi::node_declaration);
		decl.append_attribute("version") = "1.0";
		decl.append_attribute("encoding") = "UTF-8";
		def_doc.append_child("DUT_AutoMorpher");
		pugi::xml_node def_ss = def_doc.child("DUT_AutoMorpher");

		pugi::xml_node tmp_cs = def_ss.append_child("Section");
		tmp_cs.append_attribute("ID").set_value(def_section.id_);
		tmp_cs.append_attribute("SPID").set_value(def_section.spid_);

		if (def_section.isclosed_)
			tmp_cs.append_attribute("IsClosed").set_value("T");
		else
			tmp_cs.append_attribute("IsClosed").set_value("F");
		if (def_section.isdeformed_)
			tmp_cs.append_attribute("IsDeformed").set_value("T");
		else
			tmp_cs.append_attribute("IsDeformed").set_value("F");
		for (int i = 0; i < def_pt_list.size(); i++)
		{
			pugi::xml_node ctrl_pts = tmp_cs.append_child("EmbPoint");
			ctrl_pts.append_attribute("ID").set_value(i);
			ctrl_pts.append_attribute("X").set_value(def_section.x_coordinate_);
			ctrl_pts.append_attribute("Y").set_value(def_pt_list[i].y_);
			ctrl_pts.append_attribute("Z").set_value(def_pt_list[i].z_);
		}
		std::cout << "Saving Cross-sections: " << def_doc.save_file(outfname.c_str()) << std::endl;
		xml_is_empty_ = false;
	}
	else
	{
		pugi::xml_document def_doc;
		def_doc.load_file(outfname.c_str());
		pugi::xml_node def_ss = def_doc.child("DUT_AutoMorpher");

		pugi::xml_node decl = def_doc.prepend_child(pugi::node_declaration);
		decl.append_attribute("version") = "1.0";
		decl.append_attribute("encoding") = "UTF-8";
		

		pugi::xml_node tmp_cs = def_ss.append_child("Section");
		tmp_cs.append_attribute("ID").set_value(def_section.id_);
		tmp_cs.append_attribute("SPID").set_value(def_section.spid_);

		if (def_section.isclosed_)
			tmp_cs.append_attribute("IsClosed").set_value("T");
		else
			tmp_cs.append_attribute("IsClosed").set_value("F");
		if (def_section.isdeformed_)
			tmp_cs.append_attribute("IsDeformed").set_value("T");
		else
			tmp_cs.append_attribute("IsDeformed").set_value("F");
		for (int i = 0; i < def_pt_list.size(); i++)
		{
			pugi::xml_node ctrl_pts = tmp_cs.append_child("EmbPoint");
			ctrl_pts.append_attribute("ID").set_value(i);
			ctrl_pts.append_attribute("X").set_value(def_section.x_coordinate_);
			ctrl_pts.append_attribute("Y").set_value(def_pt_list[i].y_);
			ctrl_pts.append_attribute("Z").set_value(def_pt_list[i].z_);
		}
		std::cout << "Saving Cross-sections: " << def_doc.save_file(outfname.c_str()) << std::endl;

	}
}
void file_io::get_outxml(Basic_data_of_DesSection def_section, std::vector<Point2> def_pt_list, std::string outfname)
{
	Write_xml(def_section, def_pt_list,outfname);
}