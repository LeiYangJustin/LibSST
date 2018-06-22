// Copyright 2016_9 by ChenNenglun
#include"data_io.h"
#include<vector>
#include <OpenMesh/Core/IO/MeshIO.hh>

//#include"cgal_igl_converter.h"
bool CDataIO::ReadMesh(std::string fname, CMeshObject & res_mesh_obj)
{
	if (!OpenMesh::IO::read_mesh(res_mesh_obj.GetMesh(), fname))
	{
		std::cerr << "read error\n";
		return false;
	}
	COpenMeshT &mesh = res_mesh_obj.GetMesh();
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		mesh.set_color(*viter, OpenMesh::Vec3d(0.8, 0.8, 0.8));
	}
	res_mesh_obj.SetChanged();
	return true;
}

//bool CDataIO::WriteMesh(std::string fname, CMeshObject & res_mesh_obj)
//{
//	if (!OpenMesh::IO::write_mesh(res_mesh_obj.GetMesh(), fname))
//	{
//		std::cerr << "write error\n";
//		return false;
//	}
//	return true;
//}
//
//bool CDataIO::WriteMesh(std::string fname, COpenMeshT & trimesh)
//{
//	if (!OpenMesh::IO::write_mesh(trimesh, fname))
//	{
//		std::cerr << "write error\n";
//		return false;
//	}
//	return true;
//}

bool CDataIO::ReadCurve(std::string fname, CCurveObject & res_curve_obj, double scale)
{
	std::ifstream curvefile(fname);
	std::vector<OpenMesh::Vec3d>& cc = res_curve_obj.GetCurve();

	std::vector<double> tmpData;
	if (curvefile.is_open())
	{
		double a;
		while (curvefile >> a)
		{
			tmpData.push_back(a);
		}
	}
	else
	{
		std::cerr << "Error in using CDataIO::ReadCurve: File not open" << std::endl;
		return false;
	}

	if (tmpData.size() % 3 == 0)
	{
		for (int i = 0; i < tmpData.size(); i += 3)
		{
			double a, b, c;
			a = tmpData[i]* scale;
			b = tmpData[i + 1]* scale;
			c = tmpData[i + 2]* scale;
			cc.push_back(OpenMesh::Vec3d(a, b, c));
		}
	}else
	{
		std::cerr << "Error in using CDataIO::ReadCurve: Data dim does not match" << std::endl;
		return false;
	}

	return true;
}
