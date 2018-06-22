// Copyright 2016_9 by ChenNenglun
#ifndef DATA_IO_H
#define DATA_IO_H
#include"prereq.h"
#include<iostream>
#include<fstream>
#include<string>
#include"mesh_object.h"
#include"curve_object.h"

class DATACOLLE_CLASS CDataIO
{
protected:

public:
	

	static bool ReadMesh(std::string fname, CMeshObject & res_mesh_obj);
	//static bool WriteMesh(std::string fname, CMeshObject & res_mesh_obj);
	//static bool WriteMesh(std::string fname, COpenMeshT & trimesh);
	

	static bool ReadCurve(std::string fname, CCurveObject & res_curve_obj, double scale = 1.0);
};
#endif
