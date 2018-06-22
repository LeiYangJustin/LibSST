#ifndef MESH_OBJECT_H
#define MESH_OBJECT_H
#include"prereq.h"
#include "custom_openmesh_type.h"
#include<vector>
#include<set>
class DATACOLLE_CLASS CMeshObject
{
protected:
	int mesh_id_;// id of mesh object
	bool is_changed_;//if the geometry of mesh is changed , it should be true, thus the opengl render will re-computeds
	std::vector<double> bbox_;
	COpenMeshT	mesh_;

public:
	CMeshObject();
	CMeshObject(CMeshObject &b);
	~CMeshObject();

	int GetId();//get id of mesh object
	void SetId(int id);//set id of mesh object
	void CopyFrom(CMeshObject& b);
	COpenMeshT& GetMesh() { return mesh_; }

	bool IsChanged();
	void SetChanged(bool is_changed=true);
};


#endif