#ifndef CUSTOM_OPENMESH_TYPE_H
#define CUSTOM_OPENMESH_TYPE_H
#include "prereq.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Traits.hh>
#include <OpenMesh/Tools/Utils/getopt.h>


struct COMTraits : public OpenMesh::DefaultTraits
{
	// store barycenter of neighbors in this member
	typedef OpenMesh::Vec3d Point; // use double-values points
	typedef OpenMesh::Vec3d Color; 
	typedef OpenMesh::Vec3d Normal;

	VertexTraits
	{
	private:
		Point emb_coord_;
		int vlabel_;
		

	public:
		VertexT() : emb_coord_(Point(0, 0, 0)), vlabel_(-1) { }
		
		const int& get_vlabel() const { return vlabel_; }
		const Point& get_emb_coord() const { return emb_coord_; }

		void set_vlabel(const int& _right) { vlabel_ = _right; }
		void set_emb_coord(const Point& p) { emb_coord_ = p; }
	};

	HalfedgeTraits
	{
	private:
		//TexCoord2D  uv_;
		//Color color_;
	public:
		HalfedgeT() {}
		//HalfedgeT() : uv_(TexCoord2D(0.0f, 0.0f)), color_(Color(0,0, 0)) { }
		//const TexCoord2D& GetUV() const { return uv_; }
		//const Color& GetColor() const { return color_; }
		//void SetUV(const TexCoord2D& _p) { uv_ = _p; }
		//void SetColor(const TexCoord2D& c) { color_ = c; }
	};

	FaceTraits
	{
	private:
		
	public:
		FaceT() {}
		
	};
	
};
class COpenMeshT :public OpenMesh::TriMesh_ArrayKernelT<COMTraits>
{
public:
	COpenMeshT()
	{
		request_vertex_colors();
		request_face_status();
		request_edge_status();
		request_vertex_status();
		request_face_normals();
		request_vertex_normals();
		request_face_colors();
	}
};
// ---------------------------------------------------------------------------
#endif