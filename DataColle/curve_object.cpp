#include"curve_object.h"
int CCurveObject::GetId()
{
	return curve_id_;
}
void CCurveObject::SetId(int id)
{
	curve_id_ = id;
}
bool CCurveObject::IsChanged()
{
	return is_changed_;
}
void CCurveObject::SetChanged(bool is_changed)
{
	is_changed_ = is_changed;
}
std::vector<OpenMesh::Vec3d>& CCurveObject::GetCurve()
{
	return curve_;
}
CCurveObject::CCurveObject()
{
	curve_id_ = -1;
}
CCurveObject::CCurveObject(CCurveObject&b)
{
	curve_id_ = b.curve_id_;
	curve_ = b.curve_;
	is_changed_ = true;
}
CCurveObject::~CCurveObject()
{

}