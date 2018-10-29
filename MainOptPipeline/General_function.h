#ifndef GENERAL_FUNCTION_H
#define GENERAL_FUNCTION_H

class Point2
{
public:
	Point2(double z, double y) { z_ = z; y_ = y; };//以原坐标的z轴坐标为Point2的x坐标
	~Point2() {};
	double z_;
	double y_;

	// overload
	Point2 operator+(Point2 &add)//重载运算符，当两个Point2相加时，使其内部的z_和y_相加，并返回带有z_和y_最终值的Point2
	{
		z_ += add.z_;
		y_ += add.y_;
		return Point2(z_, y_);
	}
};
class Basic_data_of_DesSection
{
public:
	Basic_data_of_DesSection() {};
	~Basic_data_of_DesSection() {};
	bool isclosed_;
	bool isdeformed_;
	int id_;//截面的ID
	int spid_;//截面的SPID
	double x_coordinate_;//x坐标

	double length_;//长
	double width_;//宽
	int pts_num_on_single_edge;//截面上一边的型值点数目（取边的末端点的id）
	std::vector<double> ctrl_pt_def_list_;//该截面上全部控制点的变形量
	std::vector<int> ctrl_pid_list_;//该截面上全部控制点的ID（对应的型值点id）
};




#endif // !GENERAL_FUNCTION_H

