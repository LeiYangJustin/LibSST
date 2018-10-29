#ifndef C_CROSS_SECTION_H
#define C_CROSS_SECTION_H
#include <vector>
#include <iostream>
#include"General_function.h"
class CrossSection
{
private:
	double length_;//长
	double width_;//宽
	int num_of_edgepts_;//一边上有多少个点

	std::vector<double> ctrl_pts_def_list_; //控制点对应的变形
	std::vector<int> ctrl_pid_list_;//控制点的ID列表, 0 6 12 18 ..
	std::vector<Point2> src_pt_list_;//所有点的坐标（变形前）
	std::vector<Point2> def_pt_list_;//所有点的坐标（变形后）
	bool is_deformed_;//是否进行了变形

	void Deformation_algorithm();//变形算法
	
public:
	CrossSection() {};//构造函数
	~CrossSection() {};//析构函数
	void set_rectangular_cross_section(double length, double width);//初始化横截面
	void set_the_number_of_edgepts(int num_of_edgepts);//初始化一边的点的数目
	void set_ctrl_pid_list_and_def_list(std::vector<int> ctrl_pid_list, std::vector<double> ctrl_pts_def_list);//设置控制点ID和相应变形量
	void deform();//调用变形算法

	double get_length();//获取长
	double get_width();//获取宽
	int get_the_number_of_edgepts();//获取一边的点的数目
	std::vector<Point2> get_def_pt_list();//获取变形后所有型值点的坐标
	std::vector<Point2> get_src_pt_list();//获取截面的所有型值点的初始坐标


};



#endif // !C_CROSS_SECTION_H

