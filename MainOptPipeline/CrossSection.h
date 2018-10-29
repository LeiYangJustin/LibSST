#ifndef C_CROSS_SECTION_H
#define C_CROSS_SECTION_H
#include <vector>
#include <iostream>
#include"General_function.h"
class CrossSection
{
private:
	double length_;//��
	double width_;//��
	int num_of_edgepts_;//һ�����ж��ٸ���

	std::vector<double> ctrl_pts_def_list_; //���Ƶ��Ӧ�ı���
	std::vector<int> ctrl_pid_list_;//���Ƶ��ID�б�, 0 6 12 18 ..
	std::vector<Point2> src_pt_list_;//���е�����꣨����ǰ��
	std::vector<Point2> def_pt_list_;//���е�����꣨���κ�
	bool is_deformed_;//�Ƿ�����˱���

	void Deformation_algorithm();//�����㷨
	
public:
	CrossSection() {};//���캯��
	~CrossSection() {};//��������
	void set_rectangular_cross_section(double length, double width);//��ʼ�������
	void set_the_number_of_edgepts(int num_of_edgepts);//��ʼ��һ�ߵĵ����Ŀ
	void set_ctrl_pid_list_and_def_list(std::vector<int> ctrl_pid_list, std::vector<double> ctrl_pts_def_list);//���ÿ��Ƶ�ID����Ӧ������
	void deform();//���ñ����㷨

	double get_length();//��ȡ��
	double get_width();//��ȡ��
	int get_the_number_of_edgepts();//��ȡһ�ߵĵ����Ŀ
	std::vector<Point2> get_def_pt_list();//��ȡ���κ�������ֵ�������
	std::vector<Point2> get_src_pt_list();//��ȡ�����������ֵ��ĳ�ʼ����


};



#endif // !C_CROSS_SECTION_H

