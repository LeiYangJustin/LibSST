#ifndef GENERAL_FUNCTION_H
#define GENERAL_FUNCTION_H

class Point2
{
public:
	Point2(double z, double y) { z_ = z; y_ = y; };//��ԭ�����z������ΪPoint2��x����
	~Point2() {};
	double z_;
	double y_;

	// overload
	Point2 operator+(Point2 &add)//�����������������Point2���ʱ��ʹ���ڲ���z_��y_��ӣ������ش���z_��y_����ֵ��Point2
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
	int id_;//�����ID
	int spid_;//�����SPID
	double x_coordinate_;//x����

	double length_;//��
	double width_;//��
	int pts_num_on_single_edge;//������һ�ߵ���ֵ����Ŀ��ȡ�ߵ�ĩ�˵��id��
	std::vector<double> ctrl_pt_def_list_;//�ý�����ȫ�����Ƶ�ı�����
	std::vector<int> ctrl_pid_list_;//�ý�����ȫ�����Ƶ��ID����Ӧ����ֵ��id��
};




#endif // !GENERAL_FUNCTION_H

