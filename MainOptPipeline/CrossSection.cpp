#include"CrossSection.h"

void CrossSection::set_rectangular_cross_section(double length, double width)
{
	length_ = length;
	width_ = width;

	for (int i = 0; i < num_of_edgepts_ * 4; i++)
	{
		double y, z;
		if (i <= num_of_edgepts_)
		{
			z = length_ / 2;
			y = width_ / 2 - (width_ / num_of_edgepts_)*i;
		}
		else if (i > num_of_edgepts_ && i <= num_of_edgepts_ * 2)
		{
			z = length_ / 2 - (length_ / num_of_edgepts_)*(i - num_of_edgepts_);
			y = -width_ / 2;
		}
		else if (i > num_of_edgepts_ * 2 && i <= num_of_edgepts_ * 3)
		{
			z = -length_ / 2;
			y = -width_ / 2 + (width_ / num_of_edgepts_)*(i - num_of_edgepts_ * 2);
		}
		else
		{
			z = -length_ / 2 + (length_ / num_of_edgepts_)*(i - num_of_edgepts_ * 3);
			y = width_ / 2;
		}
		Point2 tmp_pt(z,y);
		src_pt_list_.push_back(tmp_pt);
	}
}
double CrossSection::get_length() 
{ 
	return length_; 
}
double CrossSection::get_width() 
{ 
	return width_; 
}
void CrossSection::set_the_number_of_edgepts(int ptsid) 
{ 
	num_of_edgepts_ = ptsid; 
}
int CrossSection::get_the_number_of_edgepts() 
{ 
	return num_of_edgepts_; 
}
void CrossSection::set_ctrl_pid_list_and_def_list(std::vector<int> ctrl_pid_list, std::vector<double> ctrl_pts_def_list)
{
	if (ctrl_pid_list.size() == ctrl_pts_def_list.size())
	{
		ctrl_pid_list_ = ctrl_pid_list;
		ctrl_pts_def_list_ = ctrl_pts_def_list;
		//for (int i = 0; i < edge_ctrl_pts_def_array.size(); i++)
		//{
		//	ctrl_pts_def_list_.insert(ctrl_pts_def_list_.end(), edge_ctrl_pts_def_array[i].begin(), edge_ctrl_pts_def_array[i].end());
		//}
		// make a closed loop
		ctrl_pid_list_.push_back(num_of_edgepts_*4);
		ctrl_pts_def_list_.push_back(0);
	}
	else {
		return;
	}
}
std::vector<Point2> CrossSection::get_src_pt_list()
{
	return src_pt_list_;
}
//void CrossSection::deformation_algorithm()
//{
//	std::vector<int> ctrlpts_stride_list;
//	std::vector<double> def_list;
//	for (int i = 1; i < ctrl_pid_list_.size(); i++)
//	{
//		int tmp_stride = ctrl_pid_list_[i] - ctrl_pid_list_[i - 1];
//		ctrlpts_stride_list.push_back(tmp_stride);//得到相邻两个控制点之间的点的数量
//	}
//	
//	for (int j = 0; j < 4; j++)//循环不同的边
//	{
//		for (int i = 0; i < ctrl_pid_list_.size(); i++)//每条边上不同的控制点
//		{
//			for (int pid = ctrl_pid_list_[i]; pid < ctrl_pid_list_[i + 1]; pid++)//相邻两个控制点之间点的ID
//			{
//				double tmp_def = edge_ctrl_pts_def_array_[j][i] 
//					- ((edge_ctrl_pts_def_array_[j][i] - edge_ctrl_pts_def_array_[j][i + 1]) / ctrlpts_stride_list[i])*(pid - ctrl_pid_list_[i]);
//				def_list.push_back(tmp_def);//得到每个点的变形量
//			}
//		}
//	}
//	def_pt_list_ = src_pt_list_;
//	for (int i = 0; i < def_pt_list_.size(); i++)
//	{
//		if (i < num_of_edgepts_)//位于一四象限的垂直于z轴的边（用y轴，z轴组成二维坐标系），依次顺时针选取边
//		{
//			def_pt_list_[i].y_ += def_list[i];
//		}
//		else if (i>= num_of_edgepts_ && i<num_of_edgepts_ *2)
//		{
//			def_pt_list_[i].x_ -= def_list[i];
//		}
//		else if (i >= num_of_edgepts_ *2 && i<num_of_edgepts_ * 3)
//		{
//			def_pt_list_[i].y_ -= def_list[i];
//		}
//		else if (i >= num_of_edgepts_ *3 && i<num_of_edgepts_ * 4)
//		{
//			def_pt_list_[i].x_ += def_list[i];
//		}
//	}
//	is_deformed_ = true;
//}
void CrossSection::Deformation_algorithm()
{
	// 
	std::vector<double> def_list;
	for (int i = 0; i < ctrl_pid_list_.size()-1; i++)
	{
		int ctrl_pt1 = ctrl_pid_list_[i];
		int ctrl_pt2 = ctrl_pid_list_[i+1];
		double def_at_ctrl_pt1 = ctrl_pts_def_list_[i];
		double def_at_ctrl_pt2 = ctrl_pts_def_list_[i+1];
		// do interpolation
		for (int j = ctrl_pt1; j < ctrl_pt2; j++)
		{
			// do interpolation(j, ctrl_pt1, ctrl_pt2, def_at_ctrl_pt1, def_at_ctrl_pt2)
			double tmp_def = def_at_ctrl_pt1 - ((def_at_ctrl_pt1 - def_at_ctrl_pt2) / (ctrl_pt2 - ctrl_pt1))*(j - ctrl_pt1);
			def_list.push_back(tmp_def);
		}
	}
	def_pt_list_ = src_pt_list_;
	for (int i = 0; i < def_pt_list_.size(); i++)
	{
		if (i < num_of_edgepts_)//位于一四象限的垂直于z轴的边（用y轴，z轴组成二维坐标系），依次顺时针选取边
		{
			def_pt_list_[i].z_ += def_list[i];
		}
		else if (i >= num_of_edgepts_ && i < num_of_edgepts_ * 2)
		{
			def_pt_list_[i].y_ -= def_list[i];
		}
		else if (i >= num_of_edgepts_ * 2 && i < num_of_edgepts_ * 3)
		{
			def_pt_list_[i].z_ -= def_list[i];
		}
		else if (i >= num_of_edgepts_ * 3 && i < num_of_edgepts_ * 4)
		{
			def_pt_list_[i].y_ += def_list[i];
		}
	}
	is_deformed_ = true;
}
void CrossSection::deform()
{
	Deformation_algorithm();
}
std::vector<Point2> CrossSection::get_def_pt_list()
{
	if (is_deformed_ == true)
	{
		return def_pt_list_;
	}
	else
	{
		Deformation_algorithm();
		return def_pt_list_;
	}
}
