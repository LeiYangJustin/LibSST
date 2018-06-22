#ifndef C_MEAN_SHIFT_CLASS_H
#define C_MEAN_SHIFT_CLASS_H

// this code is modified from https://github.com/mattnedrich/MeanShift_cpp

#include <vector>
#include "general_tool_prereq.h"

struct MSCluster {
	std::vector<double> mode;
	std::vector<std::vector<double> > original_points;
	std::vector<std::vector<double> > shifted_points;
};

class GENERAL_TOOLS_CLASS CMeanShift {
public:
	//MeanShift() { set_kernel(NULL); }
	//MeanShift(double(*_kernel_func)(double, double)) { set_kernel(kernel_func); }
	CMeanShift() {};
	std::vector<std::vector<double> > meanshift(const std::vector<std::vector<double>> & points, double);
	std::vector<MSCluster> cluster(const std::vector<std::vector<double>> & points, double kernel_bandwidth);
	void GetPID2Cluster(const std::vector<std::vector<double>> & shifted_points, std::vector<std::vector<int>> & pids2clusters, double cluster_epsilon);

private:
	double(*kernel_func)(double, double);
	//void set_kernel(double(*_kernel_func)(double, double));
	std::vector<double> shift_point(const std::vector<double> & point, const std::vector<std::vector<double> > & points, double kernel_bandwidth);
	std::vector<MSCluster> cluster(const std::vector<std::vector<double> > & points, const std::vector<std::vector<double> > & shifted_points);
};

#endif // C_MEAN_SHIFT_CLASS_H