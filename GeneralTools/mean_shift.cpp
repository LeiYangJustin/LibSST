#include <stdio.h>
#include <math.h>
#include "mean_shift.h"


//#define EPSILON 0.00000001
#define EPSILON 0.00001
#define CLUSTER_EPSILON 0.5

double euclidean_distance(const std::vector<double> &point_a, const std::vector<double> &point_b) {
	double total = 0;
	for (int i = 0; i<point_a.size(); i++) {
		total += (point_a[i] - point_b[i]) * (point_a[i] - point_b[i]);
	}
	return sqrt(total);
}

double gaussian_kernel(double distance, double kernel_bandwidth) {
	double temp = exp(-1.0 / 2.0 * (distance*distance) / (kernel_bandwidth*kernel_bandwidth));
	return temp;
}

//void MeanShift::set_kernel(double(*_kernel_func)(double, double)) {
//	if (!_kernel_func) {
//		kernel_func = gaussian_kernel;
//	}
//	else {
//		kernel_func = _kernel_func;
//	}
//}


// private member
std::vector<double> CMeanShift::shift_point(
	const std::vector<double> &point, 
	const std::vector<std::vector<double> > &points, 
	double kernel_bandwidth) 
{
	std::vector<double> shifted_point = point;
	for (int dim = 0; dim<shifted_point.size(); dim++) {
		shifted_point[dim] = 0;
	}
	double total_weight = 0;
	for (int i = 0; i<points.size(); i++) {
		std::vector<double> temp_point = points[i];
		double distance = euclidean_distance(point, temp_point);
		double weight = gaussian_kernel(distance, kernel_bandwidth);
		for (int j = 0; j<shifted_point.size(); j++) {
			shifted_point[j] += temp_point[j] * weight;
		}
		total_weight += weight;
	}

	for (int i = 0; i<shifted_point.size(); i++) {
		shifted_point[i] /= total_weight;
	}
	return shifted_point;
}

std::vector<std::vector<double> > CMeanShift::meanshift(const std::vector<std::vector<double> > & points, double kernel_bandwidth) 
{
	std::vector<bool> stop_moving(points.size(), false);
	std::vector<std::vector<double> > shifted_points = points;
	double max_shift_distance;
	int iter_cnt = 0;
	do {
		max_shift_distance = 0;
		for (int i = 0; i<shifted_points.size(); i++) {
			if (!stop_moving[i]) {
				std::vector<double>point_new = shift_point(shifted_points[i], points, kernel_bandwidth);
				double shift_distance = euclidean_distance(point_new, shifted_points[i]);
				//
				if (shift_distance > max_shift_distance) {
					max_shift_distance = shift_distance;
				}
				if (shift_distance <= EPSILON) {
					stop_moving[i] = true;
				}
				shifted_points[i] = point_new;
			}
		}
		if (iter_cnt%100 == 99)
			printf("iter: %d, max_shift_distance: %f\n", iter_cnt, max_shift_distance);

		iter_cnt++;
	} while (max_shift_distance > EPSILON);
	printf("iteration: %d\n", iter_cnt);
	return shifted_points;
}


std::vector<MSCluster> CMeanShift::cluster(
	const std::vector<std::vector<double>> & points,
	const std::vector<std::vector<double>> & shifted_points)
{
	std::vector<MSCluster> clusters;

	for (int i = 0; i < shifted_points.size(); i++) {

		// if distance is smaller than CLUSTER_EPSILON, this shifted point is added to the cluster mode
		int c = 0;
		for (; c < clusters.size(); c++) {
			if (euclidean_distance(shifted_points[i], clusters[c].mode) <= CLUSTER_EPSILON) {
				break;
			}
		}

		// else if the shifited point belongs to no modes, add it as a new mode
		if (c == clusters.size()) {
			MSCluster clus;
			clus.mode = shifted_points[i];
			clusters.push_back(clus);
		}

		clusters[c].original_points.push_back(points[i]);
		clusters[c].shifted_points.push_back(shifted_points[i]);
	}

	return clusters;
}

std::vector<MSCluster> CMeanShift::cluster(const std::vector<std::vector<double>> & points, double kernel_bandwidth) 
{
	std::vector<std::vector<double> > shifted_points = meanshift(points, kernel_bandwidth);
	return cluster(points, shifted_points);
}

void CMeanShift::GetPID2Cluster(
	const std::vector<std::vector<double>> & shifted_points, 
	std::vector<std::vector<int>> & pids2clusters, double cluster_eplison)
{
	std::vector<MSCluster> clusters;

	for (int i = 0; i < shifted_points.size(); i++) {

		// if distance is smaller than CLUSTER_EPSILON, this shifted point is added to the cluster mode
		int c = 0;
		for (; c < clusters.size(); c++) {
			if (euclidean_distance(shifted_points[i], clusters[c].mode) <= cluster_eplison) {
				break;
			}
		}

		// else if the shifited point belongs to no modes, add it as a new mode
		if (c == clusters.size()) {
			MSCluster clus;
			clus.mode = shifted_points[i];
			clusters.push_back(clus);
			//
			std::vector<int> p;
			pids2clusters.push_back(p);
		}

		pids2clusters[c].push_back(i);
		clusters[c].shifted_points.push_back(shifted_points[i]);
	}
}
