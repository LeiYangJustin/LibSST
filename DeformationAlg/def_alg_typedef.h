#ifndef C_DEFORM_TYPE_H
#define C_DEFORM_TYPE_H

#include "def_alg_prereq.h"
#include <iostream>

#include<Eigen/Dense>
#include<Eigen/Sparse>

typedef Eigen::Vector3d					Point3d;
typedef Eigen::Vector3d					Vector3d;
typedef Eigen::VectorXd					VectorXd;
typedef Eigen::Matrix<double, 9, 1>		Tensor3d;
typedef Eigen::MatrixXd					DenseMatrixXd;
typedef Eigen::MatrixXi					DenseMatrixXi;
typedef Eigen::Matrix3d					Mat2d;
typedef Eigen::Matrix3d					Mat3d;
typedef Eigen::Matrix4d					Mat4d;
typedef Eigen::SparseMatrix<double>		SparseMatrixXd;
typedef Eigen::SparseMatrix<int>		SparseMatrixXi;


#endif // !C_DEFORM_TYPE_H
