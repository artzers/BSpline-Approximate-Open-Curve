#pragma once 
#include <Eigen/Dense>
#include <vector>
#include "../include/types.h"

//struct Point {
//	double x, y,z;
//	Point(const double &a, const double &b, const double &c) {
//		x = a;
//		y = b;
//		z = c;
//	}
//};

struct BSpline {
	BSpline CubicApproximate(const std::vector<Point>& vecFitPoints, double alpha, double beta);

	void WMatrix(Eigen::MatrixXd & W, int n, int p, const std::vector<double>& u);

	void _SecondDerivativeCoefficient(int i, const std::vector<double>& u, std::vector<std::pair<double, double>>& a_b_array);

	double PolynomialIntegral(double quad, double linear, double con, double start, double end);

	void BesselTanget(const Point & p0, const Point & p1, const Point & p2, Vec3d & p0deriv, Vec3d & p1deriv, Vec3d & p2deriv);

	void BasisFunc(double u, int k, std::vector<double>& basis_func);

	Point Evaluate(double u);

	void InsertKnot(double u);

	double PointDistance(const Point & arg1, const Point & arg2);

	void Tesselation(double tolerance, std::vector<Point>& points);

	void _tesselation(std::vector<Point>& cv, std::vector<double>& knots, std::vector<Point>& points, double tolerance, int k_s, int k_e, int c_s, int c_e);

	void Subdeviding(double u, BSpline & sub_left, BSpline & sub_right);

	void _subdeviding(double u, int k, const std::vector<double>& knots, const std::vector<Point>& cv, std::vector<Point>& cv_left, std::vector<Point>& cv_right);

	double PointLineDistance(const Point & arg1, const Point & argS, const Point & argE);

	BSpline CubicInterpolate(const std::vector<Point>& vecFitPoints);



	int m_nDegree;
	std::vector<double> m_vecKnots;
	std::vector<Point> m_vecCVs;

};