/*
 *  main.cpp
 *  SplineFitting
 *  This file shows an example usage of the spline fitting method.
 *  This example produces a smooth sinusoid curve with just 6 anchor points.
 */
#pragma warning(disable:4996)
#include <iostream>
#include <iterator>
#include <iomanip>
#include <ostream>
#include <vector>
#include "include/BSpline.h"
using namespace std;

//1.0, 1.5, 3.7, 4.3, 5.6, 6.1
//0.0, 1.3, 0.5, -1.0, 2.1, 3.5
int main() {
	std::vector<Point> fitPt;
	fitPt.push_back(Point(0.0, 0.0, 0.0));
	fitPt.push_back(Point(1.0, 0.0, 0.0));
	fitPt.push_back(Point(1.0, 1.0, 0.0));
	fitPt.push_back(Point(0.5, 2.0, 0.0));
	fitPt.push_back(Point(0.0, 4.0, 0.0));
	fitPt.push_back(Point(-1.0, 6.0, 0.0));
	/*fitPt.push_back(Point(1.0, 0.0, 0.0));
	fitPt.push_back(Point(1.5, 1.3, 0.0));
	fitPt.push_back(Point(3.7, 0.5, 0.0));
	fitPt.push_back(Point(4.3, -1.0, 0.0));
	fitPt.push_back(Point(5.6, 2.1, 0.0));
	fitPt.push_back(Point(6.1, 3.5, 0.0));*/
	BSpline bs1;
	//BSpline bs2 = bs1.CubicInterpolate(fitPt);
	BSpline bs2 = bs1.CubicApproximate(fitPt,0.02,1);

	std::vector<Point> curve;
	bs2.Tesselation(0.05, curve);
	FILE *fp = fopen("hehe1.csv","w");
	for (size_t k = 0; k < curve.size(); ++k) {
		fprintf(fp, "%lf\n", curve[k](0));
	}
	fclose(fp);
	fp = fopen("hehe2.csv","w");
	for (size_t k = 0; k < curve.size(); ++k) {
		fprintf(fp, "%lf\n", curve[k](1));
	}
	fclose(fp);

	return 0;
}