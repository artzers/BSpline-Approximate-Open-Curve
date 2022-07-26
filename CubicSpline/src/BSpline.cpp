#include <Eigen/Dense> //Eigen
#include <iostream>
#include "..\include\BSpline.h"


using namespace std;
using namespace Eigen;
/*!
 *\brief 三次B样条整体光顺逼近
*\ param const std::vector<Point> & vecFitPoints 型值点
*\ param double alpha 光顺项权重
*\ param double beta  逼近项权重
*\ Returns:   BSpline 逼近结果
*/
BSpline BSpline::CubicApproximate(const std::vector<Point>& vecFitPoints, double alpha, double beta)
{
	const int p = 3;
	BSpline bs;
	int  x = vecFitPoints.size();
	if (x < p)
	{
		cout << "too less point !" << endl;
		return bs;
	}
	//需要的矩阵
	Eigen::MatrixXd W = Eigen::MatrixXd::Zero(x + 2, x + 2);
	Eigen::MatrixXd P = Eigen::MatrixXd::Zero(x + 2, 3);
	Eigen::MatrixXd M = Eigen::MatrixXd::Zero(x + 2, x + 2);
	Eigen::MatrixXd F = Eigen::MatrixXd::Zero(x + 2, 3);
	//参数化
	bs.m_nDegree = p;
	bs.m_vecKnots.resize(x); //x+6个节点
	//计算节点
	bs.m_vecKnots[0] = 0.0;
	for (int i = 1; i < x; ++i)
	{
		bs.m_vecKnots[i] = bs.m_vecKnots[i - 1]
			+ sqrt(PointDistance(vecFitPoints[i], vecFitPoints[i - 1]));
	}
	//节点首尾构成p+1度重复
	bs.m_vecKnots.insert(bs.m_vecKnots.begin(), p, bs.m_vecKnots.front());
	bs.m_vecKnots.insert(bs.m_vecKnots.end(), p, bs.m_vecKnots.back());

	//W 矩阵
	WMatrix(W, x + 2, p, bs.m_vecKnots);

	//M矩阵
	std::vector<double> basis_func;
	M(0, 0) = 1;
	M(x - 1, x + 1) = 1;
	for (int i = p + 1; i < x + p - 1; ++i)
	{
		//c(u)在 N_{i-p},...,N_i等p+1个基函数上非零
		bs.BasisFunc(bs.m_vecKnots[i], i, basis_func);
		for (int j = i - p, k = 0; j <= i; ++j, ++k)
		{
			M(i - p, j) = basis_func[k];
		}
	}
	//导数
	M(x, 0) = -1;
	M(x, 1) = 1;
	M(x + 1, x) = -1;
	M(x + 1, x + 1) = 1;

	//F矩阵
	for (int i = 0; i < x; ++i)
	{
		F(i, 0) = vecFitPoints[i](0);
		F(i, 1) = vecFitPoints[i](1);
		F(i, 2) = vecFitPoints[i](2);
	}

	{
		Vec3d v0, v1, v2;
		BesselTanget(vecFitPoints[0], vecFitPoints[1], vecFitPoints[2], v0, v1, v2);
		Vec3d v = v0 * (bs.m_vecKnots[p + 1] - bs.m_vecKnots[1]) / (double)p;
		F(x, 0) = v(0);
		F(x, 1) = v(1);
		F(x, 2) = v(2);
	}

	{
		Vec3d v0, v1, v2;
		BesselTanget(vecFitPoints[x - 3], vecFitPoints[x - 2], vecFitPoints[x - 1], v0, v1, v2);
		Vec3d v = v2 * (bs.m_vecKnots[x + 1 + p] - bs.m_vecKnots[x + 1]) / (double)p;
		F(x + 1, 0) = v(0);
		F(x + 1, 1) = v(1);
		F(x + 1, 2) = v(2);
	}

	//解方程
	P = (alpha*W + beta * M.transpose()*M).colPivHouseholderQr().solve(beta*M.transpose()*F);

#ifdef _DEBUG
	cout << "P------------------" << endl << P << endl;
	cout << "W------------------" << endl << W << endl;
	cout << "M------------------" << endl << M << endl;
	cout << "F------------------" << endl << F << endl;
#endif

	//将 P的值赋给 B样条
	bs.m_vecCVs.resize(x + 2);
	for (int i = 0; i < x + 2; ++i)
	{
		Point& cv = bs.m_vecCVs[i];
		cv(0) = P(i, 0);
		cv(1) = P(i, 1);
		cv(2) = P(i, 2);
	}

	return bs;
}

void BSpline::WMatrix(Eigen::MatrixXd& W, int n, int p, const std::vector<double>& u)
{
	std::vector<std::vector<std::pair<double, double>>> BasisFuncByKnot(n);
	//初始化
	for (int i = 0; i < n; ++i)
	{
		BasisFuncByKnot[i].resize(n + p);
		for (int j = 0; j < n + p; ++j)
		{
			BasisFuncByKnot[i][j].first = 0;
			BasisFuncByKnot[i][j].second = 0;
		}
	}

	std::vector<std::pair<double, double>> a_b_array;
	for (int i = 0; i < n; ++i)
	{
		_SecondDerivativeCoefficient(i, u, a_b_array);
		std::copy(a_b_array.begin(), a_b_array.end(), BasisFuncByKnot[i].begin() + i);
	}
	//基函数//其实这是一个对称矩阵，我为了方便，多算了一倍
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			double ret = 0;
			//区间
			for (int k = 0; k < n + p; ++k)
			{
				const std::pair<double, double> basis_i = BasisFuncByKnot[i][k];
				const std::pair<double, double> basis_j = BasisFuncByKnot[j][k];

				ret += PolynomialIntegral(basis_i.first * basis_j.first, basis_i.first*basis_j.second +
					basis_i.second*basis_j.first, basis_i.second*basis_j.second, u[k], u[k + 1]);
			}
			W(i, j) = ret;
		}
	}
}

void BSpline::_SecondDerivativeCoefficient(int i, const std::vector<double>& u, std::vector<std::pair<double, double>>& a_b_array)
{
	const int p = 3;
	double c1, c2, c3;
	c1 = c2 = c3 = 0;

	double div = (u[i + p - 1] - u[i])*(u[i + p] - u[i]);
	if (div != 0)c1 = p * (p - 1) / div;

	div = (u[i + p] - u[i])*(u[i + p] - u[i + 1]);
	if (div != 0)c2 -= p * (p - 1) / div;
	div = (u[i + p] - u[i + 1])*(u[i + p + 1] - u[i + 1]);
	if (div != 0)c2 -= p * (p - 1) / div;

	div = (u[i + p + 1] - u[i + 1])*(u[i + p + 1] - u[i + 2]);
	if (div != 0)c3 = p * (p - 1) / div;

	a_b_array.resize(p + 1);
	for (int i = 0; i < p + 1; ++i)
	{
		a_b_array[i].first = 0;
		a_b_array[i].second = 0;
	}

	div = u[i + p - 2] - u[i];
	if (c1 != 0 && div != 0)
	{
		a_b_array[0].first = c1 / div;
		a_b_array[0].second = -c1 * u[i] / div;
	}

	div = u[i + p - 1] - u[i + 1];
	if (div != 0)
	{
		a_b_array[1].first = (c2 - c1) / div;
		a_b_array[1].second = (c1*u[i + p - 1] - c2 * u[i + 1]) / div;
	}

	div = u[i + p] - u[i + 2];
	if (div != 0)
	{
		a_b_array[2].first = (c3 - c2) / div;
		a_b_array[2].second = (c2*u[i + p] - c3 * u[i + 2]) / div;
	}

	div = u[i + p + 1] - u[i + 3];
	if (c3 != 0 && div != 0)
	{
		a_b_array[3].first = -c3 / div;
		a_b_array[3].second = u[i + p + 1] * c3 / div;
	}
}

double BSpline::PolynomialIntegral(double quad, double linear, double con, double start, double end)
{
	if (end == start)
		return 0;

	double ret = 0;
	if (quad != 0)
	{
		ret += (end*end*end - start * start*start)*quad / 3;
	}
	if (linear != 0)
	{
		ret += (end*end - start * start)*linear / 2;
	}
	if (con != 0)
	{
		ret += (end - start)*con;
	}

	return ret;
}

void BSpline::BesselTanget(const Point& p0, const Point& p1, const Point& p2,
	Vec3d& p0deriv, Vec3d& p1deriv, Vec3d& p2deriv)
{
	double delta_t1 = PointDistance(p2, p1);
	double delta_t0 = PointDistance(p1, p0);
	Vec3d delta_p0 = (p1 - p0) / delta_t0;
	Vec3d delta_p1 = (p2 - p1) / delta_t1;
	double sum = delta_t0 + delta_t1;
	p1deriv = delta_t1 / sum * delta_p0 + delta_t0 / sum * delta_p1;
	p0deriv = 2 * delta_p0 - p1deriv;
	p2deriv = 2 * delta_p1 - p1deriv;
}

void BSpline::BasisFunc(double u, int k, std::vector<double>& basis_func)
{
	const int& p = m_nDegree;
	const std::vector<double>& knots = m_vecKnots;
	basis_func.resize(p + 1);

	//2p+1个 N_{i,0}
	int n = 2 * p + 1;
	vector<double> temp(n, 0);
	temp[p] = 1;

	//迭代p次
	for (int j = 1; j <= p; ++j)
	{
		//区间 [k-p,k+p+1)
		for (int i = k - p, h = 0; h < (n - j); ++h, ++i)
		{
			//递推公式
			double a = (u - knots[i]);
			double dev = (knots[i + j] - knots[i]);
			a = (dev != 0) ? a / dev : 0;

			double b = (knots[i + j + 1] - u);
			dev = (knots[i + j + 1] - knots[i + 1]);
			b = (dev != 0) ? b / dev : 0;

			temp[h] = a * temp[h] + b * temp[h + 1];
		}
	}
	//拷贝前 p+1个值到basis_func
	std::copy(temp.begin(), temp.begin() + p + 1, basis_func.begin());
}

Point BSpline::Evaluate(double u)
{
	const std::vector<double>& knots = m_vecKnots;
	int p = m_nDegree;

	//查找 u 所在区间
	int k = std::distance(knots.begin(), std::upper_bound(knots.begin(), knots.end(), u)) - 1;

	//申请 p+1个point,分别存放 P_{k-p},...,P_k
	Point *pArray = new Point[p + 1];
	std::copy(m_vecCVs.begin() + k - p, m_vecCVs.begin() + k + 1, pArray);

	//迭代 p次
	for (int r = 1; r <= p; ++r)
	{
		//i 从 k 到 k-p+1
		for (int i = k, j = p; i >= k - p + r; --i, --j)
		{
			double alpha = u - knots[i];
			double dev = (knots[i + p + 1 - r] - knots[i]);
			alpha = (dev != 0) ? alpha / dev : 0;

			pArray[j] = (1.0 - alpha)*pArray[j - 1] + alpha * pArray[j];
		}
	}

	//所求点是最后一个
	Point pt = pArray[p];
	delete[]pArray;
	return pt;
}

void BSpline::InsertKnot(double u)
{
	std::vector<double>& knots = m_vecKnots;
	int p = m_nDegree;

	//查找 u 所在区间
	int k = std::distance(knots.begin(), std::upper_bound(knots.begin(), knots.end(), u)) - 1;
	Point *pArray = new Point[p];

	//i \in [k,k-p+1]
	//j为数组pArray位置
	for (int i = k, j = p - 1; i > k - p; --i, --j)
	{
		double alpha = (u - knots[i]);
		double dev = (knots[i + p] - knots[i]);
		alpha = (dev != 0) ? alpha / dev : 0;

		pArray[j] = m_vecCVs[i - 1] * (1 - alpha) + m_vecCVs[i] * alpha;
	}

	//将cv [k-p+1,k-1]替换为pArray，并且在cv[k]之前插入 pArray[p-1]
	for (int i = k - p + 1, j = 0; i < k; ++i, ++j)
	{
		m_vecCVs[i] = pArray[j];
	}
	m_vecCVs.insert(m_vecCVs.begin() + k, pArray[p - 1]);

	//knots 插入u
	knots.insert(knots.begin() + k + 1, u);

	delete[] pArray;
}

double BSpline::PointDistance(const Point &arg1, const Point &arg2)
{
	return std::sqrt(std::pow(arg1(0) - arg2(0),2)
		+ std::pow(arg1(1) - arg2(1), 2)
		+ std::pow(arg1(2) - arg2(2), 2));
}

void BSpline::Tesselation(double tolerance, std::vector<Point>& points)
{
	std::vector<Point> cv_copy = m_vecCVs;
	std::vector<double> knots_copy = m_vecKnots;
	_tesselation(cv_copy, knots_copy, points, tolerance, m_nDegree, m_vecKnots.size() - m_nDegree, 0, m_vecCVs.size() - 1);
}

void BSpline::_tesselation(std::vector<Point>& cv, std::vector<double>& knots, std::vector<Point>& points, double tolerance, int k_s, int k_e, int c_s, int c_e)
{
	int p = m_nDegree;

	bool stright_enough = true;
	//计算弦高是否超过容差
	for (int i = c_s + 1; i < c_e; ++i)
	{
		//点到直线距离公式，不给出实现了
		if (PointLineDistance(cv[i], cv[c_s], cv[c_e]) > tolerance)
		{
			stright_enough = false;
			break;
		}
	}

	//满足要求，不进一步细分了
	if (stright_enough)
	{
		//为了保证控制点不重复，设计的规则为[),但是对最后一个点例外。
		//按照递归顺序，最后一段首先加入points
		int c_end = points.empty() ? c_e + 1 : c_e;
		points.insert(points.begin(), cv.begin() + c_s, cv.begin() + c_end);
		return;
	}

	//从节点中间打断
	double u_mid = knots[k_s] + (knots[k_e] - knots[k_s]) / 2.0;
	//查找 u 所在区间
	int k = std::distance(knots.begin(), std::upper_bound(knots.begin(), knots.end(), u_mid)) - 1;

	std::vector<Point> cv_left, cv_right;
	_subdeviding(u_mid, k, knots, cv, cv_left, cv_right);

	//节点区间新增p个u_mid
	knots.insert(knots.begin() + k + 1, p, u_mid);
	//控制点替换
	cv.insert(cv.begin() + k, p, Point());
	for (int i = k - p, j = 0; j < p; ++j, ++i)
		cv[i] = cv_left[j];

	for (int i = k, j = 0; j <= p; ++j, ++i)
		cv[i] = cv_right[j];

	//两部分分别递归
	//Note:后半部分在前半部分之前执行，
	//因为如果前半部分首先执行的化，后半部分的索引就发生改变了
	_tesselation(cv, knots, points, tolerance, k + 1, k_e + p, k, c_e + p);
	_tesselation(cv, knots, points, tolerance, k_s, k + 1, c_s, k);
}

void BSpline::Subdeviding(double u, BSpline& sub_left, BSpline& sub_right)
{
	//为了使用方便
	const std::vector<double>& knots = m_vecKnots;
	const std::vector<Point>& cv = m_vecCVs;
	int p = m_nDegree;

	//查找 u 所在区间
	int k = std::distance(knots.begin(), std::upper_bound(knots.begin(), knots.end(), u)) - 1;

	//保存de boor 算法迭代过程中生成的控制点
	std::vector<Point> cv_left, cv_right;
	_subdeviding(u, k, knots, cv, cv_left, cv_right);

	//曲线在u处打断，生成新的控制点序列 假设原曲线的控制点序列为 P_0,...,P_k-p-1,P_k-p,...,P_k,P_k+1,...,P_n
	//打断后的新控制点序列为：左侧：P_0,....,P_k-p-1,cv_left[...]，右侧 cv_right[...],P_k+1,...,P_n。
	//sub_left
	sub_left.m_nDegree = p;
	//knots
	sub_left.m_vecKnots.resize(k + p + 2);
	for (int i = 0; i < k + 1; ++i)
		sub_left.m_vecKnots[i] = knots[i];
	for (int i = 0, j = k + 1; i < p + 1; ++i, ++j)
		sub_left.m_vecKnots[j] = u;
	//control vertex
	sub_left.m_vecCVs.resize(k + 1);
	for (int i = 0; i < k - p; ++i)
		sub_left.m_vecCVs[i] = cv[i];
	for (int i = 0, j = k - p; i < p + 1; ++i, ++j)
		sub_left.m_vecCVs[j] = cv_left[i];

	//sub_right
	sub_right.m_nDegree = p;
	//knots
	sub_right.m_vecKnots.resize(knots.size() - k + p);
	for (int i = 0; i < p + 1; i++)
		sub_right.m_vecKnots[i] = u;
	for (int i = p + 1, j = k + 1; j < knots.size(); ++j, ++i)
		sub_right.m_vecKnots[i] = knots[j];
	//control_vertex
	sub_right.m_vecCVs.resize(cv.size() - k + p);
	for (int i = 0; i < p + 1; ++i)
		sub_right.m_vecCVs[i] = cv_right[i];
	for (int i = k + 1, j = p + 1; i < cv.size(); ++i, ++j)
		sub_right.m_vecCVs[j] = cv[i];
}

void BSpline::_subdeviding(double u, int k, const std::vector<double>& knots, const std::vector<Point>& cv,
	std::vector<Point>& cv_left, std::vector<Point>& cv_right)
{
	int p = m_nDegree;

	//保存de boor 算法迭代过程中生成的控制点
	cv_left.resize(p + 1);
	cv_right.resize(p + 1);

	//将 P_k-p,...,P_k拷贝到cv_left上面
	std::copy(cv.begin() + k - p, cv.begin() + k + 1, cv_left.begin());
	cv_right[p] = cv_left[p];

	//de-boor 迭代p次
	for (int r = 1; r <= p; ++r)
	{
		//i 从 k 到 k-p+1
		for (int i = k, j = p; i >= k - p + r; --i, --j)
		{
			double alpha = u - knots[i];
			double dev = (knots[i + p + 1 - r] - knots[i]);
			alpha = (dev != 0) ? alpha / dev : 0;

			cv_left[j] = (1.0 - alpha)*cv_left[j - 1] + alpha * cv_left[j];
		}

		cv_right[p - r] = cv_left[p];
	}
}

double BSpline::PointLineDistance(const Point &s, const Point &a, const Point &b)
{
	double ab = sqrt(pow((a(0) - b(0)), 2.0) + pow((a(1) - b(1)), 2.0) + pow((a(2) - b(2)), 2.0));
	double as = sqrt(pow((a(0) - s(0)), 2.0) + pow((a(1) - s(1)), 2.0) + pow((a(2) - s(2)), 2.0));
	double bs = sqrt(pow((s(0) - b(0)), 2.0) + pow((s(1) - b(1)), 2.0) + pow((s(2) - b(2)), 2.0));
	double cos_A = (pow(as, 2.0) + pow(ab, 2.0) - pow(bs, 2.0)) / (2 * ab*as);
	double sin_A = sqrt(1 - pow(cos_A, 2.0));
	return as * sin_A;
}

/*!
 *\brief 三次B样条插值
*\ param const std::vector<Point> & vecFitPoints待插值点集合，需要点数不小于3
*\ Returns:   BSpline 插值样条曲线
*/
BSpline BSpline::CubicInterpolate(const std::vector<Point>& vecFitPoints)
{
	const int p = 3;
	BSpline bs;
	int x = vecFitPoints.size();
	if (x < p)
	{
		cout << "too less point !" << endl;
		return bs;
	}

	//求解方程 N*P = F
	Eigen::MatrixXd N = Eigen::MatrixXd::Zero(x + 2, x + 2);
	Eigen::MatrixXd P = Eigen::MatrixXd::Zero(x + 2, 3);
	Eigen::MatrixXd F = Eigen::MatrixXd::Zero(x + 2, 3);

	bs.m_nDegree = p;
	bs.m_vecKnots.resize(x); //x+6个节点
	//计算节点
	bs.m_vecKnots[0] = 0.0;
	for (int i = 1; i < x; ++i)
	{
		bs.m_vecKnots[i] = bs.m_vecKnots[i - 1]
			+ PointDistance(vecFitPoints[i], vecFitPoints[i - 1]);
	}
	//节点首尾构成p+1度重复
	bs.m_vecKnots.insert(bs.m_vecKnots.begin(), p, bs.m_vecKnots.front());
	bs.m_vecKnots.insert(bs.m_vecKnots.end(), p, bs.m_vecKnots.back());

	//1.填写矩阵N
	std::vector<double> basis_func;
	N(0, 0) = 1;
	N(x - 1, x + 1) = 1;
	for (int i = p + 1; i < x + p - 1; ++i)
	{
		//c(u)在 N_{i-p},...,N_i等p+1个基函数上非零
		bs.BasisFunc(bs.m_vecKnots[i], i, basis_func);
		for (int j = i - p, k = 0; j <= i; ++j, ++k)
		{
			N(i - p, j) = basis_func[k];
		}
	}

	//导数
	N(x, 0) = -1;
	N(x, 1) = 1;
	N(x + 1, x) = -1;
	N(x + 1, x + 1) = 1;

	//2.填写矩阵F
	for (int i = 0; i < x; ++i)
	{
		F(i, 0) = vecFitPoints[i](0);
		F(i, 1) = vecFitPoints[i](1);
		F(i, 2) = vecFitPoints[i](2);
	}

	{
		Vec3d v0, v1, v2;
		BesselTanget(vecFitPoints[0], vecFitPoints[1], vecFitPoints[2], v0, v1, v2);
		Vec3d v = v0 * (bs.m_vecKnots[p + 1] - bs.m_vecKnots[1]) / (double)p;
		F(x, 0) = v(0);
		F(x, 1) = v(1);
		F(x, 2) = v(2);
	}

	{
		Vec3d v0, v1, v2;
		BesselTanget(vecFitPoints[x - 3], vecFitPoints[x - 2], vecFitPoints[x - 1], v0, v1, v2);
		Vec3d v = v2 * (bs.m_vecKnots[x + 1 + p] - bs.m_vecKnots[x + 1]) / (double)p;
		F(x + 1, 0) = v(0);
		F(x + 1, 1) = v(1);
		F(x + 1, 2) = v(2);
	}

	//解方程 N*P = F
	P = N.lu().solve(F);

#ifdef _DEBUG
	cout << "N--------------" << endl << N << endl;
	cout << "F--------------" << endl << F << endl;
	cout << "P--------------" << endl << P << endl;
#endif

	//将Eigen所求的结果赋给bs的control_vertex
	bs.m_vecCVs.resize(x + 2);
	for (int i = 0; i < x + 2; ++i)
	{
		Point& cv = bs.m_vecCVs[i];
		cv(0) = P(i, 0);
		cv(1) = P(i, 1);
		cv(2) = P(i, 2);
	}

	return bs;
}