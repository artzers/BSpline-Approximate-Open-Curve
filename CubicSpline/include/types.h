/*
Part of B-Splines using Eigen - a class for fitting and sampling B-splines from data.
Copyright (C) 2021  Jack Naylor

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. 
*/

#pragma once


#include "Eigen/Core"
#include "Eigen/Geometry"


#define _USE_MATH_DEFINES
#include <cmath>

using namespace Eigen;

using HPoint = Vector4d;
using Point = Vector3d;
using Vec3d = Vector3d;
using Pixel = Vector2d;
using Quat = Quaterniond;
using HomogenousMatrix = Matrix4d;
using RotationMatrix = Matrix3d;
using ProjectionMatrix = Matrix<double,3,4>;

using PixelArray = Matrix<double, 2, Dynamic>;
using PointArray = Matrix<double, 3, Dynamic>;
using HPointArray = Matrix<double, 4, Dynamic>;


