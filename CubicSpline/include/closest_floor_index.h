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

#include "types.h"
#include <algorithm>

/**
 * @brief Return closest index of a value in an nx1 Eigen vector
 * 
 * @param vector Vector to compare to
 * @param candi Candidate value
 * @return int Index
 */
int closest_floor_index(const VectorXd& vector, const double& candi) {
    VectorXd sorted_vector = vector;
    std::sort(sorted_vector.begin(),sorted_vector.end());
    VectorXd diff_vec = sorted_vector.array()-candi;

    for (int i = 0; i< diff_vec.size(); i++) {
        if (diff_vec(i)>0 && i!=0) {
            return i-1;
        } else if (diff_vec(i)>0 && i==0) {
            return i;
        }
    }

    return diff_vec.size();


}