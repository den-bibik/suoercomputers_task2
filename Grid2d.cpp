//
// Created by d.bibik on 10/12/2020.
//

#include <iostream>
#include <vector>
#include "Grid2d.h"

Grid2d::Grid2d(point2d<int> _gsize, point2d<double> _border, point2d<double> _step) {
    gsize = _gsize;
    border = _border;
    step = _step;

    int values_size = (gsize.values[0] + 1) * (gsize.values[1] + 1);
    values = new double[values_size];
}

Grid2d::Grid2d(const Grid2d &b) {
    gsize = b.gsize;
    border = b.border;
    step = b.step;

    int values_size = (gsize.values[0] + 1) * (gsize.values[1] + 1);
    values = new double[values_size];
    memcpy(values, b.values, sizeof(*b.values) * values_size);
}

Grid2d::~Grid2d() {
    delete[] values;
}

Grid2d &Grid2d::operator=(const Grid2d &b) {
    gsize = b.gsize;
    border = b.border;
    step = b.step;

    int values_size = (gsize.values[0] + 1) * (gsize.values[1] + 1);
    //values = new double[values_size];
    memcpy(values, b.values, sizeof(*b.values) * values_size);

    return *this;
}

point2d<double> Grid2d::get_coord(point2d<int> point) const {
    return border + step * point;
}

double Grid2d::get_value(point2d<int> p) const {
    if ((p.values[0] < 0) || (p.values[1] < 0)) { p.values[0] = 1 / 0; }
    return values[p.values[0] * (gsize.values[1] + 1) + p.values[1]];
}

void Grid2d::set_value(point2d<int> position, double val) {
    if ((position.values[0] < 0) || (position.values[1] < 0)) { position.values[0] = 1 / 0; }
    values[position.values[0] * (gsize.values[1] + 1) + position.values[1]] = val;
}

void Grid2d::dump() const {
    std::cout << "matrix dump" << std::endl;
    for (int j = gsize.values[1]; j >= 0; j--) {
        for (int i = 0; i <= gsize.values[0]; i++) {
            auto p = point2d<int>(i, j);
            std::cout << get_value(p) << " ";
        }
        std::cout << std::endl;
    }
}

Grid2d Grid2d::operator-(const Grid2d &b) const {
    Grid2d res = Grid2d(gsize, border, step);
    for (int i = 0; i <= gsize.values[0]; i++) {
        for (int j = 0; j <= gsize.values[1]; j++) {
            auto p = point2d<int>(i, j);
            res.set_value(p, get_value(p) - b.get_value(p));
        }
    }
    return res;
}

Grid2d Grid2d::operator*(const double C) const {
    Grid2d res = Grid2d(gsize, border, step);
    for (int i = 0; i <= gsize.values[0]; i++) {
        for (int j = 0; j <= gsize.values[1]; j++) {
            auto p = point2d<int>(i, j);
            res.set_value(p, get_value(p) * C);
        }
    }
    return res;
}

Grid2d Grid2d::operator*(const Grid2d &b) const {
    Grid2d res = Grid2d(gsize, border, step);
    for (int i = 0; i <= gsize.values[0]; i++) {
        for (int j = 0; j <= gsize.values[1]; j++) {
            auto p = point2d<int>(i, j);
            res.set_value(p, get_value(p) * b.get_value(p));
        }
    }
    return res;
}

void Grid2d::print_info() const {
    double max = -1000000000;
    double min = 10000000000;
    double sum = 0.0;
    for (int i = 0; i <= gsize.values[0]; i++) {
        for (int j = 0; j <= gsize.values[1]; j++) {
            auto p = point2d<int>(i, j);
            double val = get_value(p);
            if (val > max) max = val;
            if (val < min) min = val;
            sum += val;

        }
    }
    sum /= (gsize.values[0] + 1) * (gsize.values[1] + 1);

    std::cout << " mean " << sum << " max " << max << " min " << min << std::endl;
}

double Grid2d::dot_prod(const Grid2d &b) const {
    int values_size = gsize.values[0] * gsize.values[1];
    double res = 0;
    for (int i = 0; i <= gsize.values[0]; i++) {
        for (int j = 0; j <= gsize.values[1]; j++) {
            auto p = point2d<int>(i, j);
            double val = get_value(p) * b.get_value(p);
            if ((i == 0) | (i == gsize.values[0])) val /= 2;
            if ((j == 0) | (j == gsize.values[1])) val /= 2;
            res += val;
        }
    }
    res *= step.values[0] * step.values[1];
    return res;
}

double Grid2d::sum() const {
    double res = 0;
    for (int i = 0; i <= gsize.values[0]; i++) {
        for (int j = 0; j <= gsize.values[1]; j++) {
            auto p = point2d<int>(i, j);
            double val = get_value(p);
            res += val;
        }
    }
    return res;
}

double Grid2d::grad(point2d<int> p, int axis, bool forward_shift) const {
    auto b = point2d<int>(p.values[0], p.values[1]);
    if (forward_shift) b.values[axis] += 1; else b.values[axis] -= 1;
    double result = (get_value(p) - get_value(b)) / step.values[axis];
    if (forward_shift) result *= -1;

    return result;
}

double Grid2d::aw_grad(point2d<int> p, const Grid2d *a, int axis) const {
    auto p2 = point2d<int>(p.values[0], p.values[1]);
    p2.values[axis] += 1;
    double res = (
                         a->get_value(p2) * grad(p2, axis, 0) -
                         a->get_value(p) * grad(p, axis, 0)
                 ) / step.values[axis];
    return res;
}

double Grid2d::laplass(point2d<int> p, const Grid2d *a, const Grid2d *b) const {
    double res = aw_grad(p, a, 0) + aw_grad(p, b, 1);
    return res;
}
