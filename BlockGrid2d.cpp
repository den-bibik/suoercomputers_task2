//
// Created by d.bibik on 15/12/2020.
//

#include <vector>
#include "BlockGrid2d.h"


BlockGrid2d::BlockGrid2d(point2d<int> block_start, point2d<int> block_end, point2d<int> gsize, point2d<double> border,
                         point2d<double> step) {
    this->gsize = gsize;
    this->border = border;
    this->step = step;
    this->bstart = block_start;
    this->bend = block_end;
    this->bsize = bend - bstart;

    int block_size = (bsize.values[0] + 2) * (bsize.values[1] + 2);
    this->values = std::__1::shared_ptr<std::__1::vector<double>>(new std::__1::vector<double>(block_size));
}

BlockGrid2d &BlockGrid2d::operator=(const BlockGrid2d &b) {
    this->gsize = b.gsize;
    this->border = b.border;
    this->step = b.step;
    this->bstart = b.bstart;
    this->bsize = b.bsize;

    copy_data(b, *this);

    return *this;
}

BlockGrid2d::BlockGrid2d(const BlockGrid2d &b) : BlockGrid2d(b.bstart, b.bend, b.gsize, b.border, b.step) {
    copy_data(b, *this);
}

point2d<double> BlockGrid2d::get_coord(point2d<int> point) const {
    return this->border + this->step * point;
}

double BlockGrid2d::get_value(point2d<int> p) const {
    bool inner[2];
    inner[0] = (p.values[0] >= bstart.values[0] - 1) && (p.values[0] <= bend.values[0] + 1);
    inner[1] = (p.values[1] >= bstart.values[1] - 1) && (p.values[1] <= bend.values[1] + 1);


    if (inner[0] && inner[1]) {
        return (*values)[(p.values[0] - (bstart.values[0] - 1)) * (bsize.values[1] + 2) +
                         (p.values[1] - (bstart.values[1] - 1))];
    } else {
        throw std::runtime_error("Trying to get not internal value");
    }
}

void BlockGrid2d::set_external_border_part(block_part key, std::__1::vector<double> value) {
    int i, j, start_val;
    switch (key) {
        case TOP:
            i = bstart.values[0] - 1;
            start_val = bstart.values[1];
            for (j = start_val; j < bend.values[1]; j++)
                this->set_value(point2d<int>(i, j), value[j - start_val], false);
            break;
        case BOTTOM:
            i = bend.values[0];
            start_val = bstart.values[1];
            for (j = start_val; j < bend.values[1]; j++)
                this->set_value(point2d<int>(i, j), value[j - start_val], false);
            break;
        case LEFT:
            j = bstart.values[1] - 1;
            start_val = bstart.values[0];
            for (i = start_val; i < bend.values[0]; i++)
                this->set_value(point2d<int>(i, j), value[i - start_val], false);
            break;
        case RIGHT:
            j = bend.values[1];
            start_val = bstart.values[0];
            for (i = start_val; i < bend.values[0]; i++)
                this->set_value(point2d<int>(i, j), value[i - start_val], false);
            break;

        case LT_CORNER:
            this->set_value(point2d<int>(bstart.values[0] - 1, bstart.values[1] - 1), value[0], false);
            break;
        case RT_CORNER:
            this->set_value(point2d<int>(bstart.values[0] - 1, bend.values[1]), value[0], false);
            break;
        case LB_CORNER:
            this->set_value(point2d<int>(bend.values[0], bstart.values[1] - 1), value[0], false);
            break;
        case RB_CORNER:
            this->set_value(point2d<int>(bend.values[0], bend.values[1]), value[0], false);
            break;
    }
}

std::__1::vector<double> BlockGrid2d::get_internal_border_part(block_part key) {
    std::__1::vector<double> result;
    int i, j, start_val;
    switch (key) {
        case TOP:
            i = bstart.values[0];
            start_val = bstart.values[1];
            for (j = start_val; j < bend.values[1]; j++)
                result.push_back(this->get_value(point2d<int>(i, j)));
            break;
        case BOTTOM:
            i = bend.values[0] - 1;
            start_val = bstart.values[1];
            for (j = start_val; j < bend.values[1]; j++)
                result.push_back(this->get_value(point2d<int>(i, j)));
            break;
        case LEFT:
            j = bstart.values[1];
            start_val = bstart.values[0];
            for (i = start_val; i < bend.values[0]; i++)
                result.push_back(this->get_value(point2d<int>(i, j)));
            break;
        case RIGHT:
            j = bend.values[1] - 1;
            start_val = bstart.values[0];
            for (i = start_val; i < bend.values[0]; i++)
                result.push_back(this->get_value(point2d<int>(i, j)));
            break;

        case LT_CORNER:
            result.push_back(this->get_value(point2d<int>(bstart.values[0], bstart.values[1])));
            break;
        case RT_CORNER:
            result.push_back(this->get_value(point2d<int>(bstart.values[0], bend.values[1] - 1)));
            break;
        case LB_CORNER:
            result.push_back(this->get_value(point2d<int>(bend.values[0] - 1, bstart.values[1])));
            break;
        case RB_CORNER:
            result.push_back(this->get_value(point2d<int>(bend.values[0] - 1, bend.values[1] - 1)));
            break;
    }
    return result;
}

void BlockGrid2d::set_value(point2d<int> p, double val, bool internal_op) {
    if (check_internal(p) || !internal_op) {
        int position =
                (p.values[0] - (bstart.values[0] - 1)) * (bsize.values[1] + 2) +
                (p.values[1] - (bstart.values[1] - 1));
        (*values)[position] = val;
    } else {
        throw std::runtime_error("Trying to set not internal value");
    }
}

BlockGrid2d operator-(const BlockGrid2d &a, const BlockGrid2d &b) {
    auto res = BlockGrid2d(a.bstart, a.bend, a.gsize, a.border, a.step);
    for (int i = res.bstart.values[0] - 1; i <= res.bend.values[0]; i++) {
        for (int j = res.bstart.values[1] - 1; j <= res.bend.values[1]; j++) {
            auto p = point2d<int>(i, j);
            res.set_value(p, a.get_value(p) - b.get_value(p), false);
        }
    }
    return res;
}

BlockGrid2d operator*(const BlockGrid2d &a, const double C) {
    auto res = BlockGrid2d(a.bstart, a.bend, a.gsize, a.border, a.step);
    for (int i = res.bstart.values[0] - 1; i <= res.bend.values[0]; i++) {
        for (int j = res.bstart.values[1] - 1; j <= res.bend.values[1]; j++) {
            auto p = point2d<int>(i, j);
            res.set_value(p, a.get_value(p) * C, false);
        }
    }
    return res;
}

BlockGrid2d operator*(const double C, const BlockGrid2d &a) {
    return operator*(a, C);
}

double BlockGrid2d::dot_prod(const BlockGrid2d &b) const {
    double res = 0;
    for (int i = bstart.values[0]; i < bend.values[0]; i++) {
        for (int j = bstart.values[1]; j < bend.values[1]; j++) {
            auto p = point2d<int>(i, j);
            double val = get_value(p) * b.get_value(p);
            if ((i == 0) || (i == gsize.values[0])) val /= 2;
            if ((j == 0) || (j == gsize.values[1])) val /= 2;
            res += val;
        }
    }
    res *= step.values[0] * step.values[1];
    return res;
}

double BlockGrid2d::grad(point2d<int> p, int axis, bool forward_shift) const {
    auto b = point2d<int>(p.values[0], p.values[1]);
    if (forward_shift) b.values[axis] += 1; else b.values[axis] -= 1;
    double result = (get_value(p) - get_value(b)) / step.values[axis];
    if (forward_shift) result *= -1;

    return result;
}

double BlockGrid2d::aw_grad(point2d<int> p, const BlockGrid2d &a, int axis) const {
    auto p2 = point2d<int>(p.values[0], p.values[1]);
    p2.values[axis] += 1;
    double res = (
                         a.get_value(p2) * this->grad(p2, axis, false) -
                         a.get_value(p) * this->grad(p, axis, false)
                 ) / step.values[axis];
    return res;
}

double BlockGrid2d::laplass(point2d<int> p, const BlockGrid2d &a, const BlockGrid2d &b) const {
    double res = aw_grad(p, a, 0) + aw_grad(p, b, 1);
    return res;
}

void BlockGrid2d::copy_data(const BlockGrid2d &from, BlockGrid2d &to) {
    to.values = from.values;
    //copy(from.values->begin(), from.values->end(), to.values->begin());
}

bool BlockGrid2d::check_internal(const point2d<int> point) const {
    return (point.values[0] >= bstart.values[0]) && (point.values[0] < bend.values[0]) &&
           (point.values[1] >= bstart.values[1]) && (point.values[1] < bend.values[1]);
}
