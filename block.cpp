//
// Created by d.bibik on 30/11/2020.
//
#include <iostream>
#include <cmath>
#include <vector>
#include "Grid2d.h"

using namespace std;

enum border_part {
    TOP, BOTTOM, LEFT, RIGHT, LT_CORNER, RT_CORNER, LB_CORNER, RB_CORNER
};

class BlockGrid2d {
public:
    BlockGrid2d() = default;

    BlockGrid2d(point2d<int> block_start, point2d<int> block_end, point2d<int> gsize,
                point2d<double> border, point2d<double> step) {
        this->gsize = gsize;
        this->border = border;
        this->step = step;
        this->bstart = block_start;
        this->bend = block_end;
        this->bsize = bend - bstart;

        int block_size = (bsize.values[0] + 1) * (bsize.values[1] + 1);
        this->values = vector<double>(block_size);

        this->left_border = vector<double>(bsize.values[1] + 1);
        this->right_border = vector<double>(bsize.values[1] + 1);
        this->top_border = vector<double>(bsize.values[0] + 1);
        this->bottom_border = vector<double>(bsize.values[0] + 1);
    }

    BlockGrid2d &operator=(const BlockGrid2d &b) {
        this->gsize = b.gsize;
        this->border = b.border;
        this->step = b.step;
        this->bstart = b.bstart;
        this->bsize = b.bsize;

        copy_data(b, *this);

        return *this;
    }

    BlockGrid2d(const BlockGrid2d &b) : BlockGrid2d(b.bstart, b.bend, b.gsize, b.border, b.step) {
        copy_data(b, *this);
    }

    point2d<double> get_coord(point2d<int> point) const {
        return this->border + this->step * point;
    }

    double get_value(point2d<int> p) const {
        bool inner[2], valid[2];
        inner[0] = (p.values[0] >= bstart.values[0]) && (p.values[0] <= bend.values[0]);
        inner[1] = (p.values[1] >= bstart.values[1]) && (p.values[1] <= bend.values[1]);


        valid[0] = (p.values[0] >= std::max(bstart.values[0] - 1, 0)) &&
                   (p.values[0] <= std::min(bend.values[0] + 1, gsize.values[0]));

        valid[1] = (p.values[1] >= std::max(bstart.values[1] - 1, 0)) &&
                   (p.values[1] <= std::min(bend.values[1] + 1, gsize.values[0]));


        if (inner[0] && inner[1])
            return values[(p.values[0] - bstart.values[0]) * (gsize.values[1] + 1) +
                          (p.values[1] - bstart.values[1])];
        else if ((p.values[0] == bstart.values[0] - 1) && valid[1])
            return bottom_border[p.values[1] - bstart.values[1] + 1];
        else if ((p.values[0] == bend.values[0] + 1) && valid[1])
            return top_border[p.values[1] - bstart.values[1] + 1];
        else if ((p.values[1] == bstart.values[1] - 1) && valid[0])
            return left_border[p.values[0] - bstart.values[0] + 1];
        else if ((p.values[1] == bend.values[1] + 1) && valid[0])
            return right_border[p.values[0] - bstart.values[0] + 1];

        return NAN;
    }


    void set_border_part(border_part key, vector<double> value) {
        int i, j, start_val;
        switch (key) {
            case TOP:
                i = bstart.values[0];
                start_val = bstart.values[1];
                for (j = start_val; j <= bend.values[1]; j++)
                    this->set_value(point2d<int>(i, j), value[j - start_val]);
                break;
            case BOTTOM:
                i = bend.values[0];
                start_val = bstart.values[1];
                for (j = start_val; j <= bend.values[1]; j++)
                    this->set_value(point2d<int>(i, j), value[j - start_val]);
                break;
            case LEFT:
                j = bstart.values[1];
                start_val = bstart.values[0];
                for (i = start_val; i <= bend.values[0]; i++)
                    this->set_value(point2d<int>(i, j), value[j - start_val]);
                break;
            case RIGHT:
                j = bend.values[1];
                start_val = bstart.values[0];
                for (i = start_val; i <= bend.values[0]; i++)
                    this->set_value(point2d<int>(i, j), value[j - start_val]);
                break;


            case LT_CORNER:
                break;
            case RT_CORNER:
                break;
            case LB_CORNER:
                break;
            case RB_CORNER:
                break;
        }
    }

    void set_value(point2d<int> p, double val) {
        values[(p.values[0] - bstart.values[0]) * (bsize.values[1] + 1) + (p.values[1] - bstart.values[1])] = val;
    }

    friend BlockGrid2d operator-(const BlockGrid2d &a, const BlockGrid2d &b) {
        auto res = BlockGrid2d(a.bstart, a.bend, a.bsize, a.border, a.step);
        for (int i = res.bstart.values[0]; i <= res.bend.values[0]; i++) {
            for (int j = res.bstart.values[1]; j <= res.bend.values[1]; j++) {
                auto p = point2d<int>(i, j);
                res.set_value(p, a.get_value(p) - b.get_value(p));
            }
        }
        return res;
    }

    friend BlockGrid2d operator*(const BlockGrid2d &a, const double C) {
        auto res = BlockGrid2d(a.bstart, a.bend, a.bsize, a.border, a.step);
        for (int i = res.bstart.values[0]; i <= res.bend.values[0]; i++) {
            for (int j = res.bstart.values[1]; j <= res.bend.values[1]; j++) {
                auto p = point2d<int>(i, j);
                res.set_value(p, a.get_value(p) * C);
            }
        }
        return res;
    }

    friend BlockGrid2d operator*(const double C, const BlockGrid2d &a) {
        return operator*(a, C);
    }

    double dot_prod(const BlockGrid2d &b) const {
        double res = 0;
        for (int i = bstart.values[0]; i <= bend.values[0]; i++) {
            for (int j = bstart.values[1]; j <= bend.values[1]; j++) {
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

    double grad(point2d<int> p, int axis, bool forward_shift) const {
        auto b = point2d<int>(p.values[0], p.values[1]);
        if (forward_shift) b.values[axis] += 1; else b.values[axis] -= 1;
        double result = (get_value(p) - get_value(b)) / step.values[axis];
        if (forward_shift) result *= -1;

        return result;
    }

    double aw_grad(point2d<int> p, const BlockGrid2d &a, int axis) const {
        auto p2 = point2d<int>(p.values[0], p.values[1]);
        p2.values[axis] += 1;
        double res = (
                             a.get_value(p2) * this->grad(p2, axis, 0) -
                             a.get_value(p) * this->grad(p, axis, 0)
                     ) / step.values[axis];
        return res;
    }

    double laplass(point2d<int> p, const BlockGrid2d &a, const BlockGrid2d &b) const {
        double res = aw_grad(p, a, 0) + aw_grad(p, b, 1);
        return res;
    }

    point2d<int> bstart, bsize, bend, gsize;
    point2d<double> border, step;

private:
    vector<double> values;
    vector<double> left_border, right_border, top_border, bottom_border;
    double lt_corner, lb_corner,  rt_corner, rb_corner;


    static void copy_data(const BlockGrid2d &from, BlockGrid2d &to) {
        copy(from.values.begin(), from.values.end(), to.values.begin());
        copy(from.left_border.begin(), from.left_border.end(), to.left_border.begin());
        copy(from.right_border.begin(), from.right_border.end(), to.right_border.begin());
        copy(from.top_border.begin(), from.top_border.end(), to.top_border.begin());
        copy(from.bottom_border.begin(), from.bottom_border.end(), to.bottom_border.begin());
    }

};

class BlockContainer {
public:
    BlockContainer(point2d<int> gsize, point2d<double> gstart, point2d<double> gstep,
                   point2d<int> num_blocks) : num_blocks(num_blocks) {
        point2d<int> block_size = gsize / num_blocks;
        int xmod = gsize.values[0] % num_blocks.values[0];

        auto block_start = point2d<int>(0, 0);
        for (int i = 0; i < num_blocks.values[0]; i++) {
            int xsize = block_size.values[0];
            if (xmod > 0) {
                xmod--;
                xsize++;
            }
            block_start.values[1] = 0;

            int ymod = gsize.values[1] % num_blocks.values[1];
            for (int j = 0; j < num_blocks.values[0]; j++) {
                int ysize = block_size.values[0];
                if (ymod > 0) {
                    ymod--;
                    ysize++;
                }

                auto block = BlockGrid2d(block_start, block_start + block_size, gsize, gstart, gstep);
                matrix.push_back(block);

                block_start.values[1] += ysize;
            }
            block_start.values[0] += xsize;
        }

    }


    void set(int i, int j, const BlockGrid2d &val) {
        matrix[i * num_blocks.values[1] + j] = val;
    }

    BlockGrid2d &get(int i, int j) {
        return matrix[i * num_blocks.values[1] + j];
    }

    void sync_borders() {
        for (int i = 0; i <= num_blocks.values[0]; i++) {
            for (int j = 0; j <= num_blocks.values[0]; j++) {
                send_border_part(i, j, i - 1, j, border_part::TOP, border_part::BOTTOM);
                send_border_part(i, j, i + 1, j, border_part::BOTTOM, border_part::TOP);
                send_border_part(i, j, i, j - 1, border_part::LEFT, border_part::RIGHT);
                send_border_part(i, j, i, j + 1, border_part::RIGHT, border_part::LEFT);

                send_border_part(i, j, i - 1, j - 1, border_part::LT_CORNER, border_part::RB_CORNER);
                send_border_part(i, j, i + 1, j - 1, border_part::LB_CORNER, border_part::RT_CORNER);
                send_border_part(i, j, i - 1, j + 1, border_part::RT_CORNER, border_part::LB_CORNER);
                send_border_part(i, j, i + 1, j + 1, border_part::RB_CORNER, border_part::LT_CORNER);
            }
        }
    }


    point2d<int> num_blocks;
private:
    void send_border_part(int from_i, int from_j, int to_i, int to_j, border_part from_key, border_part to_key) {
        auto from = this->get(from_i, from_i);
        if ((to_i >= 0) && (to_i <= num_blocks.values[0]) &&
            (to_j >= 0) && (to_j <= num_blocks.values[1])) {
            auto to = this->get(to_i, to_j);
            to.set_border_part(to_key, from.get_border_part(from_key));
            this->set(to_i, to_j, to);
        }
    }

    vector<BlockGrid2d> matrix;
};

double F_func(point2d<double> point) {
    // calculated from F(x,y) = -laplass(u) + q(x, y) * u,  u = 1 + cos(pi * x * y), k = 4 + x + y, z = 0
    double x = point.values[0];
    double y = point.values[1];
    return M_PI * (M_PI * (x + y + 4) * (x * x + y * y) * cos(M_PI * x * y) + (x + y) * sin(M_PI * x * y));
}

double psi_R_func(point2d<double> point) {
    double x = point.values[0];
    double y = point.values[1];
    return -M_PI * y * (x + y + 4) * sin(M_PI * x * y) + cos(M_PI * x * y) + 1;
}

double psi_L_func(point2d<double> point) {
    double x = point.values[0];
    double y = point.values[1];
    return -M_PI * y * (x + y + 4) * sin(M_PI * x * y);
}

double psi_T_func(point2d<double> point) {
    double x = point.values[0];
    double y = point.values[1];
    return -M_PI * x * (x + y + 4) * sin(M_PI * x * y) + cos(M_PI * x * y) + 1;
}

double psi_B_func(point2d<double> point) {
    double x = point.values[0];
    double y = point.values[1];
    return -M_PI * x * (x + y + 4) * sin(M_PI * x * y);
}

inline double k_func(point2d<double> point) {
    double x = point.values[0];
    double y = point.values[1];
    return 4 + x + y;
}

inline double u_func(point2d<double> point) {
    double x = point.values[0];
    double y = point.values[1];
    return 1 + cos(M_PI * x * y);
}


double calc_Aw_val(const point2d<int> &position, const BlockGrid2d &w, const BlockGrid2d &a, const BlockGrid2d &b) {
    int N = w.gsize.values[0];
    int M = w.gsize.values[1];
    int i = position.values[0];
    int j = position.values[1];
    double xstep = w.step.values[0];
    double ystep = w.step.values[1];

    if ((i >= 1) && (i <= N - 1) && (j >= 1) && (j <= M - 1)) { // inner_point
        return -w.laplass(position, a, b);
    }

        //  Borders
    else if ((i >= 1) && (i <= N - 1) && (j == M)) { // top
        return 2 / ystep * (
                b.get_value(position) * w.grad(position, 1, 0)
                + w.get_value(position)
        ) - w.aw_grad(position, a, 0);
    } else if ((i >= 1) && (i <= N - 1) && (j == 0)) { // bottom (eq. 9)
        auto position_i1 = point2d<int>(i, 1);
        return 2 / ystep * (
                -b.get_value(position_i1) * w.grad(position_i1, 1, 0)
                + 0 * w.get_value(position)
        ) - w.aw_grad(position, a, 0);
    } else if ((i == 0) && (j >= 1) && (j <= M - 1)) { // left (eq. 8)
        auto position_1j = point2d<int>(1, j);
        return 2 / xstep * (
                -a.get_value(position_1j) * w.grad(position_1j, 0, 0)
                + 0 * w.get_value(position)
        ) - w.aw_grad(position, b, 1);
    } else if ((i == N) && (j >= 1) && (j <= M - 1)) { // right (eq. 8)
        return 2 / xstep * (
                a.get_value(position) * w.grad(position, 0, 0)
                + 1 * w.get_value(position)
        ) - w.aw_grad(position, b, 1);
    }

        //Corners (eq 10-13)
    else if ((i == 0) && (j == 0)) { // LB corner
        auto p00 = point2d<int>(0, 0);
        auto p10 = point2d<int>(1, 0);
        auto p01 = point2d<int>(0, 1);
        return
                -2 / xstep * a.get_value(p10) * w.grad(p10, 0, 0) +
                -2 / ystep * b.get_value(p01) * w.grad(p01, 1, 0) +
                (0 * 2 / xstep + 0 * 2 / ystep) * w.get_value(p00);
    } else if ((i == N) && (j == 0)) { // RB corner
        auto pM0 = point2d<int>(N, 0);
        auto pM1 = point2d<int>(N, 1);
        return 2 / xstep * a.get_value(pM0) * w.grad(pM0, 0, 0) +
               -2 / ystep * b.get_value(pM1) * w.grad(pM1, 1, 0) +
               (1 * 2 / xstep + 0 * 2 / ystep) * w.get_value(pM0) * 0;
    } else if ((i == N) && (j == M)) { // RT corner
        return 2 / xstep * a.get_value(position) * w.grad(position, 0, 0) +
               2 / ystep * b.get_value(position) * w.grad(position, 1, 0) +
               (1 * 2 / xstep + 1 * 2 / ystep) * w.get_value(position) * 0;
    } else if ((i == 0) && (j == M)) { // LT corner
        auto p0N = point2d<int>(0, M);
        auto p1N = point2d<int>(1, M);
        return -2 / xstep * a.get_value(p1N) * w.grad(p1N, 0, 0) +
               2 / ystep * b.get_value(p0N) * w.grad(p0N, 1, 0) +
               (0 * 2 / xstep + 1 * 2 / ystep) * w.get_value(p0N) * 0;

    }

    return 0;
}

BlockGrid2d &Aw(BlockGrid2d &Aw, const BlockGrid2d &w, const BlockGrid2d &a, const BlockGrid2d &b) {
    for (int i = Aw.bstart.values[0]; i <= Aw.bend.values[0]; i++) {
        for (int j = Aw.bstart.values[1]; j <= Aw.bend.values[1]; j++) {
            auto position = point2d<int>(i, j);
            double value = calc_Aw_val(position, w, a, b);
            Aw.set_value(position, value);
        }
    }
    return Aw;
}


double calc_B_val(const point2d<int> &position, const point2d<double> &coord,
                  const point2d<double> &step, const point2d<int> &gsize) {
    int N = gsize.values[0];
    int M = gsize.values[1];
    int i = position.values[0];
    int j = position.values[1];
    double xstep = step.values[0];
    double ystep = step.values[1];

    if ((i >= 1) && (i <= N - 1) && (j >= 1) && (j <= M - 1)) { // inner_point
        return F_func(coord);
    }
        //  Borders
    else if ((i >= 1) && (i <= N - 1) && (j == M)) { // top
        return F_func(coord) + 2 / ystep * psi_T_func(coord);
    } else if ((i >= 1) && (i <= N - 1) && (j == 0)) { // bottom (eq. 9)
        return F_func(coord) + 2 / ystep * psi_B_func(coord);
    } else if ((i == 0) && (j >= 1) && (j <= M - 1)) { // left (eq. 8)
        return F_func(coord) + 2 / xstep * psi_L_func(coord);
    } else if ((i == N) && (j >= 1) && (j <= M - 1)) { // right (eq. 8)
        return F_func(coord) + 2 / xstep * psi_R_func(coord);
    }
        //Corners (eq 10-13)
    else if ((i == 0) && (j == 0)) { // LB corner
        return F_func(coord);
    } else if ((i == N) && (j == 0)) { // RB corner
        double coef_x = 2 / xstep;
        double coef_y = 2 / ystep;
        return F_func(coord) + coef_x * psi_R_func(coord) + coef_y * psi_B_func(coord);
    } else if ((i == N) && (j == M)) { // RT corner
        double coef_x = 2 / xstep;
        double coef_y = 2 / ystep;
        return F_func(coord) + coef_x * psi_R_func(coord) + coef_y * psi_T_func(coord);
    } else if ((i == 0) && (j == M)) { // LT corner
        double coef_x = 2 / xstep;
        double coef_y = 2 / ystep;
        return F_func(coord) + coef_x * psi_L_func(coord) + coef_y * psi_T_func(coord);
    }

    return 0;
}


BlockGrid2d &init_B(BlockGrid2d &block) {
    for (int i = block.bstart.values[0]; i <= block.bend.values[0]; i++) {
        for (int j = block.bstart.values[1]; j <= block.bend.values[1]; j++) {
            auto position = point2d<int>(i, j);
            auto coord = block.get_coord(position);
            double value = calc_B_val(position, coord, block.step, block.gsize);
            block.set_value(position, value);
        }
    }
    return block;
}


BlockGrid2d &init_k(BlockGrid2d &block, int axis) {
    for (int i = block.bstart.values[0]; i <= block.bend.values[0]; i++) {
        for (int j = block.bstart.values[1]; j <= block.bend.values[1]; j++) {
            auto position = point2d<int>(i, j);
            auto coord = block.get_coord(position);
            coord.values[axis] -= 0.5 * block.step.values[axis];
            double value = k_func(coord);
            block.set_value(position, value);
        }
    }
    return block;
}


BlockGrid2d &init_w_test(BlockGrid2d &block) {
    for (int i = block.bstart.values[0]; i <= block.bend.values[0]; i++) {
        for (int j = block.bstart.values[1]; j <= block.bend.values[1]; j++) {
            auto position = point2d<int>(i, j);
            auto coord = block.get_coord(position);
            block.set_value(position, u_func(coord));
        }
    }
    return block;
}

BlockGrid2d &init_w(BlockGrid2d &block) {
    for (int i = block.bstart.values[0]; i <= block.bend.values[0]; i++) {
        for (int j = block.bstart.values[1]; j <= block.bend.values[1]; j++) {
            auto position = point2d<int>(i, j);
            block.set_value(position, 0);
        }
    }
    return block;
}


class AlgoVariables {
public:
    BlockContainer B, w, w_test, a, b, Aw_val, r, Ar_val;

    AlgoVariables(point2d<int> gsize, point2d<double> gstart, point2d<double> gstep,
                  point2d<int> num_blocks) :
            B(gsize, gstart, gstep, num_blocks),
            w(gsize, gstart, gstep, num_blocks),
            w_test(gsize, gstart, gstep, num_blocks),
            a(gsize, gstart, gstep, num_blocks),
            b(gsize, gstart, gstep, num_blocks),
            //don't init
            Aw_val(gsize, gstart, gstep, num_blocks),
            r(gsize, gstart, gstep, num_blocks),
            Ar_val(gsize, gstart, gstep, num_blocks) {
        for (int i = 0; i < B.num_blocks.values[0]; i++) {
            for (int j = 0; j < B.num_blocks.values[1]; j++) {
                auto B_block = B.get(i, j);
                auto w_block = w.get(i, j);
                auto w_test_block = w_test.get(i, j);
                auto a_block = a.get(i, j);
                auto b_block = b.get(i, j);

                B_block = init_B(B_block);
                w_block = init_w(w_block);
                w_test_block = init_w_test(w_test_block);
                a_block = init_k(a_block, 0);
                b_block = init_k(b_block, 1);

                B.set(i, j, B_block);
                w.set(i, j, w_block);
                w_test.set(i, j, w_test_block);
                a.set(i, j, a_block);
                b.set(i, j, b_block);
            }
        }
    }
};

void block_algo(int size, int x_blocks, int y_blocks, int num_iter) {
    auto gsize = point2d<int>(size, size);
    auto border = point2d<double>(0, 0);
    auto border2 = point2d<double>(2, 1);
    auto step = (border2 - border) / gsize;
    auto num_blocks = point2d<int>(x_blocks, y_blocks);

    AlgoVariables vars(gsize, border, step, num_blocks);

    for (int iter = 0; iter < num_iter; iter++) {
        double alpha_numerator = 0;
        double alpha_denominator = 0;
        double r_norm = 0;
        for (int i = 0; i < vars.B.num_blocks.values[0]; i++) {
            for (int j = 0; j < vars.B.num_blocks.values[0]; j++) {
                auto B_block = vars.B.get(i, j);
                auto a_block = vars.a.get(i, j);
                auto b_block = vars.b.get(i, j);
                auto Aw_val_block = vars.Aw_val.get(i, j);
                auto w_block = vars.w.get(i, j);
                auto Ar_val_block = vars.Ar_val.get(i, j);
                auto r_block = vars.r.get(i, j);

                Aw_val_block = Aw(Aw_val_block, w_block, a_block, b_block);
                r_block = Aw_val_block - B_block;
                Ar_val_block = Aw(Ar_val_block, r_block, a_block, b_block);

                vars.r.set(i, j, r_block);

                //sync: (sum) alpha_numerator, alpha_denominator, r_norm
                alpha_numerator += Ar_val_block.dot_prod(r_block);
                alpha_denominator += Ar_val_block.dot_prod(Ar_val_block);
                r_norm += r_block.dot_prod(r_block);
            }
        }
        double alpha = alpha_numerator / alpha_denominator;
        //TODO: sync (update) r_block borders
        vars.r.sync_borders();


        double stop_norm = 0;
        double diff_val = 0;
        for (int i = 0; i < vars.B.num_blocks.values[0]; i++) {
            for (int j = 0; j < vars.B.num_blocks.values[0]; j++) {
                //need r_block, w_block, w_test_block
                auto w_block = vars.w.get(i, j);
                auto w_test_block = vars.w_test.get(i, j);
                auto r_block = vars.r.get(i, j);


                stop_norm += alpha * alpha * r_norm;
                auto alg_step = r_block * alpha;

                auto diff = w_block - w_test_block;
                diff_val += diff.dot_prod(diff);

                w_block = w_block - alg_step;

                vars.w.set(i, j, w_block);

            }
        }
        std::cout << "iter: " << iter << " stop criterion val: " << stop_norm <<
                  " test: " << diff_val << std::endl;

    }

}

int main() {
    std::cout << "It's block version" << std::endl;
    block_algo(2, 2, 2, 10);
    //algo(20, 100, 1, 0);
    return 0;
}