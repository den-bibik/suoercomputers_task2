//
// Created by d.bibik on 15/12/2020.
//

#include "AlgoVariables.h"

AlgoVariables::AlgoVariables(point2d<int> gsize, point2d<double> gstart, point2d<double> gstep, point2d<int> num_blocks)
        :
        B(gsize, gstart, gstep, num_blocks),
        w(gsize, gstart, gstep, num_blocks),
        w_test(gsize, gstart, gstep, num_blocks),
        a(gsize, gstart, gstep, num_blocks),
        b(gsize, gstart, gstep, num_blocks),
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
    B.sync_borders();
    w.sync_borders();
    w_test.sync_borders();
    a.sync_borders();
    b.sync_borders();
}

double AlgoVariables::F_func(point2d<double> point) {
    // calculated from F(x,y) = -laplass(u) + q(x, y) * u,  u = 1 + cos(pi * x * y), k = 4 + x + y, z = 0
    double x = point.values[0];
    double y = point.values[1];
    return M_PI * (M_PI * (x + y + 4) * (x * x + y * y) * cos(M_PI * x * y) + (x + y) * sin(M_PI * x * y));
}

double AlgoVariables::psi_R_func(point2d<double> point) {
    double x = point.values[0];
    double y = point.values[1];
    return -M_PI * y * (x + y + 4) * sin(M_PI * x * y) + cos(M_PI * x * y) + 1;
}

double AlgoVariables::psi_L_func(point2d<double> point) {
    double x = point.values[0];
    double y = point.values[1];
    return -M_PI * y * (x + y + 4) * sin(M_PI * x * y);
}

double AlgoVariables::psi_T_func(point2d<double> point) {
    double x = point.values[0];
    double y = point.values[1];
    return -M_PI * x * (x + y + 4) * sin(M_PI * x * y) + cos(M_PI * x * y) + 1;
}

double AlgoVariables::psi_B_func(point2d<double> point) {
    double x = point.values[0];
    double y = point.values[1];
    return -M_PI * x * (x + y + 4) * sin(M_PI * x * y);
}

double AlgoVariables::k_func(point2d<double> point) {
    double x = point.values[0];
    double y = point.values[1];
    return 4 + x + y;
}

double AlgoVariables::u_func(point2d<double> point) {
    double x = point.values[0];
    double y = point.values[1];
    return 1 + cos(M_PI * x * y);
}

double AlgoVariables::calc_B_val(const point2d<int> &position, const point2d<double> &coord,
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

BlockGrid2d &AlgoVariables::init_B(BlockGrid2d &block) {
    for (int i = block.bstart.values[0]; i < block.bend.values[0]; i++) {
        for (int j = block.bstart.values[1]; j < block.bend.values[1]; j++) {
            auto position = point2d<int>(i, j);
            auto coord = block.get_coord(position);
            double value = calc_B_val(position, coord, block.step, block.gsize);
            block.set_value(position, value, true);
        }
    }
    return block;
}

BlockGrid2d &AlgoVariables::init_k(BlockGrid2d &block, int axis) {
    for (int i = block.bstart.values[0]; i < block.bend.values[0]; i++) {
        for (int j = block.bstart.values[1]; j < block.bend.values[1]; j++) {
            auto position = point2d<int>(i, j);
            auto coord = block.get_coord(position);
            coord.values[axis] -= 0.5 * block.step.values[axis];
            double value = k_func(coord);
            block.set_value(position, value, true);
        }
    }
    return block;
}

BlockGrid2d &AlgoVariables::init_w_test(BlockGrid2d &block) {
    for (int i = block.bstart.values[0]; i < block.bend.values[0]; i++) {
        for (int j = block.bstart.values[1]; j < block.bend.values[1]; j++) {
            auto position = point2d<int>(i, j);
            auto coord = block.get_coord(position);
            block.set_value(position, u_func(coord), true);
        }
    }
    return block;
}

BlockGrid2d &AlgoVariables::init_w(BlockGrid2d &block) {
    for (int i = block.bstart.values[0]; i < block.bend.values[0]; i++) {
        for (int j = block.bstart.values[1]; j < block.bend.values[1]; j++) {
            auto position = point2d<int>(i, j);
            block.set_value(position, 0, true);
        }
    }
    return block;
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
                b.get_value(position) * w.grad(position, 1, false)
                + w.get_value(position)
        ) - w.aw_grad(position, a, 0);
    } else if ((i >= 1) && (i <= N - 1) && (j == 0)) { // bottom (eq. 9)
        auto position_i1 = point2d<int>(i, 1);
        return 2 / ystep * (
                -b.get_value(position_i1) * w.grad(position_i1, 1, false)
                + 0 * w.get_value(position)
        ) - w.aw_grad(position, a, 0);
    } else if ((i == 0) && (j >= 1) && (j <= M - 1)) { // left (eq. 8)
        auto position_1j = point2d<int>(1, j);
        return 2 / xstep * (
                -a.get_value(position_1j) * w.grad(position_1j, 0, false)
                + 0 * w.get_value(position)
        ) - w.aw_grad(position, b, 1);
    } else if ((i == N) && (j >= 1) && (j <= M - 1)) { // right (eq. 8)
        return 2 / xstep * (
                a.get_value(position) * w.grad(position, 0, false)
                + 1 * w.get_value(position)
        ) - w.aw_grad(position, b, 1);
    }

        //Corners (eq 10-13)
    else if ((i == 0) && (j == 0)) { // LB corner
        auto p00 = point2d<int>(0, 0);
        auto p10 = point2d<int>(1, 0);
        auto p01 = point2d<int>(0, 1);
        return
                -2 / xstep * a.get_value(p10) * w.grad(p10, 0, false) +
                -2 / ystep * b.get_value(p01) * w.grad(p01, 1, false) +
                (0 * 2 / xstep + 0 * 2 / ystep) * w.get_value(p00);
    } else if ((i == N) && (j == 0)) { // RB corner
        auto pM0 = point2d<int>(N, 0);
        auto pM1 = point2d<int>(N, 1);
        return 2 / xstep * a.get_value(pM0) * w.grad(pM0, 0, false) +
               -2 / ystep * b.get_value(pM1) * w.grad(pM1, 1, false) +
               (1 * 2 / xstep + 0 * 2 / ystep) * w.get_value(pM0) * 0;
    } else if ((i == N) && (j == M)) { // RT corner
        return 2 / xstep * a.get_value(position) * w.grad(position, 0, false) +
               2 / ystep * b.get_value(position) * w.grad(position, 1, false) +
               (1 * 2 / xstep + 1 * 2 / ystep) * w.get_value(position) * 0;
    } else if ((i == 0) && (j == M)) { // LT corner
        auto p0N = point2d<int>(0, M);
        auto p1N = point2d<int>(1, M);
        return -2 / xstep * a.get_value(p1N) * w.grad(p1N, 0, false) +
               2 / ystep * b.get_value(p0N) * w.grad(p0N, 1, false) +
               (0 * 2 / xstep + 1 * 2 / ystep) * w.get_value(p0N) * 0;

    }

    return 0;
}

BlockGrid2d &Aw(BlockGrid2d &Aw, const BlockGrid2d &w, const BlockGrid2d &a, const BlockGrid2d &b) {
    for (int i = Aw.bstart.values[0]; i < Aw.bend.values[0]; i++) {
        for (int j = Aw.bstart.values[1]; j < Aw.bend.values[1]; j++) {
            auto position = point2d<int>(i, j);
            double value = calc_Aw_val(position, w, a, b);
            Aw.set_value(position, value, true);
        }
    }
    return Aw;
}