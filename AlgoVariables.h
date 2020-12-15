//
// Created by d.bibik on 15/12/2020.
//



#ifndef SUPERCOMP_MATPHIS_ALGOVARIABLES_H
#define SUPERCOMP_MATPHIS_ALGOVARIABLES_H

#include "BlockContainer.h"
#include "BlockGrid2d.h"
#include <cmath>

class AlgoVariables {
public:
    BlockContainer B, w, w_test, a, b, Aw_val, r, Ar_val;

    AlgoVariables(point2d<int> gsize, point2d<double> gstart, point2d<double> gstep,
                  point2d<int> num_blocks);

private:
    static double F_func(point2d<double> point);

    static double psi_R_func(point2d<double> point);

    static double psi_L_func(point2d<double> point);

    static double psi_T_func(point2d<double> point);

    static double psi_B_func(point2d<double> point);

    static double k_func(point2d<double> point);

    static double u_func(point2d<double> point);

    static double calc_B_val(const point2d<int> &position, const point2d<double> &coord,
                             const point2d<double> &step, const point2d<int> &gsize);

    static BlockGrid2d &init_B(BlockGrid2d &block);

    static BlockGrid2d &init_k(BlockGrid2d &block, int axis);

    static BlockGrid2d &init_w_test(BlockGrid2d &block);

    static BlockGrid2d &init_w(BlockGrid2d &block);
};

double calc_Aw_val(const point2d<int> &position, const BlockGrid2d &w, const BlockGrid2d &a, const BlockGrid2d &b);

BlockGrid2d &Aw(BlockGrid2d &Aw, const BlockGrid2d &w, const BlockGrid2d &a, const BlockGrid2d &b);

#endif //SUPERCOMP_MATPHIS_ALGOVARIABLES_H
