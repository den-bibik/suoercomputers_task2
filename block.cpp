//
// Created by d.bibik on 30/11/2020.
//
#include <iostream>

#include "BlockGrid2d.h"
#include "AlgoVariables.h"


void easy_iteration_method(AlgoVariables vars, int num_iter){
    for (int iter = 0; iter < num_iter; iter++) {
        double alpha_numerator = 0;
        double alpha_denominator = 0;
        double r_norm = 0;
        for (int i = 0; i < vars.B.num_blocks.values[0]; i++) {
            for (int j = 0; j < vars.B.num_blocks.values[1]; j++) {
                auto B_block = vars.B.get(i, j);
                auto a_block = vars.a.get(i, j);
                auto b_block = vars.b.get(i, j);
                auto Aw_val_block = vars.Aw_val.get(i, j);
                auto w_block = vars.w.get(i, j);
                auto r_block = vars.r.get(i, j);

                Aw_val_block = Aw(Aw_val_block, w_block, a_block, b_block);
                r_block = Aw_val_block - B_block;
                vars.r.set(i, j, r_block);
            }
        }
        vars.r.sync_borders();

        for (int i = 0; i < vars.B.num_blocks.values[0]; i++) {
            for (int j = 0; j < vars.B.num_blocks.values[1]; j++) {
                auto Ar_val_block = vars.Ar_val.get(i, j);
                auto r_block = vars.r.get(i, j);
                auto a_block = vars.a.get(i, j);
                auto b_block = vars.b.get(i, j);

                Ar_val_block = Aw(Ar_val_block, r_block, a_block, b_block);
                //sync: (sum) alpha_numerator, alpha_denominator, r_norm
                alpha_numerator += Ar_val_block.dot_prod(r_block);
                alpha_denominator += Ar_val_block.dot_prod(Ar_val_block);
                r_norm += r_block.dot_prod(r_block);
            }
        }
        double alpha = alpha_numerator / alpha_denominator;
        double stop_norm = alpha * alpha * r_norm;
        double diff_val = 0;
        for (int i = 0; i < vars.B.num_blocks.values[0]; i++) {
            for (int j = 0; j < vars.B.num_blocks.values[1]; j++) {
                //need r_block, w_block, w_test_block
                auto w_block = vars.w.get(i, j);
                auto w_test_block = vars.w_test.get(i, j);
                auto r_block = vars.r.get(i, j);

                auto alg_step = r_block * alpha;

                auto diff = w_block - w_test_block;
                diff_val += diff.dot_prod(diff);

                w_block = w_block - alg_step;

                vars.w.set(i, j, w_block);

            }
        }
        vars.w.sync_borders();
        std::cout << "iter: " << iter << " stop criterion val: " << stop_norm <<
                  " test: " << diff_val << std::endl;

    }
}

void block_algo(int size, int x_blocks, int y_blocks, int num_iter) {
    auto gsize = point2d<int>(size - 1, size - 1);
    auto border = point2d<double>(0, 0);
    auto border2 = point2d<double>(2, 1);
    auto step = (border2 - border) / gsize;
    auto num_blocks = point2d<int>(x_blocks, y_blocks);

    AlgoVariables vars(gsize, border, step, num_blocks);

    easy_iteration_method(vars, num_iter);
}

int main() {
    std::cout << "It's block version" << std::endl;
    block_algo(100, 50, 50, 100);
    //algo(20, 100, 1, 0);
    return 0;
}