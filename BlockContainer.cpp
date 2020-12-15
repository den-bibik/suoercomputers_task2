//
// Created by d.bibik on 15/12/2020.
//

#include "BlockGrid2d.h"
#include "BlockContainer.h"

BlockContainer::BlockContainer(point2d<int> gsize, point2d<double> gstart, point2d<double> gstep,
                               point2d<int> num_blocks) : num_blocks(num_blocks) {
    point2d<int> block_size = (gsize + point2d<int>(1, 1)) / num_blocks;

    auto block_start = point2d<int>(0, 0);
    for (int i = 0; i < num_blocks.values[0]; i++) {
        block_start.values[1] = 0;
        for (int j = 0; j < num_blocks.values[1]; j++) {
            auto block_end = block_start + block_size;
            auto block = BlockGrid2d(block_start, block_end, gsize, gstart, gstep);
            matrix.push_back(block);

            block_start.values[1] += block_size.values[1];
        }
        block_start.values[0] += block_size.values[0];
    }

}

void BlockContainer::set(int i, int j, const BlockGrid2d &val) {
    matrix[i * num_blocks.values[1] + j] = val;
}

BlockGrid2d &BlockContainer::get(int i, int j) {
    return matrix[i * num_blocks.values[1] + j];
}

void BlockContainer::sync_borders() {
    for (int i = 0; i <= num_blocks.values[0]; i++) {
        for (int j = 0; j <= num_blocks.values[0]; j++) {
            send_border_part(i, j, i - 1, j, block_part::TOP, block_part::BOTTOM);
            send_border_part(i, j, i + 1, j, block_part::BOTTOM, block_part::TOP);
            send_border_part(i, j, i, j - 1, block_part::LEFT, block_part::RIGHT);
            send_border_part(i, j, i, j + 1, block_part::RIGHT, block_part::LEFT);

            send_border_part(i, j, i - 1, j - 1, block_part::LT_CORNER, block_part::RB_CORNER);
            send_border_part(i, j, i + 1, j - 1, block_part::LB_CORNER, block_part::RT_CORNER);
            send_border_part(i, j, i - 1, j + 1, block_part::RT_CORNER, block_part::LB_CORNER);
            send_border_part(i, j, i + 1, j + 1, block_part::RB_CORNER, block_part::LT_CORNER);
        }
    }
}

void BlockContainer::send_border_part(int from_i, int from_j, int to_i, int to_j,
                                      block_part from_key, block_part to_key) {
    if (check_index(from_i, from_j) && check_index(to_i, to_j)) {
        auto from = this->get(from_i, from_j);
        auto sync_data = from.get_internal_border_part(from_key);
        this->get(to_i, to_j).set_external_border_part(to_key, sync_data);
    }
}

bool BlockContainer::check_index(int i, int j) {
    return (i >= 0) && (i < num_blocks.values[0]) && (j >= 0) && (j < num_blocks.values[1]);
}
