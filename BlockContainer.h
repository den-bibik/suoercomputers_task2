//
// Created by d.bibik on 15/12/2020.
//


#ifndef SUPERCOMP_MATPHIS_BLOCKCONTAINER_H
#define SUPERCOMP_MATPHIS_BLOCKCONTAINER_H

#include "BlockGrid2d.h"
#include <vector>
class BlockContainer {
public:
    BlockContainer() = default;

    BlockContainer(point2d<int> gsize, point2d<double> gstart, point2d<double> gstep,
                   point2d<int> num_blocks);


    void set(int i, int j, const BlockGrid2d &val);

    BlockGrid2d &get(int i, int j);

    void sync_borders();

    point2d<int> num_blocks;
private:
    void send_border_part(int from_i, int from_j, int to_i, int to_j, block_part from_key, block_part to_key);

    inline bool check_index(int i, int j);

    std::vector<BlockGrid2d> matrix;
};

#endif //SUPERCOMP_MATPHIS_BLOCKCONTAINER_H

