//
// Created by d.bibik on 15/12/2020.
//

#ifndef SUPERCOMP_MATPHIS_BLOCKGRID2D_H
#define SUPERCOMP_MATPHIS_BLOCKGRID2D_H

#include <vector>




template <class T>
class point2d {
public:
    T values [2];
    point2d(){ values[0]=0; values[1]=0;}
    point2d (T x, T y){ values[0]=x; values[1]=y;}

    point2d<T> operator+(const point2d b) const{
        return point2d<T>(values[0] + b.values[0], values[1] + b.values[1]);
    }
    point2d<T>  operator-(const point2d b) const{
        return point2d<T>(values[0] - b.values[0], values[1] - b.values[1]);
    }
    template < typename U >
    point2d<T>  operator/(const point2d<U> b) const{
        return point2d<T>(values[0] / b.values[0], values[1] / b.values[1]);
    }
    template < typename U >
    point2d<T>  operator*(const point2d<U> b) const{
        return point2d<T>(values[0] * b.values[0], values[1] * b.values[1]);
    }
};



enum block_part {
    TOP, BOTTOM, LEFT, RIGHT, LT_CORNER, RT_CORNER, LB_CORNER, RB_CORNER
};

class BlockGrid2d {
public:
    BlockGrid2d(point2d<int> block_start, point2d<int> block_end, point2d<int> gsize,
                point2d<double> border, point2d<double> step);

    BlockGrid2d &operator=(const BlockGrid2d &b);

    BlockGrid2d(const BlockGrid2d &b);

    point2d<double> get_coord(point2d<int> point) const;

    double get_value(point2d<int> p) const;


    void set_external_border_part(block_part key, std::__1::vector<double> value);

    std::__1::vector<double> get_internal_border_part(block_part key);

    void set_value(point2d<int> p, double val, bool internal_op);

    friend BlockGrid2d operator-(const BlockGrid2d &a, const BlockGrid2d &b);

    friend BlockGrid2d operator*(const BlockGrid2d &a, const double C);

    friend BlockGrid2d operator*(const double C, const BlockGrid2d &a);

    double dot_prod(const BlockGrid2d &b) const;

    double grad(point2d<int> p, int axis, bool forward_shift) const;

    double aw_grad(point2d<int> p, const BlockGrid2d &a, int axis) const;

    double laplass(point2d<int> p, const BlockGrid2d &a, const BlockGrid2d &b) const;

    point2d<int> bstart, bsize, bend, gsize;
    point2d<double> border, step;
private:
    std::__1::shared_ptr<std::__1::vector<double>> values;
    static void copy_data(const BlockGrid2d &from, BlockGrid2d &to);

    bool check_internal(const point2d<int> point) const;

};

#endif //SUPERCOMP_MATPHIS_BLOCKGRID2D_H