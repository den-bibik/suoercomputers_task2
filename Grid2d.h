//
// Created by d.bibik on 10/12/2020.
//

#ifndef SUPERCOMP_MATPHIS_GRID2D_H
#define SUPERCOMP_MATPHIS_GRID2D_H

#endif //SUPERCOMP_MATPHIS_GRID2D_H

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




class Grid2d {
public:
    Grid2d(point2d<int> _gsize, point2d<double> _border, point2d<double> _step);
    Grid2d(const Grid2d &b);
    ~Grid2d();


    point2d<double> get_coord(point2d<int> point) const;
    virtual double get_value(point2d<int> p) const;
    virtual void set_value(point2d<int> position, double val);

    Grid2d &operator=(const Grid2d &b);
    void dump() const;
    virtual Grid2d operator-(const Grid2d &b) const;
    virtual Grid2d operator*(const double C) const;
    virtual Grid2d operator*(const Grid2d &b) const;


    void print_info() const;

    virtual double dot_prod(const Grid2d &b) const;
    double sum() const;

    double grad(point2d<int> p, int axis, bool forward_shift) const;
    double aw_grad(point2d<int> p, const Grid2d *a, int axis) const;
    double laplass(point2d<int> p, const Grid2d *a, const Grid2d *b) const;

    point2d<int> gsize;
    point2d<double> border, step;
protected:
    double *values;

};