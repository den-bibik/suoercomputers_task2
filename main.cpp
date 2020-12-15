//
// Created by d.bibik on 27/11/2020.
//

#include <iostream>
#include <cmath>
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


double F_func(point2d<double> point){
    // calculated from F(x,y) = -laplass(u) + q(x, y) * u,  u = 1 + cos(pi * x * y), k = 4 + x + y, z = 0
    double x = point.values[0]; double y = point.values[1];
    return M_PI * (M_PI * (x + y + 4) * (x * x + y * y) * cos(M_PI * x * y) + (x + y)*sin(M_PI * x * y));
}

double psi_R_func(point2d<double> point){
    double x = point.values[0]; double y = point.values[1];
    return -M_PI * y * (x + y + 4) * sin(M_PI * x * y) + cos(M_PI * x * y) + 1;
}

double psi_L_func(point2d<double> point){
    double x = point.values[0]; double y = point.values[1];
    return -M_PI * y * (x + y + 4) * sin(M_PI * x * y);
}

double psi_T_func(point2d<double> point){
    double x = point.values[0]; double y = point.values[1];
    return -M_PI * x * (x + y + 4) * sin(M_PI * x * y) + cos(M_PI * x * y) + 1;
}

double psi_B_func(point2d<double> point){
    double x = point.values[0]; double y = point.values[1];
    return -M_PI * x * (x + y + 4) * sin(M_PI * x * y);
}

double k_func(point2d<double> point){
    double x = point.values[0]; double y = point.values[1];
    return 4 + x + y;
}

double u_func(point2d<double> point){
    double x = point.values[0]; double y = point.values[1];
    return 1 + cos(M_PI*x*y);
}


class Grid2d{
public:
    Grid2d(point2d<int> _gsize, point2d<double> _border, point2d<double>_step)
    {
        gsize = _gsize;
        border = _border;
        step = _step;

        int values_size = (gsize.values[0] + 1) * (gsize.values[1] + 1);
        values = new double[values_size];
    }

    Grid2d(const Grid2d& b)
    {
        gsize = b.gsize;
        border = b.border;
        step = b.step;

        int values_size = (gsize.values[0] + 1) * (gsize.values[1] + 1);
        values = new double[values_size];
        memcpy(values, b.values, sizeof(*b.values)*values_size);
    }



    ~Grid2d(){
        //std::cout << "destructor val " << values << std::endl;
        delete[] values;
    }

    Grid2d& operator=(const Grid2d& b)
    {
        gsize = b.gsize;
        border = b.border;
        step = b.step;

        int values_size = (gsize.values[0] + 1) * (gsize.values[1] + 1);
        //values = new double[values_size];
        memcpy(values, b.values, sizeof(*b.values)*values_size);

        return *this;
    }

    point2d<double> get_coord(point2d<int> point) const{
        return border + step * point;
    }

    double get_value(point2d<int> p) const{
        if((p.values[0] < 0)||(p.values[1] < 0)){p.values[0] = 1/0;}
        return values[p.values[0] * (gsize.values[1] + 1) + p.values[1]];
    }

    void set_value(point2d<int> position, double val){
        if((position.values[0] < 0)||(position.values[1] < 0)){position.values[0] = 1/0;}
        values[position.values[0] * (gsize.values[1] + 1) + position.values[1]] = val;
    }

    void dump() const{
        std::cout << "matrix dump" << std::endl;
        for(int j = gsize.values[1]; j >= 0; j--) {
            for (int i = 0; i <= gsize.values[0]; i++) {
                auto p = point2d<int>(i, j);
                std::cout << get_value(p) << " ";
            }
            std::cout << std::endl;
        }
    }

    Grid2d operator-(const Grid2d& b) const{
        Grid2d res = Grid2d(gsize, border, step);
        for(int i = 0; i <= gsize.values[0]; i++) {
            for (int j = 0; j <= gsize.values[1]; j++) {
                auto p = point2d<int>(i, j);
                res.set_value(p, get_value(p) - b.get_value(p));
            }
        }
        return res;
    }


    Grid2d operator*(const double C) const{
        Grid2d res = Grid2d(gsize, border, step);
        for(int i = 0; i <= gsize.values[0]; i++) {
            for (int j = 0; j <= gsize.values[1]; j++) {
                auto p = point2d<int>(i, j);
                res.set_value(p, get_value(p) * C);
            }
        }
        return res;
    }

    Grid2d operator*(const Grid2d& b) const{
        Grid2d res = Grid2d(gsize, border, step);
        for(int i = 0; i <= gsize.values[0]; i++) {
            for (int j = 0; j <= gsize.values[1]; j++) {
                auto p = point2d<int>(i, j);
                res.set_value(p, get_value(p) * b.get_value(p));
            }
        }
        return res;
    }


    void print_info() const{
        double max = -1000000000;
        double min = 10000000000;
        double sum = 0.0;
        for(int i = 0; i <= gsize.values[0]; i++) {
            for (int j = 0; j <= gsize.values[1]; j++) {
                auto p = point2d<int>(i, j);
                double val = get_value(p);
                if(val > max) max = val;
                if(val < min) min = val;
                sum += val;

            }
        }
        sum /= (gsize.values[0] + 1) * (gsize.values[1] + 1);

        std::cout << " mean " << sum << " max " << max << " min " << min << std::endl;
    }

    double dot_prod(const Grid2d& b) const{
        int values_size = gsize.values[0] * gsize.values[1];
        double res = 0;
        for(int i = 0; i <= gsize.values[0]; i++) {
            for (int j = 0; j <= gsize.values[1]; j++) {
                auto p = point2d<int>(i, j);
                double val = get_value(p) * b.get_value(p);
                if ((i == 0) | (i == gsize.values[0])) val /= 2;
                if ((j == 0) | (j == gsize.values[1])) val /= 2;
                res += val;
            }
        }
        res *=  step.values[0] * step.values[1];
        return res;
    }

    double sum() const{
        double res = 0;
        for(int i = 0; i <= gsize.values[0]; i++) {
            for (int j = 0; j <= gsize.values[1]; j++) {
                auto p = point2d<int>(i, j);
                double val = get_value(p);
                res += val;
            }
        }
        res *=  step.values[0] * step.values[1];
        return res;
    }

    double grad(point2d<int> p, int axis, bool forward_shift) const{
        auto b = point2d<int>(p.values[0], p.values[1]);
        if(forward_shift) b.values[axis] += 1; else b.values[axis] -= 1;
        double result = (get_value(p) - get_value(b)) / step.values[axis];
        if(forward_shift) result *= -1;

        return result;
    }

    double aw_grad(point2d<int> p, const Grid2d& a, int axis) const{
        auto p2 = point2d<int>(p.values[0], p.values[1]);
        p2.values[axis] += 1;
        double res =(
                    a.get_value(p2) * grad(p2, axis, 0) -
                    a.get_value(p) * grad(p, axis, 0)
                ) / step.values[axis];
        return res;
    }

    double laplass(point2d<int> p, const Grid2d& a, const Grid2d& b) const{
        double res = aw_grad(p, a, 0) + aw_grad(p, b, 1);
        return res;
    }

    point2d<int> gsize;
    point2d<double> border, step;
private:
    double* values;

};

Grid2d Aw(const Grid2d& w, const Grid2d& a, const Grid2d& b){
    Grid2d res = Grid2d(w.gsize, w.border, w.step);
    // internal points
    for(int i=1; i<=res.gsize.values[0]-1; i++){
        for(int j=1; j<=res.gsize.values[1]-1; j++){
            auto position = point2d<int>(i, j);
            double value = -w.laplass(position, a, b);
            res.set_value(position, value);
        }
    }

    double mul = 1.0;
    //top (eq. 9)
    int j = res.gsize.values[1];
    for(int i=1; i<=res.gsize.values[0]-1; i++){
        auto position = point2d<int>(i, j);
        double step = w.step.values[1];
        double value = 2/step * (
                    b.get_value(position) * w.grad(position, 1, 0)
                    + w.get_value(position)
                ) - w.aw_grad(position, a, 0);
        res.set_value(position, value * mul);
    }

    //bottom (eq. 9)
    j = 0;
    for(int i=1; i<=res.gsize.values[1]-1; i++){
        auto position = point2d<int>(i, j);
        auto position_i1 = point2d<int>(i, j + 1);
        double step = w.step.values[1];
        double value = 2/step * (
                -b.get_value(position_i1) * w.grad(position_i1, 1, 0)
                + 0 * w.get_value(position)
        ) - w.aw_grad(position, a, 0);
        res.set_value(position, value * mul);
    }

    //left (eq. 8)
    int i = 0;
    for(j=1; j<=res.gsize.values[0]-1; j++){
        auto position = point2d<int>(i, j);
        auto position_1j = point2d<int>(i + 1, j);
        double step = w.step.values[0];
        double value = 2/step * (
                -a.get_value(position_1j) * w.grad(position_1j, 0, 0)
                + 0 * w.get_value(position)
        ) - w.aw_grad(position, b, 1);
        res.set_value(position, value* mul);
    }
    //right (eq. 8)

    i = res.gsize.values[1];
    for(j=1; j<=res.gsize.values[0]-1; j++){
        auto position = point2d<int>(i, j);
        double step = w.step.values[0];
        double value = 2/step * (
                a.get_value(position) * w.grad(position, 0, 0)
                + 1 * w.get_value(position)
        ) - w.aw_grad(position, b, 1);
        res.set_value(position, value* mul);
    }

    //eq 10-13
    double val;
    auto p00 = point2d<int>(0, 0);
    auto p10 = point2d<int>(1,0);
    auto p01 = point2d<int>(0,1);
    val =
            -2/res.step.values[0] * a.get_value(p10) * w.grad(p10, 0, 0) +
            -2/res.step.values[1] * b.get_value(p01) * w.grad(p01, 1, 0) +
            (0* 2/res.step.values[0] + 0*2/res.step.values[1]) * w.get_value(p00);
    res.set_value(p00, val* mul);


    auto pM0 = point2d<int>(res.gsize.values[0], 0);
    auto pM1 = point2d<int>(res.gsize.values[0],1);
    val =
            2/res.step.values[0] * a.get_value(pM0) * w.grad(pM0, 0, 0) +
            -2/res.step.values[1] * b.get_value(pM1) * w.grad(pM1, 1, 0) +
            (1* 2/res.step.values[0] + 0*2/res.step.values[1]) * w.get_value(pM0)*0;
    res.set_value(pM0, val* mul);

    auto pMN = point2d<int>(res.gsize.values[0],res.gsize.values[1]);
    val =
            2/res.step.values[0] * a.get_value(pMN) * w.grad(pMN, 0, 0) +
            2/res.step.values[1] * b.get_value(pMN) * w.grad(pMN, 1, 0) +
            (1* 2/res.step.values[0] + 1*2/res.step.values[1]) * w.get_value(pMN)*0;
    res.set_value(pMN, val* mul);

    auto p0N = point2d<int>(0,res.gsize.values[1]);
    auto p1N = point2d<int>(1,res.gsize.values[1]);
    val =
            -2/res.step.values[0] * a.get_value(p1N) * w.grad(p1N, 0, 0) +
            2/res.step.values[1] * b.get_value(p0N) * w.grad(p0N, 1, 0) +
            (0*2/res.step.values[0] + 1*2/res.step.values[1]) * w.get_value(p0N)*0;
    res.set_value(p0N, val* mul);



    return res;
}

Grid2d init_B(point2d<int> _gsize, point2d<double> _border, point2d<double> _step){
    Grid2d res = Grid2d(_gsize, _border, _step);

    // internal points
    for(int i=1; i<=res.gsize.values[0]-1; i++){
        for(int j=1; j<=res.gsize.values[1]-1; j++){
            auto position = point2d<int>(i, j);
            double value = F_func(res.get_coord(position));
            res.set_value(position, value );
        }
    }

    double mul = 1.0; //res.step.values[1] / 100;
    //top (eq. 9)
    int j = res.gsize.values[1];
    for(int i=1; i<=res.gsize.values[0]-1; i++){
        auto position = point2d<int>(i, j);
        double step = res.step.values[1];
        auto coord = res.get_coord(position);
        double value = F_func(coord) + 2 / step * psi_T_func(coord);
        res.set_value(position, value * mul);
    }

    //bottom (eq. 9)
    j = 0;
    for(int i=1; i<=res.gsize.values[1]-1; i++){
        auto position = point2d<int>(i, j);
        double step = res.step.values[1];
        auto coord = res.get_coord(position);
        double value = F_func(coord) + 2 / step * psi_B_func(coord);
        res.set_value(position, value * mul);
    }

    //left (eq. 8)
    int i = 0;
    for(j=1; j<=res.gsize.values[0]-1; j++){
        auto position = point2d<int>(i, j);
        double step = res.step.values[0];
        auto coord = res.get_coord(position);
        double value = F_func(coord) + 2 / step * psi_L_func(coord);
        res.set_value(position, value * mul);
    }

    //right (eq. 8)
    i = res.gsize.values[1];
    for(j=1; j<=res.gsize.values[0]-1; j++){
        auto position = point2d<int>(i, j);
        double step = res.step.values[0];
        auto coord = res.get_coord(position);
        double value = F_func(coord) + 2 / step * psi_R_func(coord);
        res.set_value(position, value * mul);
    }

    // TODO: eq 10-13
    double coef_x = 2 / res.step.values[0];
    double coef_y =  2 / res.step.values[1];

    auto position = point2d<int>(0, 0);
    auto coord = res.get_coord(position);
    double value = F_func(coord); // psi_L(0,0) ==  psi_B(0,0) == 0
    res.set_value(position, value * mul);

    position = point2d<int>(res.gsize.values[0], 0);
    coord = res.get_coord(position);
    value = F_func(coord) + coef_x * psi_R_func(coord) + coef_y * psi_B_func(coord);
    res.set_value(position, value * mul);

    position = point2d<int>(0, res.gsize.values[1]);
    coord = res.get_coord(position);
    value = F_func(coord) + coef_x * psi_L_func(coord) + coef_y * psi_T_func(coord);
    res.set_value(position, value * mul);

    position = point2d<int>(res.gsize.values[0], res.gsize.values[1]);
    coord = res.get_coord(position);
    value = F_func(coord) + coef_x * psi_R_func(coord) + coef_y * psi_T_func(coord);
    res.set_value(position, value * mul);


    return res;
}

Grid2d init_k_grid(point2d<int> _gsize, point2d<double> _border, point2d<double> _step, int axis){
    Grid2d res = Grid2d(_gsize, _border, _step);
    for(int i=0; i<=res.gsize.values[0]; i++){
        for(int j=0; j<=res.gsize.values[1]; j++){
            auto position = point2d<int>(i, j);
            auto coord = res.get_coord(position);
            coord.values[axis] -= 0.5 * res.step.values[axis];
            double value = k_func(coord);
            res.set_value(position, value);
        }
    }
    return res;
}

Grid2d init_w_test_grid(point2d<int> _gsize, point2d<double> _border, point2d<double> _step){
    Grid2d res = Grid2d(_gsize, _border, _step);
    for(int i=0; i<=res.gsize.values[0]; i++){
        for(int j=0; j<=res.gsize.values[1]; j++){
            auto position = point2d<int>(i, j);
            auto coord = res.get_coord(position);
            double value = u_func(coord);
            res.set_value(position, value);
        }
    }
    return res;
}

Grid2d init_w_grid(point2d<int> _gsize, point2d<double> _border, point2d<double> _step){
    Grid2d res = Grid2d(_gsize, _border, _step);
    for(int i=0; i<=res.gsize.values[0]; i++){
        for(int j=0; j<=res.gsize.values[1]; j++){
            auto position = point2d<int>(i, j);
            auto coord = res.get_coord(position);
            double value = u_func(coord) * 0.1;
            res.set_value(position, 0 );

        }
    }
    return res;
}


void algo(int size, int num_iter, bool use_CG, bool debug){
    auto gsize = point2d<int>(size, size);
    auto border = point2d<double>(0, 0);
    auto border2 = point2d<double>(2, 1);
    auto step = (border2 - border) / gsize;

    auto w = init_w_grid(gsize, border, step);
    if(debug){
        w = init_w_test_grid(gsize, border, step);
        std::cout << "w" << std::endl; w.dump(); std::cout << std::endl;
    }
    auto w_test = init_w_test_grid(gsize, border, step);
    auto a = init_k_grid(gsize, border, step, 0);
    auto b = init_k_grid(gsize, border, step, 1);
    auto B = init_B(gsize, border, step);
    auto Aw_val = Aw(w, a, b);
    if(debug){
        auto B = Aw(w_test, a, b);
        std::cout << "B" << std::endl; B.dump(); std::cout << std::endl;
    }

    if(use_CG){

        auto r = B - Aw_val;
        auto z = r;
        auto diff = w - w_test;
        if(debug) std::cout << "diff_norm" <<diff.dot_prod(diff) << std::endl;

        for(int i = 0; i<num_iter; i++){
            auto Az = Aw(z, a, b);
            double alpha = r.dot_prod(r) / Az.dot_prod(z);
            auto alg_step = z * (-alpha);
            if(debug) {
                //std::cout << "z norm " << z.dot_prod(z) << "z sum " << z.sum() <<  " "; z.print_info();
                auto err = alg_step * (alg_step - (w - w_test) * 2);
                //std::cout << "err norm " << err.dot_prod(err) << "err sum " << err.sum() <<  " "; err.print_info();
                //err.dump();
                double lr = 1.0;
                /*auto diff_new = w - alg_step * lr -w_test;
                r.print_info();
                while (diff_new.dot_prod(diff_new) > diff.dot_prod(diff)){
                    lr /= 2.0;
                    diff_new = w - alg_step * lr -w_test;
                }*/
                auto lin_diff = Aw(w - w_test, a, b);
                //std::cout << "lin_diff "; lin_diff.print_info(); (w - w_test).print_info();
            }
            w = w - alg_step;

            auto r_next = r - Az * alpha;
            double betta = r_next.dot_prod(r_next) / r.dot_prod(r);
            if(debug) {
                //if(betta > 1) betta = 1;
            }
            r = r_next;
            z = r - z * (-betta);

            diff = w - w_test;
            if(debug){
                //std::cout << i  << " lr: " << lr << " betta: " << betta <<" alg_step: " << alg_step.dot_prod(alg_step) ;
            }
            std::cout << i  <<" alg_step: " << alg_step.dot_prod(alg_step) <<" test: " << diff.dot_prod(diff) <<std::endl;
        }
    }
    else{
        double stop_norm = 10000000;
        for(int i = 0; i<num_iter; i++){
            Aw_val = Aw(w, a, b);
            auto r = Aw_val - B;
            auto Ar_val = Aw(r, a, b);


            double alpha_numerator = Ar_val.dot_prod(r);
            double alpha_denominator = Ar_val.dot_prod(Ar_val);
            double alpha = alpha_numerator / alpha_denominator;
            auto alg_step = r * alpha;


            if(debug){
                std::cout << " r_norm " << r.dot_prod(r)  << " w_norm " << w.dot_prod(w)<< std::endl;
                std::cout << "Ar_val norm " << Ar_val.dot_prod(Ar_val)
                          << " Ar_val.dot_prod(r)  " << Ar_val.dot_prod(r) << std::endl;

                // ~=  ||(w - alg_step) - w_test||^2 -||w - w_test||^2
                auto err = alg_step * (alg_step -  (w - w_test)  * 2);
                std::cout << "err norm " << err.dot_prod(err) << "err sum " << err.sum() <<  " "; err.print_info();
                err.dump();
            }


            stop_norm = alg_step.dot_prod(alg_step);
            auto diff = w - w_test;
            w = w - alg_step;

            std::cout << "stop criterion val: " << stop_norm <<
                      " test: " << diff.dot_prod(diff) <<std::endl;
            if(debug) {
                std::cout << "diff info";
                diff.print_info();


                /*std::cout << "Aw_val" << std::endl; Aw_val.dump(); std::cout << std::endl;
                std::cout << "w" << std::endl; w.dump(); std::cout << std::endl;
                std::cout << "B" << std::endl; B.dump(); std::cout << std::endl;
                std::cout << "r" << std::endl; r.dump(); std::cout << std::endl;
                std::cout << "Ar" << std::endl; Ar_val.dump(); std::cout << std::endl;*/
            }
        }
    }
}

int main() {
    std::cout << "Hello, World!" << std::endl;
    algo(100, 100, 0, 0);
    return 0;
}

