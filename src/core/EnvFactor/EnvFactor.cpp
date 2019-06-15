#include "EnvFactor.h"

EnvFactor::EnvFactor(MPFM* run, int size, double h_val){
    this->mpfm = run;

    // check n patches


    this->fractal = this->generate_fractal(size, h_val);

    this->map = new double*[size];
    for (int i = 0; i < size; i++){
        this->map[i] = new double[size];
    }

    this->resolution = size;
    this->range_center = size/2;
    this->gradient_strength = 1.0;
    this->fractal_weight = 0.3;

    this->shift(0.5);
}

void EnvFactor::shift(double center_of_range){
    int n = this->resolution;

    double center_row = n*center_of_range;

    double E_C = 1.0;
    double dist, E_L, val;
    double F = this->fractal_weight;
    double G = this->gradient_strength;

    double max = 0.0;
    double min = 1.0;

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            dist = double(center_row - i) / double(n);
            E_L = E_C - G*abs(dist);
            val = E_L + F*(this->fractal[i][j]);

            if (val > max){
                max = val;
            }
            if (val < min){
                min = val;
            }
            this->map[i][j] = val;
        }
    }

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            this->map[i][j] = double(this->map[i][j] - min)/double(max-min);
        }
    }
}

double EnvFactor::get_env_factor_value(double x, double y){

    int res = this->mpfm->EF_GRANULARITY;
    int x_cell = int(double(res)*(x));
    int y_cell = int(double(res)*(y));
    return this->map[x_cell][y_cell];
}

double EnvFactor::get_cell_value(int x, int y){
    return this->map[x][y];
}

double** EnvFactor::generate_fractal(int size, double H_VAL){

    int n = size+1;
    double** tmp = new double*[n];
    for (int i = 0; i < n; i++){
        tmp[i] = new double[n];
    }

    double delta = 1.0;
    int N = size;

    int maxlevel = log2(size);

    double H = H_VAL;

    int D = N;
    int d = N/2;

    for (int i = 0; i <= N; i++){
        for (int j = 0; j <= N; j++){
            tmp[i][j] = 0.0;
        }
    }

    // Init corner values
    tmp[0][0] = delta*this->mpfm->normal(this->mpfm->ef_gen, 0.0, 1.0);
    tmp[0][N] = delta*this->mpfm->normal(this->mpfm->ef_gen, 0.0, 1.0);
    tmp[N][0] = delta*this->mpfm->normal(this->mpfm->ef_gen, 0.0, 1.0);
    tmp[N][N] = delta*this->mpfm->normal(this->mpfm->ef_gen, 0.0, 1.0);

    for (int stage = 0; stage < maxlevel; stage++){
        delta = delta * double(pow(double(0.5), double(0.5*H)));

        // interpolate and offset points
        for (int x = d; x <= N-d; x+=D){
            for (int y = d; y <= N-d; y+=D){
                tmp[x][y] = this->f4(delta, tmp[x+d][y+d], tmp[x+d][y-d], tmp[x-d][y+d],tmp[x-d][y-d]);
            }
        }

        delta = delta * double(pow(double(0.5), double(0.5*H)));

        // boundary grid points
        for (int x = d; x <= N-d; x+=D){
            tmp[x][0] = this->f3(delta, tmp[x+d][0], tmp[x-d][0], tmp[x][d]);
            tmp[x][N] = this->f3(delta, tmp[x+d][N], tmp[x-d][N], tmp[x][N-d]);
            tmp[0][x] = this->f3(delta, tmp[0][x+d], tmp[0][x-d], tmp[d][x]);
            tmp[N][x] = this->f3(delta, tmp[N][x+d], tmp[N][x-d], tmp[N-d][x]);
        }

        // interpolate and offset interior grid points
        for (int x = d; x <= N-d; x+=D){
            for (int y = D; y <= N-d; y+=D){
                tmp[x][y] = this->f4(delta, tmp[x][y+d], tmp[x][y-d], tmp[x+d][y], tmp[x-d][y]);
            }
        }

        for (int x = D; x <= N-d; x+=D){
            for (int y = d; y <= N-d; y+=D){
                tmp[x][y] = this->f4(delta, tmp[x][y+d], tmp[x][y-d], tmp[x+d][y], tmp[x-d][y]);
            }
        }

        D = D/2;
        d = d/2;
    }


    return tmp;
    /*
    double** grid = new double*[n];
    for (int i = 0; i < n; i++){
        grid[i] = new double[n];
    }


    double min = tmp[0][0];
    double max = tmp[0][0];


    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (tmp[i][j] > max){
                max = tmp[i][j];
            }
            if (tmp[i][j] < min){
                min = tmp[i][j];
            }
        }
    }

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            grid[i][j] = double(max - tmp[i][j])/(max-min);
        }
    }


    for (int i = 0; i < n; i++){
        delete tmp[i];
    }
    delete[] tmp;

    return grid;*/
}

double EnvFactor::f4(double delta, double a, double b, double c, double d){
    double sum = a+b+c+d;
    double avg = sum/4.0;
    double val =  avg + double(delta) * double(this->mpfm->normal(this->mpfm->ef_gen, 0.0, 1.0));
    return val;
}

double EnvFactor::f3(double delta, double a, double b, double c){
    double sum = a+b+c;
    double avg = sum/3.0;
    double val =  avg + double(delta) * double(this->mpfm->normal(this->mpfm->ef_gen, 0.0, 1.0));
    return val;
}
