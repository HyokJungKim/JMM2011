#include "NumIntegrate.h"
#include <random>
#include <iostream>

typedef std::vector<double> vd;
typedef std::vector<vd> v2d;

double test_fun(const vd& in_vec) {
    // Situation similar to the final period of consumption and savings problem
    return std::log(in_vec[0] + 5.0 * (1.0+in_vec[1]));
}

int main() {

    std::normal_distribution<double> dist_normal(0.0, 1.0);
    std::default_random_engine gen_random;

    const vd vcv {1.0*1.0, 0.0,
                  0.0, 0.01*0.01};
    const int Nd = 2;

    // -----------------------------------------------------------------------------------------------------------------
    // Monte Carlo
    // -----------------------------------------------------------------------------------------------------------------
    const int Nsim = 1000000; // Number of points for Monte-Carlo
    double MC_val = 0.0;
    vd tempvec(Nd, 0.0);
    for (int nn = 0; nn < Nsim; nn++) {
        for (int qq = 0; qq < Nd; qq++) {
            tempvec[qq] = std::sqrt(vcv[qq + qq*Nd]) * dist_normal(gen_random);
        }
        MC_val += test_fun(tempvec) / static_cast<double>(Nsim);
    }

    std::cout << "Numerical Integration using Monte Carlo Method: " << MC_val <<"\n";

    // -----------------------------------------------------------------------------------------------------------------
    // Gauss-Hermite Quadrature
    // -----------------------------------------------------------------------------------------------------------------
    v2d epsi;
    vd wght;
    const int Qn = 5; // Number of integration node per dimension
    double GH_val = 0.0; // Result using the Gauss-Hermite quadrature

    vd epsi_flat;

    std::tie(epsi, wght) = GH_Quadrature(Qn, Nd, vcv);

    for (int ii = 0; ii < static_cast<int>(wght.size()); ii++) {
        GH_val += test_fun(epsi[ii]) * wght[ii];
    }

    std::cout << "Numerical Integration using Gauss-Hermite Quadrature: " << GH_val <<"\n";

    // -----------------------------------------------------------------------------------------------------------------
    // Monomial 1
    // -----------------------------------------------------------------------------------------------------------------
    double M1_val = 0.0; // Result using Monomials 1

    std::tie(epsi, wght) = Monomials_1(Nd, vcv);

    for (int ii = 0; ii < static_cast<int>(wght.size()); ii++) {
        M1_val += test_fun(epsi[ii]) * wght[ii];
    }

    std::cout << "Numerical Integration using Monomial 1: " << M1_val <<"\n";

    // -----------------------------------------------------------------------------------------------------------------
    // Monomial 2
    // -----------------------------------------------------------------------------------------------------------------
    double M2_val = 0.0; // Result using Monomials 2

    std::tie(epsi, wght) = Monomials_2(Nd, vcv);

    for (int ii = 0; ii < static_cast<int>(wght.size()); ii++) {

        M2_val += test_fun(epsi[ii]) * wght[ii];
    }

    std::cout << "Numerical Integration using Monomial 2: " << M2_val <<"\n";

	return 0;
}