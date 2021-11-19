//
// Created by hyokzzang on 12/27/20.
//

#include<vector>
#include<cmath>
#include<tuple>

#ifndef CPP_EXP_NUMINTEGRATE_H
#define CPP_EXP_NUMINTEGRATE_H

const double PI_VAL = 3.141592653589793;

/*
 *  Custom built Cholesky decomposition
 */
std::vector<double> chol(const int nDim, const std::vector<double> &vcv){

    int matsize = nDim * nDim;

    std::vector<double> LL(matsize, 0.0);

    double ss;

    for (int ii = 0; ii < nDim; ii++){
        for (int jj = 0; jj <= ii; jj++){
            ss = 0.0;

            if (jj == ii){
                for (int kk = 0; kk < jj; kk++) {
                    ss += LL[jj*nDim + kk] * LL[jj*nDim + kk];
                }

                LL[ii*nDim + jj] = sqrt(vcv[jj*nDim + jj] - ss);
            }
            else {
                for (int kk = 0; kk < jj; kk++) {
                    ss += LL[ii*nDim + kk] * LL[jj*nDim + kk];
                }

                LL[ii*nDim + jj] = (vcv[ii*nDim + jj] - ss)/LL[jj*nDim + jj];
            }
        }
    }

    // Returns upper triangular matrix!!!
    std::vector<double> UU(matsize, 0.0);
    for (int ii = 0; ii < nDim; ii++){
        for (int jj = ii; jj < nDim; jj++){
            UU[ii*nDim + jj] = LL[jj*nDim + ii];
        }
    }

    return UU;
}

/*
 *  Matrix multiplication
 */
std::vector<double> matmul(const int r1, const int c1, const std::vector<double> &Mat1,
                           const int c2, const std::vector<double> &Mat2){

    std::vector<double> outMat(r1*c2, 0.0);

    for (int ii = 0; ii < r1; ii++){
        for (int jj = 0; jj < c2; jj++){
            for (int idx = 0; idx < c1; idx++){
                outMat[ii*c2 + jj] += Mat1[ii*c1+idx] * Mat2[idx*c2+jj];
            }
        }
    }

    return outMat;
}

/*
 * Monomial rule version 1
 *
 * Algorithm adapted from
 * Solving Dynamic Economic Models" by Kenneth L. Judd, Lilia Maliar and Serguei Maliar, (2011)
 */
std::tuple<std::vector<std::vector<double>>, std::vector<double>> Monomials_1(const int nDim,
                                                                              const std::vector<double> &vcv){
    int n_nodes = 2*nDim;
    auto nDimD = static_cast<double>(nDim);

    std::vector<double> zl(n_nodes*nDim);

    for (int ii = 0; ii < nDim; ii++){
        zl[2*ii*nDim + ii] = 1.0;
        zl[(2*ii+1)*nDim + ii] = -1.0;
    }

    std::vector<double> RR = chol(nDim, vcv);
    int totSize = nDim*nDim;
    double tempval = sqrt(nDimD);

    for (int ii = 0; ii <totSize; ii++){
        RR[ii] = tempval*RR[ii];
    }

    std::vector<double> temp_epsi_nodes = matmul(n_nodes, nDim, zl, nDim, RR);
    std::vector<std::vector<double>> epsi_nodes = std::vector<std::vector<double>>(n_nodes, std::vector<double>(nDim, 0.0));

    for (int ii = 0; ii < n_nodes; ii++) {
        for (int jj = 0; jj < nDim; jj++) {
            epsi_nodes[ii][jj] = temp_epsi_nodes[ii * nDim + jj];
        }
    }

    std::vector<double> weight_nodes(n_nodes, 1.0 / static_cast<double>(n_nodes));

    return std::make_tuple(epsi_nodes, weight_nodes);
}

/*
 * Monomial rule version 2
 *
 * Algorithm adapted from
 * Solving Dynamic Economic Models" by Kenneth L. Judd, Lilia Maliar and Serguei Maliar, (2011)
 */
std::tuple<std::vector<std::vector<double>>, std::vector<double>> Monomials_2(const int nDim,
                                                                              const std::vector<double> &vcv){

    int n_nodes = 2*nDim*nDim + 1;
    auto nDimD = static_cast<double>(nDim);

    //std::vector<double> z0(nDim, 0.0);    // Exists in original codes, but not needed.
    std::vector<double> zl(2*nDim*nDim);

    for (int ii = 0; ii < nDim; ii++){
        zl[2*ii*nDim + ii] = 1.0;
        zl[(2*ii+1)*nDim + ii] = -1.0;
    }

    std::vector<double> z2(2*nDim*(nDim-1)*nDim, 0.0);

    int idx = 0;

    for (int pp = 0; pp < nDim-1; pp++){
        for (int qq = pp+ 1; qq < nDim; qq++){
            z2[ 4*idx*nDim      + pp] =  1.0;
            z2[(4*idx + 1)*nDim + pp] = -1.0;
            z2[(4*idx + 2)*nDim + pp] =  1.0;
            z2[(4*idx + 3)*nDim + pp] = -1.0;

            z2[ 4*idx*nDim      + qq] =  1.0;
            z2[(4*idx + 1)*nDim + qq] = -1.0;
            z2[(4*idx + 2)*nDim + qq] =  1.0;
            z2[(4*idx + 3)*nDim + qq] = -1.0;

            idx += 1;
        }
    }

    std::vector<double> RR = chol(nDim, vcv);
    std::vector<double> SS = RR;

    int totSize = nDim*nDim;
    double tempval1 = sqrt(nDimD+2.0);
    double tempval2 = sqrt((nDimD+2.0)/2.0);

    for (int ii = 0; ii <totSize; ii++){
        RR[ii] = tempval1*RR[ii];
        SS[ii] = tempval2*SS[ii];
    }

    std::vector<std::vector<double>> epsi_nodes = std::vector<std::vector<double>>(n_nodes, std::vector<double>(nDim, 0.0));
    std::vector<double> weight_nodes(n_nodes, 0.0);

    weight_nodes[0] = 2.0/(nDimD + 2.0);

    zl = matmul(2*nDim, nDim, zl, nDim, RR);
    z2 = matmul(2*nDim*(nDim-1), nDim, z2, nDim, SS);

    /*
     * Copy zl and weights for zl
     */
    int staidx = 1;
    int endidx = staidx + 2*nDim;
    for (int ii = staidx; ii < endidx; ii++){
        weight_nodes[ii] = ((4.0 - nDimD) / 2.0) / ((nDimD+2.0)*(nDimD+2.0));

        for (int jj = 0; jj < nDim; jj++) {
            epsi_nodes[ii][jj] = zl[(ii-staidx)*nDim + jj];
        }
    }

    staidx = endidx;
    endidx = staidx + 2*nDim*(nDim-1);
    for (int ii = staidx; ii < endidx; ii++) {
        weight_nodes[ii] = 1.0 / ((nDimD + 2.0) * (nDimD + 2.0));

        for (int jj = 0; jj < nDim; jj++){
            epsi_nodes[ii][jj] = z2[(ii-staidx)*nDim + jj];
        }
    }

    return std::make_tuple(epsi_nodes, weight_nodes);
}

/*
 * Gauss-Hermite Quadrature
 *
 * Algorithm adapted from
 * Solving Dynamic Economic Models" by Kenneth L. Judd, Lilia Maliar and Serguei Maliar, (2011)
 */
std::tuple<std::vector<std::vector<double>>, std::vector<double>> GH_Quadrature(const int Qn, const int nDim,
                                                                                const std::vector<double> &vcv) {
    std::vector<double> eps_vec(Qn, 0.0);
    std::vector<double> weight(Qn, 0.0);

    switch (Qn){
        case 1: {
            eps_vec = {0.0};
            weight = {sqrt(PI_VAL)};
            break;
        }
        case 2: {
            eps_vec = {0.7071067811865475, -0.7071067811865475};
            weight = {0.8862269254527580,  0.8862269254527580};
            break;
        }
        case 3: {
            eps_vec = {1.224744871391589, 0.0, -1.224744871391589};
            weight = {0.2954089751509193, 1.181635900603677, 0.2954089751509193};
            break;
        }
        case 4: {
            eps_vec = {1.650680123885785, 0.5246476232752903, -0.5246476232752903, -1.650680123885785};
            weight = {0.08131283544724518, 0.8049140900055128, 0.8049140900055128, 0.08131283544724518};
            break;
        }
        case 5: {
            eps_vec = {2.020182870456086, 0.9585724646138185, 0.0, -0.9585724646138185, -2.020182870456086};
            weight = {0.01995324205904591, 0.3936193231522412, 0.9453087204829419, 0.3936193231522412,
                      0.01995324205904591};
            break;
        }
        case 6: {
            eps_vec = {2.350604973674492, 1.335849074013697, 0.4360774119276165, -0.4360774119276165,
                       -1.335849074013697, -2.350604973674492};
            weight = {0.004530009905508846, 0.1570673203228566, 0.7246295952243925, 0.7246295952243925,
                      0.1570673203228566, 0.004530009905508846};
            break;
        }
        case 7: {
            eps_vec = {2.651961356835233, 1.673551628767471, 0.8162878828589647, 0.0,
                       -0.8162878828589647, -1.673551628767471, -2.651961356835233};
            weight = {0.0009717812450995192, 0.05451558281912703, 0.4256072526101278, 0.8102646175568073,
                      0.4256072526101278, 0.05451558281912703, 0.0009717812450995192};
            break;
        }
        case 8: {
            eps_vec = {2.930637420257244, 1.981656756695843, 1.157193712446780, 0.3811869902073221,
                       -0.3811869902073221, -1.157193712446780, -1.981656756695843, -2.930637420257244};
            weight = {0.0001996040722113676, 0.01707798300741348, 0.2078023258148919, 0.6611470125582413,
                      0.6611470125582413, 0.2078023258148919, 0.01707798300741348, 0.0001996040722113676};
            break;
        }
        case 9: {
            eps_vec = {3.190993201781528, 2.266580584531843, 1.468553289216668, 0.7235510187528376, 0.0,
                       -0.7235510187528376, -1.468553289216668, -2.266580584531843, -3.190993201781528};
            weight = {0.00003960697726326438, 0.004943624275536947, 0.08847452739437657, 0.4326515590025558,
                      0.7202352156060510, 0.4326515590025558, 0.08847452739437657, 0.004943624275536947,
                      0.00003960697726326438};
            break;
        }
        default: {
            eps_vec = {3.436159118837738, 2.532731674232790, 1.756683649299882, 1.036610829789514, 0.3429013272237046,
                       -0.3429013272237046, -1.036610829789514, -1.756683649299882, -2.532731674232790,
                       -3.436159118837738};
            weight = {7.640432855232621e-06, 0.001343645746781233, 0.03387439445548106, 0.2401386110823147,
                      0.6108626337353258, 0.6108626337353258, 0.2401386110823147, 0.03387439445548106,
                      0.001343645746781233, 7.640432855232621e-06};
        }
    }

    int n_nodes = 1;
    double normed_weight = 1.0;

    for (int ii = 0; ii < nDim; ii++){
        n_nodes *= Qn;
        normed_weight /= sqrt(PI_VAL);
    }

    std::vector<double> z1(n_nodes*nDim, sqrt(2.0));
    std::vector<double> weight_nodes(n_nodes, normed_weight);

    int mod_div = 1;

    for (int ii = 0; ii < nDim; ii++){
        for (int jj = 0; jj < n_nodes; jj++){
            int idx = (jj/mod_div) % Qn;

            z1[jj*nDim + ii] *= eps_vec[idx];
            weight_nodes[jj] *= weight[idx];
        }

        mod_div *= Qn;
    }

    std::vector<double> RR = chol(nDim, vcv);
    std::vector<double> temp_epsi_nodes = matmul(n_nodes, nDim, z1, nDim, RR);

    std::vector<std::vector<double>> epsi_nodes = std::vector<std::vector<double>>(n_nodes, std::vector<double>(nDim, 0.0));
    for (int ii = 0; ii < n_nodes; ii++) {
        for (int jj = 0; jj < nDim; jj++) {
            epsi_nodes[ii][jj] = temp_epsi_nodes[ii * nDim + jj];
        }
    }

    return std::make_tuple(epsi_nodes, weight_nodes);
}

#endif //CPP_EXP_NUMINTEGRATE_H
