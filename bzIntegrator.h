/*
This code is part of DIETA

Authored by Kirill Alpin
*/

#pragma once
#include <Eigen/Dense>
#include "parser.h"
#include "cpt.h"
#include <complex>
#include "gauss_legendre/gauss_legendre.h"

using Eigen::VectorXd;
using Eigen::VectorXcd;
using Eigen::MatrixXd;
using Eigen::MatrixXcd;

template<class ITYPE>
class BZData
{
public:
    void* data;
    CPT* cpt;
    ITYPE (*f)(double,double,void*);
    int poly_bz_start;
};


/*
 * this class hosts methods, which can integrate an arbitrary function over the Brillouin Zone
 */
class BZIntegrator
{
public:
    static double nodal(double eta, double nu, unsigned int i)
    {
        return 0.25 * (1.0 - (1.0 - 2.0 * (i & 1)) * eta) * (1.0 - (1.0 - 2.0 * ((i & 2) >> 1)) * nu);
    };

    static double nodal_deta(double eta, double nu, unsigned int i)
    {
        return 0.25 * (- (1.0 - 2.0 * (i & 1))) * (1.0 - (1.0 - 2.0 * ((i & 2) >> 1)) * nu);
    };

    static double nodal_dnu(double eta, double nu, unsigned int i)
    {
        return 0.25 * (1.0 - (1.0 - 2.0 * (i & 1)) * eta) * (- (1.0 - 2.0 * ((i & 2) >> 1)));
    };

    static double jacobian(double eta, double nu, std::vector<Vector3d> poly_bz, int poly_bz_start)
    {
        double x_deta = 0.0;
        double y_deta = 0.0;
        double x_dnu = 0.0;
        double y_dnu = 0.0;
        for(int i = 0; i < 4; ++i)
        {
            x_deta += poly_bz[i + poly_bz_start][0]*nodal_deta(eta, nu, i);
            y_deta += poly_bz[i + poly_bz_start][1]*nodal_deta(eta, nu, i);
            x_dnu += poly_bz[i + poly_bz_start][0]*nodal_dnu(eta, nu, i);
            y_dnu += poly_bz[i + poly_bz_start][1]*nodal_dnu(eta, nu, i);
        }
        return x_deta * y_dnu - x_dnu * y_deta;
    };
    
    template<class ITYPE>
    static ITYPE bzIntegrate(CPT* cpt, ITYPE (*f)(double,double,void*), void* data, ITYPE init, bool parallel);
};
