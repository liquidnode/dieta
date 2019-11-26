/*
This code is part of DIETA

Authored by Kirill Alpin
*/

#pragma once
#include <Eigen/Dense>
#include "parser.h"
#include "hamiltonian.h"
#include "greens.h"
#include "vmatrix.h"
#include "cpt.h"
#include "gauss_legendre/gauss_legendre.h"
#include "bzIntegrator.h"

using Eigen::VectorXcd;
using Eigen::VectorXd;
using Eigen::MatrixXcd;

//using senechal integration
//only difference in single site integration
#define SENECHAL

//using self energy periodization
//only difference in mean integration
//untested
//#define PERIODIZATION

struct InData
{
public:
    std::complex<double> start;
    std::complex<double> end;
    std::complex<double> omega;
    MatrixXcd G;
    CPT* cpt;
    bool real;
    int site;
    bool mean;
    bool onedim;
    bool nambu;
    bool is_quadri;
    bool with_hole_dens;
    bool inverse_omega = false;
    int poly_bz_start;
    
    InData()
    {
    }
    
    InData(InData* other)
    {
        //make deep copy
        start = other->start;
        end = other->end;
        omega = other->omega;
        G = MatrixXcd(other->G);
        cpt = other->cpt;
        real = other->real;
        site = other->site;
        mean = other->mean;
        onedim = other->onedim;
        nambu = other->nambu;
        is_quadri = other->is_quadri;
        with_hole_dens = other->with_hole_dens;
        inverse_omega = other->inverse_omega;
        poly_bz_start = other->poly_bz_start;
    }
};

/*
 * this class is used to calculate operator expectation values
 */
class Density
{
public:
    static double calcDensity(CPT* cpt, int site, bool mean=false, bool nambu=true);
    static VectorXd calcAllDensity(CPT* cpt, bool nambu=true, bool with_hole_dens=false);
    template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
    static VectorXd calcAllClusterDensity(Hamiltonian* H, VECTYPE& groundstate);
    static MatrixXcd calcEverything(CPT* cpt, bool nambu, bool with_hole_dens=false);
    template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
    static VECTYPE calcEverythingCluster(Hamiltonian* H);
    template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
    static void calcClusterSingleOperatorDensity(CPT* cpt);
};
