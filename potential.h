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
#include "density.h"
#include "gauss_legendre/gauss_legendre.h"
#include "bzIntegrator.h"

using Eigen::VectorXcd;
using Eigen::VectorXd;
using Eigen::MatrixXcd;
using Eigen::ComplexEigenSolver;

//does not work properly. reason unknown
//#define SENECHAL_

//integration using a finite contour
#define CONTOUR_

//use self energy periodization. experimental
//#define PERIODIZED

struct InDataPotential
{
public:
    std::complex<double> start;
    std::complex<double> end;
    std::complex<double> omega;
    MatrixXcd G;
    MatrixXcd Ginv;
    std::complex<double> detG;
    CPT* cpt;
    bool nambu = false;
    bool is_quadri = false;
    bool onedim = false;
    bool inverse_omega = false;
    bool real = true;
    int poly_bz_start;
    double zero_freq_contrib;
    bool cmplx_eigenval = false;
    double max_imag_eigenval = 0.0;
    std::vector<double> neg_pols;
    InDataPotential()
    {
    }
    InDataPotential(InDataPotential* cop)
    {
        //make deep copy
        start = cop->start;
        end = cop->end;
        omega = cop->omega;
        G = cop->G;
        Ginv = cop->Ginv;
        detG = cop->detG;
        cpt = cop->cpt;
        nambu = cop->nambu;
        is_quadri = cop->is_quadri;
        onedim = cop->onedim;
        inverse_omega = cop->inverse_omega;
        real = cop->real;
        poly_bz_start = cop->poly_bz_start;
        zero_freq_contrib = cop->zero_freq_contrib;
        cmplx_eigenval = cop->cmplx_eigenval;
        max_imag_eigenval = cop->max_imag_eigenval;
        neg_pols = std::vector<double>(cop->neg_pols);
    }
};

/*
 * this class is used to compute the Grand Potential using either "Bosonic Self-Energy Theory" 
 * for bosons/hardcore bosons or "Variational Cluster Approach" for fermions
 */
class Potential
{
public:
    static bool cmplx_eigenval;
    static double max_imag_eigenval;
    static double calcGrandPotential(CPT* cpt);
    template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
    static double calcClusterGrandPotential(Hamiltonian H);
    template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
    static double calcClusterFreeEnergy(Hamiltonian H);
};
