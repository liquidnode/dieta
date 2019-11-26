/*
This code is part of DIETA

Authored by Kirill Alpin
*/

#pragma once
#include <iostream>
#include <vector>
#include <random>
#include <Eigen/Dense>
#include "hamiltonian.h"

using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using Eigen::SelfAdjointEigenSolver;

template <typename T> 
struct LanczosParams
{
    std::vector<FLTYPE> bsqr;
    std::vector<T> a;
};

/*
 * this class hosts the Lanczos code for Exact Diagonalization of Hamiltonians
 */
template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
class Lanczos
{
private:
    std::vector<FLTYPE> bsqr;
    std::vector<VALTYPE> a;
    Hamiltonian* H;
public:
    FARRAY energies;
    MATTYPE evecs;
public:
    Lanczos(Hamiltonian* hamiltonian);
    MATTYPE get_tridiagonal(int seed);
    VECTYPE get_groundstate(int seed, MATTYPE tridiagonal, FLTYPE &genergy);
    template <typename T, typename U> 
    LanczosParams<T> get_params(U init_state, bool hole=false);
    void get_Qmatrix(VECTYPE groundstate, MATTYPE& T, MATTYPE& Q_elec, MATTYPE& Q_hole, bool hole, MATTYPE phi_elec, MATTYPE phi_hole);
};
