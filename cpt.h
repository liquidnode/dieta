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

using Eigen::MatrixXcd;
using Eigen::Vector3d;

/*
 * The main "Cluster Perturbation Theory" class
 */
class CPT
{
public:
    Hamiltonian* H;
    Dummy* greens;
    VMatrix* vmatrix;
    CPT(std::string graph_filename, std::string params_filename, std::string lanczos_params_filename, bool load, bool precalc=true);
    ~CPT();
    std::complex<FLTYPE> G(std::complex<FLTYPE> omega, Vector3d k);
    std::complex<FLTYPE> F(std::complex<FLTYPE> omega, Vector3d k);
    std::complex<FLTYPE> nambu_G(std::complex<FLTYPE> omega, Vector3d k);
    MatrixXcd cluster_nambu_G(std::complex<FLTYPE> omega, Vector3d k);
    MatrixXcd cluster_nambu_G_precalc(std::complex<FLTYPE> omega, Vector3d k, MatrixXcd G);
    MatrixXcd cluster_G(std::complex<FLTYPE> omega, Vector3d k);
    MatrixXcd cluster_G_precalc(std::complex<FLTYPE> omega, Vector3d k, MatrixXcd G);
    MatrixXcd cluster_noninteracting_G(std::complex<FLTYPE> omega, Vector3d k);
    MatrixXcd cluster_noninteracting_inverse_G(std::complex<FLTYPE> omega, Vector3d k);
    std::vector<std::complex<FLTYPE>> nambu_G_k_sweep(std::complex<FLTYPE> omega, Vector3d k_start, Vector3d k_end, int k_num_steps, bool with_self_energy_periodization = false);
    std::vector<std::complex<FLTYPE>> G_k_sweep(std::complex<FLTYPE> omega, Vector3d k_start, Vector3d k_end, int k_num_steps, bool with_self_energy_periodization = false);
    std::vector<std::complex<FLTYPE>> noninteracting_G_k_sweep(std::complex<FLTYPE> omega, Vector3d k_start, Vector3d k_end, int k_num_steps);
    void recalculate_greens(bool no_interaction = false);
    std::complex<FLTYPE> self_energy(std::complex<FLTYPE> omega, Vector3d k);
    std::complex<FLTYPE> G_with_self_energy(std::complex<FLTYPE> omega, Vector3d k);
    std::complex<FLTYPE> self_energy_precalc(std::complex<FLTYPE> omega, Vector3d k, MatrixXcd Ginv);
    std::complex<FLTYPE> G_with_self_energy_precalc(std::complex<FLTYPE> omega, Vector3d k, MatrixXcd Ginv);
    MatrixXcd partial_ft_self_energy_precalc(std::complex<FLTYPE> omega, Vector3d k, MatrixXcd Ginv);
    MatrixXcd partial_ft_G_with_self_energy_precalc(std::complex<FLTYPE> omega, Vector3d k, MatrixXcd Ginv);
};
