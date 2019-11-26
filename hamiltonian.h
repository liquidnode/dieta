/*
This code is part of DIETA

Authored by Kirill Alpin
*/

#pragma once
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <bitset>
#include <sstream> 
#include <Eigen/Dense>
#include <stdint.h>
#include "parser.h"

using Eigen::VectorXd;
using Eigen::VectorXcd;
using Eigen::MatrixXcd;

#define KA(X) K(X, Tag<T>())
#define PA(X) P(X, Tag<T>())

template<typename T>
struct Tag{};

class Hamiltonian
{
public:
    //stores the expectation value of a single operator <a>
    VectorXcd single_operator_densities;
    //is true if single_operator_densities is set
    bool single_operator_densities_set = false;
    //stores particle densities <n>
    VectorXd densities;
    //stores hole densities (derived using the Nambu CPT Greensfunction)
    VectorXd h_densities;
    //stores the U(1) symmetry breaking field strength
    VectorXcd F_field;
    //stores the U(1) symmetry breaking field strength caused by mean fields
    VectorXcd meanF_field;
    //stores the chemical potential
    VectorXd mu_field;
    //bath site stuff
    //holds the bath sites chemical potential
    VectorXd bath_mu;
    //holds the bath sites hybridization with lattice sites
    MatrixXcd bath_hopping;
    
    //given the Hamiltonian, this function set the field variables
    void setFields();
    
    //is true if current densities are set. 
    bool current_set = false;
    //is true if hardcore commutation relations are used. is set by extra flags (see parser.cxx)
    bool hardcore = false;
    //is true if the user wants to set interactions to zero. is set by extra flags (see parser.cxx)
    bool no_interaction = false;
    //holds current densities <a^dag_i a_j>
    MatrixXcd current_densities;
    //holds current densities <a_i a^dag_j> derived from G CPT
    MatrixXcd h_current_densities;
    //holds anomalous current densities <a_i a_j>
    MatrixXcd anom_current_densities;
    //holds information about the entire system
    LatticeGraph graph;
    //is true if Hamiltonian has complex parameters
    bool is_complex;
    //is true if Hilbert space is truncated using particle number conservation. is set by extra flags (see parser.cxx)
    bool is_particle_conserving;
    //is true if Hilbert space is truncated using parity conservation. is set by extra flags (see parser.cxx)
    bool is_parity_conserving;
    //is true if the user wants to apply a particle hole transformation. is set by extra flags (see parser.cxx)
    bool ph_transform;
    //holds the Hilbert space dimension. if Hilbert space truncation is applied, dim is the dimension of the groundstate sector
    unsigned int dim;
    //holds dimension of the groundstate+particle sector. only used when Hilbert space truncation is applied
    unsigned int dim_part;
    //holds dimension of the groundstate+hole sector. only used when Hilbert space truncation is applied
    unsigned int dim_hole;
    //used to determine the current particle space the code is operating on
    unsigned int current_particle_space;
    //holds the number of sites
    unsigned int numsites;
    //holds the groundstate particle number. only used when Hilbert space truncation is applied
    unsigned int particle_number;
    //number of space dimensions of the lattice
    unsigned int space_dims;
    //precalulated sqrt(n) coefficients for boson operators
    FLTYPE* sqrt_coeff;
    //precalulated binomial coefficients for Hilbert space truncation code
    unsigned int* binom_coeff;
    //for particle number conserving hamiltonian
    unsigned int* occup_list;
    unsigned int* occup_list_part;
    unsigned int* occup_list_hole;
    //for parity conserving hamiltonian
    unsigned int** occup_lists;
    unsigned int* dims;
    
    unsigned int occup_to_index(unsigned int occup);
    unsigned int index_to_occup(unsigned int index);
    unsigned int get_occup(unsigned int occup, unsigned int site);
    FLTYPE n_operator(unsigned int occup, unsigned int site);
    FLTYPE holprim_expansion(int n);
    bool add_particle(unsigned int &occup, unsigned int site, FLTYPE& coeff);
    bool remove_particle(unsigned int &occup, unsigned int site, FLTYPE& coeff);
    void set_particle_space(unsigned int p);
    FARRAY single_particle_state(unsigned int site);
    std::complex<FLTYPE> P(int i, Tag<VectorXcd>);
    FLTYPE P(int i, Tag<VectorXd>);
    std::complex<FLTYPE> K(std::complex<FLTYPE> val, Tag<VectorXcd>);
    FLTYPE K(std::complex<FLTYPE> val, Tag<VectorXd>);
    template<typename T>
    T operator*(const T &s);
    template<typename T>
    T nop_site(T s, int site_num);
    template<typename T, typename U=FLTYPE, typename G>
    void add_particle_state(G s, int site, T &ts, U factor=U(1.0));
    template<typename T, typename U=FLTYPE, typename G>
    void remove_particle_state(G s, int site, T &ts, U factor=U(1.0));
    Hamiltonian(std::string graph_filename, std::string params_filename);
    Hamiltonian(std::string graph_filename, std::string params_filename, unsigned int particle_number);
    ~Hamiltonian();
private:
    unsigned int get_dimension();
    void compute_occup_lookup_table(unsigned int* &occup_lookup, unsigned int particle_number_now, unsigned int dim_now);
    void check_hilbert_dim(unsigned int dim_now);
};
