/*
This code is part of DIETA

Authored by Kirill Alpin
*/

#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "./exprtk/exprtk.hpp"


using Eigen::VectorXd;
using Eigen::VectorXcd;
using Eigen::Vector3d;
using FLTYPE = double;
using FARRAY = VectorXd;
/*
 * for an automatic Holstein Primakoff expansion
 * only usable with hardcore bosons
 */
constexpr int MAX_HOLPRIM_ORDER = 0;
constexpr FLTYPE HOLPRIM_COEFF[10] = {1.0/2.0, 1.0/8.0, 1.0/16.0, 5.0/128.0, 7.0/256.0, 21.0/1024.0, 33.0/2048.0, 429.0/32768.0, 715.0/65536.0, 2431.0/262144.0};

/* FOREGROUND */
#define RST  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

#define FRED(x) KRED x RST
#define FGRN(x) KGRN x RST
#define FYEL(x) KYEL x RST
#define FBLU(x) KBLU x RST
#define FMAG(x) KMAG x RST
#define FCYN(x) KCYN x RST
#define FWHT(x) KWHT x RST

#define BOLD(x) "\x1B[1m" x RST
#define UNDL(x) "\x1B[4m" x RST

struct Vec
{
    Vector3d pos;
    char type;
    int num;
    int ind;
    bool is_edge = false;
};

struct Edge
{
    unsigned int a;
    unsigned int b;
    char type;
    int trans_index;
    int ind;
};

enum class GPIntegrationMode {AUTOMATIC, EXACT, CONTOUR};
enum class ParticleType {BOSON, FERMION, HARDCORE};

/*
 * LatticeGraph holds information of the entire system and also extra flags
 */
class LatticeGraph
{
public:
    //for a description of these variables see parser.cxx
    std::vector<std::string> param_names;
    std::vector<std::complex<FLTYPE>> params;
    std::vector<Vec> site_pos;
    std::vector<Edge> edges;
    std::vector<Edge> inter_edges;
    std::vector<Edge> mean_edges;
    std::vector<Vec> translations;
    std::vector<Vector3d> sym_points;
    std::vector<Vector3d> poly_bz;
    std::vector<FLTYPE> poly_weight;
    double eps;
    int k_steps;
    bool hardcore_current;
    bool is_complex;
    bool is_particle_conserving;
    bool is_half_filling;
    bool is_parity_conserving;
    bool is_nambu;
    bool is_onedim;
    bool cluster_dens;
    bool better_densities;
    bool ph_transform;
    bool use_current_decoupling;
    bool F_field_zero_in_continuum_limit;
    bool real_F_field_is_mean_field;
    ParticleType particleType;
    GPIntegrationMode intMode;
    int max_particle;
    unsigned int particle_number;
    unsigned int num_bath;
    double omega_start;
    double omega_end;
    int omega_steps;
    FLTYPE parse_float(std::string expression_string);
    LatticeGraph(std::string graph_filename, std::string params_filename);
}; 
