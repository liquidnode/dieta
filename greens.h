/*
This code is part of DIETA

Authored by Kirill Alpin
*/

#pragma once
#include <Eigen/Dense>
#include <vector>
#include "lanczos.h"
#include "hamiltonian.h"


using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;

class Dummy {   
public:
    virtual void prepare_greens(bool no_interaction=false) {};
    virtual MatrixXcd get_greens(std::complex<FLTYPE> omega) {};
    virtual MatrixXcd get_nambu_greens(std::complex<FLTYPE> omega) {};
    virtual MatrixXcd get_noninteracting_inverse_greens(std::complex<FLTYPE> omega, bool nambu, bool with_mean_field=true) {};
    virtual FARRAY clusterDensities() {};
    virtual std::vector<FLTYPE> get_negative_poles() {};
    virtual MatrixXcd omegaInfLimit() {};
    virtual void load_params(std::string file) {};
    virtual void save_params(std::string file) {};
    virtual FLTYPE get_low_energy() {};
    virtual FLTYPE get_high_energy() {};
    virtual FLTYPE get_groundstate_energy() {};
    virtual void set_eps(FLTYPE eps) {};
    virtual void ph_symmetry(bool ph_sym) {};
    virtual bool is_full_diag() {};
    virtual Hamiltonian* get_H() {};
    virtual MatrixXcd bfh_sign(MatrixXcd g) {};
    virtual MatrixXcd get_Q() {};
    virtual MatrixXcd get_Lambda() {};
    virtual MatrixXd get_Omega() {};
};

/*
 * this class hosts the code to construct and handle cluster Greens functions
 */
template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
class Greens : public Dummy
{
public:
    //TODO dieser falg macht noch nichts
    bool force_ph_symmetry;
    Hamiltonian* H;
    Lanczos<VALTYPE, VECTYPE, MATTYPE> lanczos;
    int deg = 1;
    VECTYPE groundstate;
    std::complex<FLTYPE> groundenergy;
    FLTYPE highenergy;
    FLTYPE eps;
    bool full_diag;
    MATTYPE Q_matrix_elec;
    MATTYPE Q_matrix_elec_h;
    MATTYPE T_elec;
    MATTYPE Q_matrix_hole;
    MATTYPE Q_matrix_hole_e;
    MATTYPE T_hole;
    
    VECTYPE Q_hardcore_proxy;
    VectorXd T_elec_energies;
    VectorXd T_hole_energies;
    
    MatrixXcd get_Q()
    {
        if(H->graph.is_nambu)
        {
            MatrixXcd ret(H->numsites * 2, Q_matrix_elec.cols() + Q_matrix_hole.cols());
            ret.block(0,0,H->numsites,Q_matrix_elec.cols()) << Q_matrix_elec;
            ret.block(0,Q_matrix_elec.cols(),H->numsites,Q_matrix_hole.cols()) << Q_matrix_hole.conjugate();
            if(full_diag)
            {
                ret.block(H->numsites,0,H->numsites,Q_matrix_elec.cols()) << Q_matrix_hole;
                ret.block(H->numsites,Q_matrix_elec.cols(),H->numsites,Q_matrix_hole.cols()) << Q_matrix_elec.conjugate();
            }
            else
            {
                ret.block(H->numsites,0,H->numsites,Q_matrix_elec.cols()) << Q_matrix_elec_h;
                ret.block(H->numsites,Q_matrix_elec.cols(),H->numsites,Q_matrix_hole.cols()) << Q_matrix_hole_e.conjugate();
            }
            return ret;
        }
        else
        {
            MatrixXcd ret(H->numsites, Q_matrix_elec.cols() + Q_matrix_hole.cols());
            ret.block(0,0,H->numsites,Q_matrix_elec.cols()) << Q_matrix_elec;
            ret.block(0,Q_matrix_elec.cols(),H->numsites,Q_matrix_hole.cols()) << Q_matrix_hole.conjugate();
            return ret;
        }
    }
    MatrixXcd get_Lambda()
    {
        if(full_diag)
        {
            VectorXcd c(lanczos.energies.size()*2);
            c << VectorXcd(lanczos.energies) - groundenergy * VectorXcd::Ones(lanczos.energies.size()), 
                groundenergy * VectorXcd::Ones(lanczos.energies.size()) - VectorXcd(lanczos.energies);
            return c.asDiagonal();
        }
        else
        {
            VectorXcd c(T_elec_energies.size() + T_hole_energies.size());
            c << VectorXcd(T_elec_energies) - groundenergy * VectorXcd::Ones(T_elec_energies.size()), 
                groundenergy * VectorXcd::Ones(T_hole_energies.size()) - VectorXcd(T_hole_energies);
            return c.asDiagonal();
        }
    }
    MatrixXd get_Omega()
    {
        if(H->graph.particleType == ParticleType::FERMION)
        {
            if(full_diag)
            {
                VectorXd c(lanczos.energies.size()*2);
                c << VectorXd::Ones(lanczos.energies.size()), VectorXd::Ones(lanczos.energies.size());
                return c.asDiagonal();
            }
            else
            {
                VectorXd c(T_elec_energies.size() + T_hole_energies.size());
                c << VectorXd::Ones(T_elec_energies.size()), VectorXd::Ones(T_hole_energies.size());
                return c.asDiagonal();
            }
        }
        else
        {
            if(full_diag)
            {
                VectorXd c(lanczos.energies.size()*2);
                c << VectorXd::Ones(lanczos.energies.size()), -VectorXd::Ones(lanczos.energies.size());
                return c.asDiagonal();
            }
            else
            {
                VectorXd c(T_elec_energies.size() + T_hole_energies.size());
                c << VectorXd::Ones(T_elec_energies.size()), -VectorXd::Ones(T_hole_energies.size());
                return c.asDiagonal();
            }
        }
    }
    std::vector<FLTYPE> get_negative_poles()
    {
        //TODO filter Q and eigs so that the ifs below are not necessary, ie remove eigenvalues
        std::vector<FLTYPE> ret;
        if(full_diag)
        {
            for(int i = 0; i < lanczos.energies.size(); ++i)
            {
                if(std::real(lanczos.energies[i] - groundenergy) < 0.0)
                //if((Q_matrix_elec.col(i)*Q_matrix_elec.col(i).adjoint()).cwiseAbs().maxCoeff() > 1e-10 || 
                //    (H->graph.is_nambu && (Q_matrix_elec.col(i)*Q_matrix_elec_h.col(i).adjoint()).cwiseAbs().maxCoeff() > 1e-10))
                    ret.push_back(std::real(lanczos.energies[i] - groundenergy));
                //if(std::real(lanczos.energies[i] - groundenergy) > 0.0 && H->graph.is_nambu)
                //    ret.push_back(-std::real(lanczos.energies[i] - groundenergy));
            }
            for(int i = 0; i < lanczos.energies.size(); ++i)
            {
                if(std::real(groundenergy - lanczos.energies[i]) < 0.0)
                //if((Q_matrix_hole.col(i)*Q_matrix_hole.col(i).adjoint()).cwiseAbs().maxCoeff() > 1e-10 || 
                //    (H->graph.is_nambu && (Q_matrix_hole.col(i)*Q_matrix_hole_e.col(i).adjoint()).cwiseAbs().maxCoeff() > 1e-10))
                    ret.push_back(std::real(groundenergy - lanczos.energies[i]));
                //if(std::real(groundenergy - lanczos.energies[i]) > 0.0 && H->graph.is_nambu)
                //    ret.push_back(-std::real(groundenergy - lanczos.energies[i]));
            }
        }
        else
        {
            for(int i = 0; i < T_elec_energies.size(); ++i)
            {
                if(std::real(T_elec_energies[i] - groundenergy) < 0.0)
                //if((Q_matrix_elec.col(i)*Q_matrix_elec.col(i).adjoint()).cwiseAbs().maxCoeff() > 1e-10 || 
                //    (H->graph.is_nambu && (Q_matrix_elec.col(i)*Q_matrix_hole.col(i).adjoint()).cwiseAbs().maxCoeff() > 1e-10))
                    ret.push_back(std::real(T_elec_energies[i] - groundenergy));
                if(std::real(T_elec_energies[i] - groundenergy) > 0.0)
                    ret.push_back(-std::real(T_elec_energies[i] - groundenergy));
            }
            for(int i = 0; i < T_hole_energies.size(); ++i)
            {
                if(std::real(groundenergy - T_hole_energies[i]) < 0.0)
                //if((Q_matrix_hole.col(i)*Q_matrix_hole.col(i).adjoint()).cwiseAbs().maxCoeff() > 1e-10 || 
                //    (H->graph.is_nambu && (Q_matrix_hole.col(i)*Q_matrix_elec.col(i).adjoint()).cwiseAbs().maxCoeff() > 1e-10))
                    ret.push_back(std::real(groundenergy - T_hole_energies[i]));
                if(std::real(groundenergy - T_hole_energies[i]) > 0.0)
                    ret.push_back(-std::real(groundenergy - T_hole_energies[i]));
            }
        }
        return ret;
    }
    FLTYPE get_low_energy()
    {
        return std::abs(std::real(groundenergy));
    }
    FLTYPE get_groundstate_energy()
    {
        return std::real(groundenergy);
    }
    FLTYPE get_high_energy()
    {
        return std::abs(highenergy);
    }
    void set_eps(FLTYPE eps_p)
    {
        eps = eps_p;
        groundenergy = std::complex<FLTYPE>(std::real(groundenergy),eps);
    }
    bool is_full_diag()
    {
        return full_diag;
    }
    Hamiltonian* get_H()
    {
        return H;
    }
    void ph_symmetry(bool ph_sym)
    {
        force_ph_symmetry = ph_sym;
    }
    void cleanMatrix(MATTYPE& x)
    {
        for(int n = 0; n < x.cols(); ++n)
            for(int m = 0; m < x.rows(); ++m)
                if(std::abs(x(m,n)) < 1e-11)
                    x(m,n) = 0.0;
    }
    std::vector<LanczosParams<VALTYPE>> gparams_elec;
    std::vector<LanczosParams<VALTYPE>> gparams_hole;
    
    std::vector<LanczosParams<VALTYPE>> gparams_P;
    std::vector<LanczosParams<std::complex<FLTYPE>>> gparams_I;
    
    Greens(Hamiltonian* hamiltonian, bool precalc=true);
    Greens(Hamiltonian* hamiltonian, std::string file);
    void prepare_greens(bool no_interaction=false);
    void load_params(std::string file);
    void save_params(std::string file);
    void load_vector(std::ifstream& infile, FARRAY& vec);
    void load_vector(std::ifstream& infile, VectorXcd& vec);
    void load_matrix(std::ifstream& infile, MATTYPE& matrix);
    std::string next_line(std::ifstream& infile);
    void parse_vector(std::ofstream& ofile, FARRAY vec);
    void parse_vector(std::ofstream& ofile, VectorXcd vec);
    void parse_matrix(std::ofstream& ofile, MATTYPE matrix);
    template <typename T> 
    std::complex<FLTYPE> calc_fraction(std::complex<FLTYPE> omega, LanczosParams<T> param, bool hole);
    //template <typename T>
    MatrixXcd bfh_sign(MatrixXcd g);
    MatrixXcd get_greens(std::complex<FLTYPE> omega);
    MatrixXcd get_anom_greens(std::complex<FLTYPE> omega, bool conjugate=false);
    MatrixXcd get_nambu_greens(std::complex<FLTYPE> omega);
    MatrixXcd get_noninteracting_inverse_greens(std::complex<FLTYPE> omega, bool nambu, bool with_mean_field=true);
    MatrixXcd get_bath_gamma(std::complex<FLTYPE> omega);
    
    FARRAY clusterDensities();
    MatrixXcd omegaInfLimit();
};
