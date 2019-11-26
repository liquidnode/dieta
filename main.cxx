/*
 * DIETA (<D>IETA <I>st <E>in <T>olles <A>kronym) was created by Kirill Alpin. 
 * This code is able to calculate excitation spectra of spin/bosonic/fermionic systems using 
 * Cluster Perturbation Theory (CPT) and also implements Bosonic Self-Energy Functional Theory(BSFT)/
 * Variational Cluster Approach (VCA). See the examples folder for a few use cases of this program.
 * 
 * 
 * 
 * 
 */


#include <iostream>
#include <iomanip>
#include <fstream>
#include <bitset>
#include <Eigen/Dense>
#include "parser.h"
#include "hamiltonian.h"
#include "lanczos.h"
#include "greens.h"
#include "cpt.h"
#include "meanfield.h"
#include "potential.h"
#include "gauss_legendre/gauss_legendre.h"

using Eigen::Vector3d;

int main(int argc, char* argv[])
{
    const std::complex<FLTYPE> img(0, 1);
    if(argc < 3)
    {
        std::cerr << "Two arguments needed." << std::endl;
        return 1;
    }
    //flags are described in the code below
    bool load = false;
    bool calcGrandPotential = false;
    bool calcParticleDensity = false;
    bool calcMuFields = false;
    bool cluster = false;
    bool calcFreeEnergy = false;
    bool calcSingleOperatorDensity = false;
    bool withMeanField = false;
    bool periodize_self_energy = false;
    bool plot_self_energy = false;
    bool get_G_at_k_omega = false;
    bool self_energy_period = false;
    Vector3d sk;
    double somega;
    bool witheps = false;
    std::string lanczos_params_file = "";
    std::string lanczos_params_file_tmp = "";
    bool spec_output = false;
    std::string spec_output_file = "";
    //iterate over all flags
    for(int fargc = 3; fargc < argc; ++fargc)
    {
        std::string flag = argv[fargc];
        if(flag == "-l" && argc>fargc+1)
        {
            /*
             * -l <file_path> is used to load a Greens function stored in <file_path>
             */
            load = true;
            lanczos_params_file = argv[fargc+1];
            fargc++;
        }
        if(flag == "-s" && argc>fargc+1)
        {
            /*
             * -s <file_path> is used to save the calculated Greens function to <file_path>
             */
            load = false;
            lanczos_params_file = argv[fargc+1];
            fargc++;
        }
        if(flag == "-gp")
        {
            /*
             * -gp calculates and returns the grand potential of the system
             * (if -c is set, the cluster grand potential is returned)
             */
            calcGrandPotential = true;
        }
        if(flag == "-pd")
        {
            /*
             * -pd calculates and returns the particle density of the system 
             * using either the CPT Greens function or the cluster groundstate
             * if the flag -c is also set
             */
            calcParticleDensity = true;
        }
        if(flag == "-mf")
        {
            /*
             * -mf performes a mean field decoupling
             */
            withMeanField = true;
        }
        if(flag == "-mu")
        {
            /*
             * -mf returns the chemical potential
             */
            calcMuFields = true;
        }
        if(flag == "-sp")
        {
            /*
             * -sp returns the expectation value <a>
             */
            calcSingleOperatorDensity = true;
        }
        if(flag == "-c")
        {
            /*
             * -c specifies if cluster values are returned
             */
            cluster = true;
        }
        if(flag == "-cgp")
        {
            /*
             * is equivalent to -c -gp
             */
            calcGrandPotential = true;
            cluster = true;
        }
        if(flag == "-cpd")
        {
            /*
             * is equivalent to -c -pd
             */
            calcParticleDensity = true;
            cluster = true;
        }
        if(flag == "-cfe")
        {
            /*
             * returnes the cluster free energy
             */
            calcFreeEnergy = true;
            cluster = true;
        }
        if(flag == "-csp")
        {
            /*
             * is equivalent to -c -sp
             */
            calcSingleOperatorDensity = true;
            cluster = true;
        }
        if(flag == "-sf")
        {
            /*
             * instead of the Greens function spectrum, the self energy one is calculated
             */
            plot_self_energy = true;
        }
        if(flag == "-G" && argc>fargc+4)
        {
            /*
             * -G returns the Greens function at the specified frequency and wave vector
             */
            get_G_at_k_omega = true;
            somega = std::stod(argv[fargc+1]);
            sk[0] = std::stod(argv[fargc+2]);
            sk[1] = std::stod(argv[fargc+3]);
            sk[2] = std::stod(argv[fargc+4]);
            fargc += 4;
        }
        if(flag == "-Geps" && argc>fargc+4)
        {
            /*
             * -Geps returns the Greens function at the specified frequency and wave vector
             * with artificial peak broadening
             */
            witheps = true;
            get_G_at_k_omega = true;
            somega = std::stod(argv[fargc+1]);
            sk[0] = std::stod(argv[fargc+2]);
            sk[1] = std::stod(argv[fargc+3]);
            sk[2] = std::stod(argv[fargc+4]);
            fargc += 4;
        }
        if(flag == "-sfp")
        {
            /*
             * -sfp specifies the use of self-energy periodization when spectra are calculated
             */
            self_energy_period = true;
        }
        if(flag == "-o" && argc>fargc+1)
        {
            /*
             * -o specifies the file name/path the spectrum file is created
             */
            spec_output_file = argv[fargc+1];
            spec_output = true;
            fargc++;
        }
    }
    if(lanczos_params_file != "" && !load)
    {
        lanczos_params_file_tmp = lanczos_params_file;
        lanczos_params_file = "";
    }
    //create an instance of CPT using the graph and parameter file. 
    //load the Greens function from <lanczos_params_file> if required
    CPT cpt(argv[1], argv[2], lanczos_params_file, load, false);
    if(cpt.H->graph.eps == 0.0)
    {
        std::cerr << "Eps is 0." << std::endl;
        return 1;
    }
    //perform a mean field decoupling if required
    if(cpt.H->graph.mean_edges.size() != 0 && withMeanField)
    {
        MeanField meanf(&cpt);
        if(cpt.H->graph.use_current_decoupling)
            meanf.findFixpointWithCurrent();
        else
            meanf.findFixpoint();
    }
    //return the free energy
    if(calcFreeEnergy)
    {
        double freeEnergy = 0.0;
        if(cluster)
            if(cpt.H->is_complex)
                freeEnergy = Potential::calcClusterFreeEnergy<std::complex<FLTYPE>, VectorXcd, MatrixXcd>(*cpt.H);
            else
                freeEnergy = Potential::calcClusterFreeEnergy<FLTYPE, FARRAY, MatrixXd>(*cpt.H);
        else
        {
            std::cerr << "Non cluster free energy not implemented." << std::endl;
        }
        std::cout << "FE\t" << std::setprecision(16) << freeEnergy << std::endl;
        return 0;
    }
    //return the grand potential
    if(calcGrandPotential)
    {
        double grandPotential = 0.0;
        if(cluster)
            if(cpt.H->is_complex)
                grandPotential = Potential::calcClusterGrandPotential<std::complex<FLTYPE>, VectorXcd, MatrixXcd>(*cpt.H);
            else
                grandPotential = Potential::calcClusterGrandPotential<FLTYPE, FARRAY, MatrixXd>(*cpt.H);
        else
            grandPotential = Potential::calcGrandPotential(&cpt);
        if(Potential::cmplx_eigenval)
            std::cout << "Complex eigenvalue encountered with imag:\t" << Potential::max_imag_eigenval << std::endl;
        std::cout << "GP\t" << std::setprecision(16) << grandPotential << std::endl;
        return 0;
    }
    //return the single operator expectation value
    if(calcSingleOperatorDensity)
    {
        //TODO not efficient // just calc groundstate only
        cpt.recalculate_greens();
        
        if(cpt.H->is_complex)
            Density::calcClusterSingleOperatorDensity<std::complex<FLTYPE>, VectorXcd, MatrixXcd>(&cpt);
        else
            Density::calcClusterSingleOperatorDensity<FLTYPE, FARRAY, MatrixXd>(&cpt);
        cpt.H->single_operator_densities_set = false;
        
        if(!cluster)
        {
            if(!cpt.H->graph.F_field_zero_in_continuum_limit)
            {
                cpt.H->setFields();
                VectorXcd F_field = cpt.H->F_field;
                
                VectorXcd real_F_field = cpt.H->single_operator_densities;
                if(cpt.H->graph.real_F_field_is_mean_field)
                    real_F_field = cpt.H->meanF_field;
                
                MatrixXcd G_0prime_inv = cpt.greens->get_noninteracting_inverse_greens(0, cpt.greens->get_H()->graph.is_nambu);
                MatrixXcd G_0_inv = cpt.cluster_noninteracting_inverse_G(0, Vector3d(0,0,0));
                if(cpt.greens->get_H()->graph.is_nambu)
                {
                    VectorXcd tmp = cpt.H->single_operator_densities;
                    cpt.H->single_operator_densities = VectorXcd(cpt.H->numsites * 2);
                    cpt.H->single_operator_densities << tmp, tmp.conjugate();
                    tmp = F_field;
                    F_field = VectorXcd(cpt.H->numsites * 2);
                    F_field << tmp, tmp.conjugate();
                    tmp = real_F_field;
                    real_F_field = VectorXcd(cpt.H->numsites * 2);
                    real_F_field << tmp, tmp.conjugate();
                }
                VectorXcd FminusSigma = G_0prime_inv * cpt.H->single_operator_densities;
                VectorXcd Sigma = F_field - FminusSigma;
                cpt.H->single_operator_densities = G_0_inv.inverse() * (real_F_field - Sigma);
            }
            else
            {
                cpt.H->setFields();
                VectorXcd F_field = cpt.H->F_field;
                
                MatrixXcd G_0prime_inv = cpt.greens->get_noninteracting_inverse_greens(0, cpt.greens->get_H()->graph.is_nambu);
                MatrixXcd G_0_inv = cpt.cluster_noninteracting_inverse_G(0, Vector3d(0,0,0));
                
                if(cpt.greens->get_H()->graph.is_nambu)
                {
                    VectorXcd tmp = cpt.H->single_operator_densities;
                    cpt.H->single_operator_densities = VectorXcd(cpt.H->numsites * 2);
                    cpt.H->single_operator_densities << tmp, tmp.conjugate();
                    tmp = F_field;
                    F_field = VectorXcd(cpt.H->numsites * 2);
                    F_field << tmp, tmp.conjugate();
                }
                VectorXcd FminusSigma = G_0prime_inv * cpt.H->single_operator_densities;
                VectorXcd Sigma = F_field - FminusSigma;
                cpt.H->single_operator_densities = -G_0_inv.inverse() * Sigma;
            }
        }
        
        std::cout << "SOD" << std::setprecision(16) << std::endl;
        if(cpt.H->is_complex)
            std::cout << "C";
        for(int i = 0; i < cpt.H->numsites; ++i)
        {
            if(cpt.H->is_complex)
                std::cout << std::real(cpt.H->single_operator_densities[i]) << "\t" << std::imag(cpt.H->single_operator_densities[i]);
            else
                std::cout << std::real(cpt.H->single_operator_densities[i]);
            if(i != cpt.H->numsites - 1)
                std::cout << "\t";
            else
                std::cout << "\n";
        }
        return 0;
    }
    //return the particle density
    if(calcParticleDensity)
    {
        if(cluster)
        {
            VectorXd density = VectorXd::Zero(cpt.H->numsites);
            if(cpt.H->is_complex)
                density = Density::calcAllClusterDensity<std::complex<FLTYPE>, VectorXcd, MatrixXcd>(cpt.H, ((Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>*)cpt.greens)->groundstate);
            else
                density = Density::calcAllClusterDensity<FLTYPE, FARRAY, MatrixXd>(cpt.H, ((Greens<FLTYPE, FARRAY, MatrixXd>*)cpt.greens)->groundstate);
            
            std::cout << "PD" << std::setprecision(16) << std::endl;
            for(int i = 0; i < cpt.H->numsites; ++i)
            {
                std::cout << density[i];
                if(i != cpt.H->numsites - 1)
                    std::cout << "\t";
                else
                    std::cout << "\n";
            }
            return 0;
        }
        
        cpt.greens->set_eps(0.0);
        cpt.recalculate_greens();
        
        if(periodize_self_energy)
        {
            double dens = Density::calcDensity(&cpt, 0, true, cpt.H->graph.is_nambu);
            std::cout << "PD" << std::setprecision(16) << std::endl;
            std::cout << dens << std::endl;
            return 0;
        }
        
        VectorXd density;
        if(cpt.H->graph.better_densities)
        {
            if(cpt.H->graph.max_particle != 1)
            {
                std::stringstream ss;
                ss << "Corrected densities are only possible when max_particle is equal to 1.";
                throw std::invalid_argument(ss.str());
            }
            VectorXd alldensities = Density::calcAllDensity(&cpt, cpt.H->graph.is_nambu, true);
            density = alldensities.head(cpt.H->numsites);
            VectorXd ph_densities = alldensities.tail(cpt.H->numsites);
            std::cout << "Particle densities: " << std::endl << density << std::endl << std::endl;
            std::cout << "Hole densities: " << std::endl << ph_densities << std::endl;
            density = (density + (VectorXd::Ones(cpt.H->numsites) - ph_densities)) / 2.0;
        }
        else
            density = Density::calcAllDensity(&cpt, cpt.H->graph.is_nambu);
        std::cout << "PD" << std::setprecision(16) << std::endl;
        for(int i = 0; i < cpt.H->numsites; ++i)
        {
            std::cout << density[i];
            if(i != cpt.H->numsites - 1)
                std::cout << "\t";
            else
                std::cout << "\n";
        }
        return 0;
    }
    //return the chemical potential
    if(calcMuFields)
    {
        cpt.H->setFields();
        std::cout << "MU\t" << FRED("ACHTUNG!") << " Per Definition H=H_0" << FRED("-") << "mu*N" << std::setprecision(16) << std::endl;
        for(int i = 0; i < cpt.H->numsites; ++i)
        {
            std::cout << cpt.H->mu_field[i];
            if(i != cpt.H->numsites - 1)
                std::cout << "\t";
            else
                std::cout << "\n";
        }
        return 0;
    }
    FLTYPE epstmp = cpt.H->graph.eps;
    
    //prepare the Greens function
    cpt.recalculate_greens();
    
    if(lanczos_params_file_tmp != "")
        cpt.greens->save_params(lanczos_params_file_tmp);
    cpt.greens->set_eps(0.0);
    
    //if only one point of the Greens function is needed, return said point
    if(get_G_at_k_omega)
    {
        if(!witheps)
            cpt.H->graph.eps = 0.0;
        std::complex<double> ret = 0;
        if(self_energy_period)
        {
            //use self energy periodization when neede
            MatrixXcd G;
            if(cpt.H->graph.is_nambu)
                G = cpt.greens->get_nambu_greens(somega+cpt.H->graph.eps*img);
            else
                G = cpt.greens->get_greens(somega+cpt.H->graph.eps*img);
            G = G.inverse();
            ret = cpt.partial_ft_G_with_self_energy_precalc(somega+cpt.H->graph.eps*img, sk, G)(0,0)+
                cpt.partial_ft_G_with_self_energy_precalc(somega+cpt.H->graph.eps*img, sk, G)(1,0)+
                cpt.partial_ft_G_with_self_energy_precalc(somega+cpt.H->graph.eps*img, sk, G)(0,1)+
                cpt.partial_ft_G_with_self_energy_precalc(somega+cpt.H->graph.eps*img, sk, G)(1,1);
        }
        else
        {
            if(cpt.H->graph.is_nambu)
                ret = cpt.nambu_G(somega+cpt.H->graph.eps*img, sk);
            else
                ret =  cpt.G(somega+cpt.H->graph.eps*img, sk);
        }
        std::cout << "G" << std::setprecision(16) << "\t";
        std::cout << std::real(ret) << "\t" << std::imag(ret) << std::endl;
        return 0;
    }
    
    //if no physical observables are returned, calculate the exciation spectrum
    const int num_sym = cpt.H->graph.sym_points.size();
    
    if(!spec_output)
    {
        spec_output_file = "cpt_res/spec_output.txt";
    }
    //open spectral output file
    std::ofstream results(spec_output_file);
    //set up bounds
    FLTYPE omega_start = cpt.H->graph.omega_start;
    FLTYPE omega_end = cpt.H->graph.omega_end;
    int max_omega_step = cpt.H->graph.omega_steps;
    FLTYPE eps = cpt.H->graph.eps;
    int k_num_steps = cpt.H->graph.k_steps;
    results << max_omega_step << "\t" << k_num_steps << "\t" << num_sym << std::endl;
    //perform a sweep over frequencies omega
    for(int p = 0; p < max_omega_step; ++p)
    {
        FLTYPE omega = (omega_end * (FLTYPE)p + omega_start * (FLTYPE)(max_omega_step - 1 - p)) / (FLTYPE)(max_omega_step - 1);
        for(int i = 0; i < num_sym - 1; ++i)
        {
            std::vector<std::complex<FLTYPE>> b;
            if(plot_self_energy)
            {
                //calculate the self-energy instead of the Greens function
                MatrixXcd Ginv;
                if(cpt.H->graph.is_nambu)
                   Ginv = cpt.greens->get_nambu_greens(omega);
                else
                   Ginv = cpt.greens->get_greens(omega);
                for(int j = 0; j < k_num_steps; ++j)
                {
                    Vector3d k = (cpt.H->graph.sym_points[(i+1)%num_sym] * (FLTYPE)j + cpt.H->graph.sym_points[i]*((FLTYPE)k_num_steps - (FLTYPE)j))/((FLTYPE)k_num_steps);
                    double s = std::real(cpt.self_energy_precalc(omega, k, Ginv));
                    results << omega << "\t" << k[0] << "\t" << k[1] << "\t" << k[2] << "\t" << s << std::endl;
                }
            }
            else
            {
                //perform a k sweep at this omega
                if(cpt.H->graph.is_nambu)
                    b= cpt.nambu_G_k_sweep(omega+eps*img, cpt.H->graph.sym_points[i], cpt.H->graph.sym_points[(i+1)%num_sym], k_num_steps, self_energy_period);
                else
                    b= cpt.G_k_sweep(omega+eps*img, cpt.H->graph.sym_points[i], cpt.H->graph.sym_points[(i+1)%num_sym], k_num_steps, self_energy_period);
                int j = 0;
                //write results to file
                for(const std::complex<FLTYPE>& s: b)
                {
                    Vector3d k = (cpt.H->graph.sym_points[(i+1)%num_sym] * (FLTYPE)j + cpt.H->graph.sym_points[i]*((FLTYPE)k_num_steps - (FLTYPE)j))/((FLTYPE)k_num_steps);
                    results << omega << "\t" << k[0] << "\t" << k[1] << "\t" << k[2] << "\t" << -std::imag(s) << std::endl;
                    j++;
                }
            }
        }
    }
    results.close();
}
