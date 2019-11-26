/*
This code is part of DIETA

Authored by Kirill Alpin
*/

#include "potential.h"
#include <unsupported/Eigen/MatrixFunctions>


bool Potential::cmplx_eigenval = false;
double Potential::max_imag_eigenval = 0.0;

//for SENECHAL_ integration
double grandPotentialInt(double x, double y, void* data)
{
    InDataPotential* sdata = (InDataPotential*)data;
    CPT* cpt = sdata->cpt;
    std::complex<double> omega = ((InDataPotential*)data)->omega;
    MatrixXcd G = ((InDataPotential*)data)->G;
    
    Vector3d k = Vector3d(x, y, 0);
    
    if(sdata->nambu)
    {
        //1-VG
        //return std::log(std::abs((MatrixXcd::Identity(cpt->H->numsites*2,cpt->H->numsites*2) - cpt->vmatrix->get_nambu_V(k) * G).determinant()))/2.0;
        
        G = cpt->cluster_nambu_G_precalc(omega, k, G);
    }
    else
    {
        //1-VG
        //return std::log(std::abs((MatrixXcd::Identity(cpt->H->numsites,cpt->H->numsites) - cpt->vmatrix->getV(k) * G).determinant()));
        
        G = cpt->cluster_G_precalc(omega, k, G);
    }
    
    return std::log(std::abs(sdata->detG / G.determinant())) / (sdata->nambu ? 2.0 : 1.0);
}

//for SENECHAL_ integration
double senechal_integration_grandPotential(double omega, void* data)
{
    const double pi = std::acos(-1);
    const std::complex<double> img(0, 1);
	InDataPotential* sdata = new InDataPotential((InDataPotential*)data);
    CPT* cpt = sdata->cpt;
    if(sdata->inverse_omega)
        omega = 1.0 / omega;
    sdata->omega = img*omega;
    
    if(sdata->nambu)
        sdata->G = cpt->greens->get_nambu_greens(sdata->omega);
    else
        sdata->G = cpt->greens->get_greens(sdata->omega);
    sdata->detG = sdata->G.determinant();
    
    //integrate G over brillouin zone
    double init = 0.0;
    double val = BZIntegrator::bzIntegrate<double>(cpt, grandPotentialInt, (void*)sdata, init, false);
    
    if(sdata->inverse_omega)
        val *= omega * omega;
    
    if(std::real(val) != std::real(val))
        return 0;
    delete sdata;
    
    return val;
}





//function used in contour integration
double grandPotentialInt_contour(double x, double y, void* data)
{
    InDataPotential* sdata = new InDataPotential((InDataPotential*)data);
    CPT* cpt = sdata->cpt;
    std::complex<double> omega = sdata->omega;
    MatrixXcd G = sdata->G;
    
    Vector3d k = Vector3d(x, y, 0);
    
    std::complex<double> val = 0.0;
    
#ifdef PERIODIZED
    //experimental code!!
    
    MatrixXcd Ginv = sdata->Ginv;
    //Full periodization
    /*
    std::complex<double> cGperiod = 0.0;
    for(int m = 0; m < (sdata->nambu ? cpt->H->numsites * 2 : cpt->H->numsites); ++m)
        for(int n = 0; n < (sdata->nambu ? cpt->H->numsites * 2 : cpt->H->numsites); ++n)
        {
            cGperiod += G(m,n) * std::exp(-std::complex<FLTYPE>(0,1)*k.dot(cpt->H->graph.site_pos[m%cpt->H->numsites].pos-cpt->H->graph.site_pos[n%cpt->H->numsites].pos));
        }
    cGperiod /= (sdata->nambu ? cpt->H->numsites * 2 : cpt->H->numsites);
    val -= std::log(cpt->G_with_self_energy_precalc(omega, k, Ginv)/cGperiod);
    */
    
    //partial periodization
    /*MatrixXcd cGperiod = MatrixXcd::Zero((sdata->nambu ? 2 : 1),(sdata->nambu ? 2 : 1));
    for(int m = 0; m < (sdata->nambu ? cpt->H->numsites * 2 : cpt->H->numsites); ++m)
        for(int n = 0; n < (sdata->nambu ? cpt->H->numsites * 2 : cpt->H->numsites); ++n)
        {
            cGperiod(m/cpt->H->numsites,n/cpt->H->numsites) += G(m,n) * std::exp(-std::complex<FLTYPE>(0,1)*k.dot(cpt->H->graph.site_pos[m%cpt->H->numsites].pos-cpt->H->graph.site_pos[n%cpt->H->numsites].pos));
        }
    cGperiod /= cpt->H->numsites;
    val -= std::log(cpt->partial_ft_G_with_self_energy_precalc(omega, k, Ginv).determinant()/sdata->detG);*/
    
    //self energy periodization of cluster greens
    /*
    MatrixXcd G0prime = cpt->greens->get_noninteracting_inverse_greens(omega, cpt->H->graph.is_nambu, false);
    std::complex<FLTYPE> lret(0,0);
    for(int m = 0; m < (cpt->H->graph.is_nambu ? cpt->H->numsites * 2 : cpt->H->numsites); ++m)
        for(int n = 0; n < (cpt->H->graph.is_nambu ? cpt->H->numsites * 2 : cpt->H->numsites); ++n)
        {
            lret += G0prime(m,n) * std::exp(-std::complex<FLTYPE>(0,1)*k.dot(cpt->H->graph.site_pos[m%cpt->H->numsites].pos-cpt->H->graph.site_pos[n%cpt->H->numsites].pos));
        }
    lret /= (cpt->H->graph.is_nambu ? cpt->H->numsites * 2 : cpt->H->numsites);
    std::complex<double> cGperiod = 1.0/((1.0/lret) - cpt->self_energy_precalc(omega, k, Ginv));
    val -= std::log(cpt->G_with_self_energy_precalc(omega, k, Ginv)/cGperiod);*/
    
    //partial self energy periodization of cluster greens
    MatrixXcd G0prime = cpt->greens->get_noninteracting_inverse_greens(omega, cpt->H->graph.is_nambu, false);
    MatrixXcd cGperiod = MatrixXcd::Zero((sdata->nambu ? 2 : 1),(sdata->nambu ? 2 : 1));
    for(int m = 0; m < (sdata->nambu ? cpt->H->numsites * 2 : cpt->H->numsites); ++m)
        for(int n = 0; n < (sdata->nambu ? cpt->H->numsites * 2 : cpt->H->numsites); ++n)
        {
            cGperiod(m/cpt->H->numsites,n/cpt->H->numsites) += G0prime(m,n) * std::exp(-std::complex<FLTYPE>(0,1)*k.dot(cpt->H->graph.site_pos[m%cpt->H->numsites].pos-cpt->H->graph.site_pos[n%cpt->H->numsites].pos));
        }
    cGperiod /= cpt->H->numsites;
    cGperiod = (cGperiod.inverse() - cpt->partial_ft_self_energy_precalc(omega, k, Ginv)).inverse();
    val -= std::log(cpt->partial_ft_G_with_self_energy_precalc(omega, k, Ginv).determinant()/cGperiod.determinant());
#else
    if(sdata->nambu)
        G = cpt->cluster_nambu_G_precalc(omega, k, G);
    else
        G = cpt->cluster_G_precalc(omega, k, G);
    
    val -= std::log(G.determinant()/sdata->detG);
#endif
    
    delete sdata;
    
    if(sdata->real)
        return std::real(val);
    else
        return std::imag(val);
}

//function used in contour integration
double contour_integration_grandPotential(double x, void* data)
{
    const double pi = std::acos(-1);
    const std::complex<double> img(0, 1);
	InDataPotential* sdata = new InDataPotential((InDataPotential*)data);
    std::complex<double> start = sdata->start;
    std::complex<double> end = sdata->end;
    std::complex<double> omega = start + (end - start) * x;
    CPT* cpt = sdata->cpt;
    sdata->omega = omega;
    
    if(sdata->nambu)
        sdata->G = cpt->greens->get_nambu_greens(sdata->omega);
    else
        sdata->G = cpt->greens->get_greens(sdata->omega);
    
    sdata->detG = sdata->G.determinant();
    
#ifdef PERIODIZED
    sdata->Ginv = sdata->G.inverse();
#endif
    
    std::complex<double> val = 0.0;
    if(std::abs(std::imag(end - start)) > 1e-11)
        sdata->real = true;
    else
        sdata->real = false;
    
    double init = 0.0;
    val = BZIntegrator::bzIntegrate<double>(cpt, grandPotentialInt_contour, (void*)sdata, init, false);
    
    if(!sdata->real)
        val *= img;
    
        
    val *= -img / (2.0*pi);
    val *= end - start; //jacobi
    if(sdata->nambu)
        val /= 2.0;
    if(std::real(val) != std::real(val))
    {
        std::cout << FRED("ALARM") << std::endl;
        return 0;
    }
    
    delete sdata;
    
    return std::real(val);
}


//function used in exact integration
double grandPotentialInt_ex(double x, double y, void* data)
{
    InDataPotential* sdata = new InDataPotential((InDataPotential*)data);
    CPT* cpt = sdata->cpt;
    
    Vector3d k = Vector3d(x, y, 0);
    
    double ret = 0.0;
    
    MatrixXcd V;
    if(sdata->nambu)
        V = sdata->cpt->vmatrix->get_nambu_V(k);
    else
        V = sdata->cpt->vmatrix->getV(k);
    MatrixXcd Q = sdata->cpt->greens->get_Q();
    MatrixXd Om = sdata->cpt->greens->get_Omega();
    MatrixXcd L = sdata->cpt->greens->get_Lambda() + Om * Q.adjoint() * V * Q;
    
    ComplexEigenSolver<MatrixXcd> es(L);
    VectorXcd eigs = es.eigenvalues();
    for(int i = 0; i < eigs.size(); ++i)
    {
        if(std::abs(std::imag(eigs[i])) > 1e-8)
        {
            //warning: not thread safe
            sdata->cmplx_eigenval = true;
            sdata->max_imag_eigenval = std::max(sdata->max_imag_eigenval, std::abs(std::imag(eigs[i])));
        }
        if(std::real(eigs[i]) < 0.0)
            ret += std::real(eigs[i]);
    }
        
    for(auto d : sdata->neg_pols)
        ret -= d;
    
    delete sdata;
    return ret;
}

/*
 * calculates the Grand Potential Omega of a system
 */
double Potential::calcGrandPotential(CPT* cpt)
{
    cmplx_eigenval = false;
    bool nambu = cpt->H->graph.is_nambu;
    //check if zero frequency mode is needed
    bool zero_freq = false;
    MatrixXcd mu_prime;
    if(nambu)
        mu_prime = MatrixXcd::Zero(cpt->H->numsites * 2, cpt->H->numsites * 2);
    else
        mu_prime = MatrixXcd::Zero(cpt->H->numsites, cpt->H->numsites);
    //get mu_prime
    for(auto const& edge: cpt->H->graph.inter_edges)
    {
        Vector3d translation = cpt->H->graph.translations[edge.trans_index].pos;
        if(translation == Vector3d(0,0,0) && std::abs(cpt->H->graph.params[edge.ind]) > 1e-10)
        {
            zero_freq = true;
            //integration doesn't give 0
            switch(edge.type)
            {
                case 'A': //hopping term c^dag*c+c*c^dag
                    {
                        mu_prime(edge.b,edge.a) += cpt->H->graph.params[edge.ind];
                        if(edge.a != edge.b || translation.cwiseAbs().maxCoeff() > 1e-10)
                            mu_prime(edge.a,edge.b) += std::conj(cpt->H->graph.params[edge.ind]);
                        if(nambu)
                        {
                            mu_prime(edge.b+cpt->H->numsites,edge.a+cpt->H->numsites) += cpt->H->graph.params[edge.ind];
                            if(edge.a != edge.b || translation.cwiseAbs().maxCoeff() > 1e-10)
                                mu_prime(edge.a+cpt->H->numsites,edge.b+cpt->H->numsites) += std::conj(cpt->H->graph.params[edge.ind]);
                        }
                    }
                    break;
                case 'B': //hopping term c*c+c^dag*c^dag
                    {
                        if(nambu)
                        {
                            mu_prime(edge.b,edge.a+cpt->H->numsites) += cpt->H->graph.params[edge.ind];
                            if(edge.a != edge.b || translation.cwiseAbs().maxCoeff() > 1e-10)
                                mu_prime(edge.a,edge.b+cpt->H->numsites) += std::conj(cpt->H->graph.params[edge.ind]);
                            mu_prime(edge.b+cpt->H->numsites,edge.a) += cpt->H->graph.params[edge.ind];
                            if(edge.a != edge.b || translation.cwiseAbs().maxCoeff() > 1e-10)
                                mu_prime(edge.a+cpt->H->numsites,edge.b) += std::conj(cpt->H->graph.params[edge.ind]);
                        }
                    }
                    break;
                case 'F':
                    {
                        mu_prime(edge.a,edge.a) += cpt->H->graph.params[edge.ind];
                        if(nambu)
                            mu_prime(edge.a+cpt->H->numsites,edge.a+cpt->H->numsites) += cpt->H->graph.params[edge.ind];
                    }
                    break;
                case 'G':
                    {
                        mu_prime(edge.a,edge.a) += cpt->H->graph.params[edge.ind];
                    }
                    break;
                case 'H':
                    {
                        if(nambu)
                            mu_prime(edge.a+cpt->H->numsites,edge.a+cpt->H->numsites) += cpt->H->graph.params[edge.ind];
                    }
                    break;
                default:
                    std::stringstream ss;
                    ss << "Edge type " << edge.type << " not implemented.";
                    throw std::invalid_argument(ss.str());
                    break;
            }
        }
    }
    
    //calculate greens function
    cpt->greens->set_eps(0.0);
    cpt->recalculate_greens();
    
    //set integration mode
    if(cpt->H->graph.intMode == GPIntegrationMode::AUTOMATIC)
    {
        if(cpt->greens->is_full_diag())
        {
            if(cpt->H->dim <= 64)
                cpt->H->graph.intMode = GPIntegrationMode::EXACT;
            else
                cpt->H->graph.intMode = GPIntegrationMode::CONTOUR;
        }
        else
        {
            cpt->H->graph.intMode = GPIntegrationMode::CONTOUR;
        }
    }
    
    
    double clusterGrandPotential = cpt->greens->get_groundstate_energy(); //-mu*n ist schon enthalten, da Hamiltonian mu*adag*a terme enthaelt
    
    //TrG' of bath sites
    for(int i = 0; i < cpt->H->graph.num_bath;++i)
        if(cpt->H->bath_mu[i] < 0.0)
            clusterGrandPotential += cpt->H->bath_mu[i];
        
    //add normal ordering term if Nambu formalism is used
    if(nambu)
    {
#ifdef PERIODIZED
        MatrixXcd Ginv = cpt->greens->get_nambu_greens(1.0e8);
        Ginv = Ginv.inverse();
        clusterGrandPotential += 0.5*std::real(mu_prime.block(cpt->H->numsites, cpt->H->numsites, cpt->H->numsites, cpt->H->numsites).diagonal().mean()*
            1.0e8*(-cpt->partial_ft_G_with_self_energy_precalc(1.0e8, Vector3d(0,0,0), Ginv)(1,1)));
        std::cout << 1.0e8*(-cpt->partial_ft_G_with_self_energy_precalc(1.0e8, Vector3d(0,0,0), Ginv)(1,1)) << std::endl;
#else
        clusterGrandPotential += 0.5*std::real((mu_prime.block(cpt->H->numsites, cpt->H->numsites, cpt->H->numsites, cpt->H->numsites)*cpt->greens->omegaInfLimit().block(cpt->H->numsites, cpt->H->numsites, cpt->H->numsites, cpt->H->numsites)).trace());
        std::cout << "Normal ordering term:\t"<< 0.5*std::real((mu_prime.block(cpt->H->numsites, cpt->H->numsites, cpt->H->numsites, cpt->H->numsites)*cpt->greens->omegaInfLimit().block(cpt->H->numsites, cpt->H->numsites, cpt->H->numsites, cpt->H->numsites)).trace()) << std::endl;
#endif
    }
    

    if(cpt->H->graph.particleType != ParticleType::FERMION)
    {
        //add "Bosonic Self-energy Functional Theory" boson condensation terms
        if(!cpt->H->graph.F_field_zero_in_continuum_limit)
        {
            if(cpt->H->graph.real_F_field_is_mean_field)
            {
                std::cout << "Real F field is mean field." << std::endl;
                //phi=phi_mean
                cpt->H->setFields();
                VectorXcd F_field = cpt->H->F_field;
                VectorXcd real_F_field = cpt->H->meanF_field;
                
                if(cpt->H->is_complex)
                    Density::calcClusterSingleOperatorDensity<std::complex<FLTYPE>, VectorXcd, MatrixXcd>(cpt);
                else
                    Density::calcClusterSingleOperatorDensity<FLTYPE, FARRAY, MatrixXd>(cpt);
                cpt->H->single_operator_densities_set = false;
                
                MatrixXcd G_0prime_inv = cpt->greens->get_noninteracting_inverse_greens(0, cpt->greens->get_H()->graph.is_nambu, false);
                MatrixXcd G_0_inv = cpt->cluster_noninteracting_inverse_G(0, Vector3d(0,0,0));
                if(cpt->greens->get_H()->graph.is_nambu)
                {
                    VectorXcd tmp = cpt->H->single_operator_densities;
                    cpt->H->single_operator_densities = VectorXcd(cpt->H->numsites * 2);
                    cpt->H->single_operator_densities << tmp, tmp.conjugate();
                    tmp = F_field;
                    F_field = VectorXcd(cpt->H->numsites * 2);
                    F_field << tmp, tmp.conjugate();
                    tmp = real_F_field;
                    real_F_field = VectorXcd(cpt->H->numsites * 2);
                    real_F_field << tmp, tmp.conjugate();
                }
                VectorXcd FminusSigma = G_0prime_inv * cpt->H->single_operator_densities;
                VectorXcd Sigma = F_field - FminusSigma;
                if(std::abs(G_0_inv.determinant()) > 1e-10 && std::abs(G_0prime_inv.determinant()) > 1e-10)
                {
                    clusterGrandPotential += (nambu ? 0.5 : 1.0) * std::real(((real_F_field - Sigma).adjoint() * G_0_inv.inverse() * (real_F_field - Sigma) - FminusSigma.adjoint() * G_0prime_inv.inverse() * FminusSigma)(0,0));
                }
                else
                {
                    if((real_F_field - Sigma).cwiseAbs().maxCoeff() > 1e-9 || FminusSigma.cwiseAbs().maxCoeff() > 1e-9)
                    {
                        std::stringstream ss;
                        ss << "Gapless phase detected. Please break the symmetry!";
                        throw std::invalid_argument(ss.str());
                    }
                }
            }
            else
            {
                //phi=phi'
                if(cpt->H->is_complex)
                    Density::calcClusterSingleOperatorDensity<std::complex<FLTYPE>, VectorXcd, MatrixXcd>(cpt);
                else
                    Density::calcClusterSingleOperatorDensity<FLTYPE, FARRAY, MatrixXd>(cpt);
                cpt->H->single_operator_densities_set = false;
                
                MatrixXcd G_0prime_inv = cpt->greens->get_noninteracting_inverse_greens(0, cpt->greens->get_H()->graph.is_nambu, false);
                MatrixXcd G_0_inv = cpt->cluster_noninteracting_inverse_G(0, Vector3d(0,0,0));
                if(cpt->greens->get_H()->graph.is_nambu)
                {
                    VectorXcd tmp = cpt->H->single_operator_densities;
                    cpt->H->single_operator_densities = VectorXcd(cpt->H->numsites * 2);
                    cpt->H->single_operator_densities << tmp, tmp.conjugate();
                }
                VectorXcd FminusSigma = G_0prime_inv * cpt->H->single_operator_densities;
                if(std::abs(G_0_inv.determinant()) > 1e-10 && std::abs(G_0prime_inv.determinant()) > 1e-10)
                {
                    clusterGrandPotential += (nambu ? 0.5 : 1.0) * std::real((FminusSigma.adjoint() * (G_0_inv.inverse() - G_0prime_inv.inverse()) * FminusSigma)(0,0));
                }
                else
                {
                    if(FminusSigma.cwiseAbs().maxCoeff() > 1e-9)
                    {
                        std::stringstream ss;
                        ss << "Gapless phase detected. Please break the symmetry!";
                        throw std::invalid_argument(ss.str());
                    }
                }
            }
        }
        else
        {
            std::cout << "Real F field is zero." << std::endl;
            //phi=0
            cpt->H->setFields();
            VectorXcd F_field = cpt->H->F_field;
            if(cpt->H->is_complex)
                Density::calcClusterSingleOperatorDensity<std::complex<FLTYPE>, VectorXcd, MatrixXcd>(cpt);
            else
                Density::calcClusterSingleOperatorDensity<FLTYPE, FARRAY, MatrixXd>(cpt);
            cpt->H->single_operator_densities_set = false;
            
            MatrixXcd G_0prime_inv = cpt->greens->get_noninteracting_inverse_greens(0, cpt->greens->get_H()->graph.is_nambu);
            MatrixXcd G_0_inv = cpt->cluster_noninteracting_inverse_G(0, Vector3d(0,0,0));
            if(cpt->greens->get_H()->graph.is_nambu)
            {
                VectorXcd tmp = cpt->H->single_operator_densities;
                cpt->H->single_operator_densities = VectorXcd(cpt->H->numsites * 2);
                cpt->H->single_operator_densities << tmp, tmp.conjugate();
                tmp = F_field;
                F_field = VectorXcd(cpt->H->numsites * 2);
                F_field << tmp, tmp.conjugate();
            }
            VectorXcd FminusSigma = G_0prime_inv * cpt->H->single_operator_densities;
            VectorXcd Sigma = F_field - FminusSigma;
            if(std::abs(G_0_inv.determinant()) > 1e-13 && std::abs(G_0prime_inv.determinant()) > 1e-13)
            {
                clusterGrandPotential += (nambu ? 0.5 : 1.0) * std::real((Sigma.adjoint() * G_0_inv.inverse() * Sigma - FminusSigma.adjoint() * G_0prime_inv.inverse() * FminusSigma)(0,0));
            }
            else
            {
                if(Sigma.cwiseAbs().maxCoeff() > 1e-9 || FminusSigma.cwiseAbs().maxCoeff() > 1e-9)
                {
                    std::stringstream ss;
                    ss << "Gapless phase detected. Please break the symmetry!";
                    throw std::invalid_argument(ss.str());
                }
            }
        }
    }
    
    const double pi = std::acos(-1);
    InDataPotential* data = new InDataPotential();
    data->cpt = cpt;
    data->inverse_omega = false;
    data->onedim = cpt->H->graph.is_onedim;
    data->nambu = nambu;
    data->is_quadri = (cpt->H->graph.poly_bz.size() != 0);
    
    if(cpt->H->graph.intMode == GPIntegrationMode::EXACT)
    {
#ifdef PERIODIZED
        std::stringstream ss;
        ss << "Exact method not possible with self-energy periodization.";
        throw std::invalid_argument(ss.str());
#endif
        if(cpt->H->graph.num_bath > 0)
        {
            std::stringstream ss;
            ss << "Exact method not possible with bath sites.";
            throw std::invalid_argument(ss.str());
        }
        
        //trln(-G')
        double sum_neg_pols = 0.0;
        std::vector<double> neg_pols = cpt->greens->get_negative_poles();
        for(auto d : neg_pols)
        {
            if(d > -100.0)
                sum_neg_pols += d;
            else
                data->neg_pols.push_back(d);
        }
            
        //trln(-G)
        double init = 0.0;
        double sum_neg_pols_cpt = BZIntegrator::bzIntegrate<double>(cpt, grandPotentialInt_ex, data, init, true);
        
        cmplx_eigenval = data->cmplx_eigenval;
        max_imag_eigenval = data->max_imag_eigenval;
        
        if(cpt->H->graph.particleType == ParticleType::FERMION)
            return clusterGrandPotential + (sum_neg_pols_cpt - sum_neg_pols) / (nambu ? 2.0 : 1.0);
        else
            return clusterGrandPotential + (sum_neg_pols - sum_neg_pols_cpt) / (nambu ? 2.0 : 1.0);
    }
    
    
    if(cpt->H->graph.intMode == GPIntegrationMode::CONTOUR)
    {
#ifdef CONTOUR_
        double zero_freq_contrib = 0.0;
        if(zero_freq)
        {
            //G(w->inf)_a,b=???/w
            MatrixXcd greensLimit = cpt->greens->omegaInfLimit();
            zero_freq_contrib = std::real((mu_prime * greensLimit).trace());
            
            std::cout << "Zero energy contribution:\t" << zero_freq_contrib << std::endl;
        }
        data->zero_freq_contrib = zero_freq_contrib;
            
        double width = cpt->greens->get_high_energy();
        double eps_int = 1e-10;
        double freq = cpt->greens->get_high_energy()*4.0;
        
        data->start = std::complex<FLTYPE>(-eps_int,-width);
        data->end = std::complex<FLTYPE>(-eps_int,width);
        double A = gauss_legendre(256,contour_integration_grandPotential,(void*)data,0.0,1.0, true);
        
        data->start = std::complex<FLTYPE>(-eps_int,width);
        data->end = std::complex<FLTYPE>(-freq,width);
        double B = gauss_legendre(512,contour_integration_grandPotential,(void*)data,0.0,1.0, true);
        
        data->start = std::complex<FLTYPE>(-freq,width);
        data->end = std::complex<FLTYPE>(-freq,-width);
        double C = gauss_legendre(512,contour_integration_grandPotential,(void*)data,0.0,1.0, true);
        
        data->start = std::complex<FLTYPE>(-freq,-width);
        data->end = std::complex<FLTYPE>(-eps_int,-width);
        double D = gauss_legendre(256,contour_integration_grandPotential,(void*)data,0.0,1.0, true);
    
        if(cpt->H->graph.particleType == ParticleType::FERMION)
            return clusterGrandPotential - A - B - C - D;
        else
            return clusterGrandPotential + A + B + C + D;
#endif
    
    
#ifdef SENECHAL_
        //does not work properly. reason unknown
        
        double freq = cpt->greens->get_high_energy()*4.0;
        double low_energy = gauss_legendre(256,senechal_integration_grandPotential,(void*)data,0.0,freq, true);
        //low_energy += gauss_legendre(128,senechal_integration_grandPotential,(void*)data,-freq,0.0);
        data->inverse_omega = true;
        double hi_energy = gauss_legendre(256,senechal_integration_grandPotential,(void*)data,0.0,1.0/freq, true);
        //hi_energy += gauss_legendre(128,senechal_integration_grandPotential,(void*)data,-1.0/freq,0.0);
        delete data;
        
        double zero_freq_contrib = 0.0;
        if(zero_freq)
        {
            //G(w->inf)_a,b=???/w
            MatrixXcd greensLimit = cpt->greens->omegaInfLimit();
            zero_freq_contrib -= std::real((mu_prime * greensLimit).trace());
            if(nambu)
                zero_freq_contrib /= 2.0;
            
            std::cout << "Zero energy contribution:\t" << zero_freq_contrib << std::endl;
        }
        
        return clusterGrandPotential - ((low_energy + hi_energy) / (pi)) + zero_freq_contrib;
#endif
    }
    
    std::stringstream ss;
    ss << "This should not happen. Check intMode.";
    throw std::invalid_argument(ss.str());
    return 42;
}



//cluster stuff




template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
double Potential::calcClusterGrandPotential(Hamiltonian H)
{
    Lanczos<VALTYPE, VECTYPE, MATTYPE> lanczos(&H);

    if(H.is_parity_conserving || H.is_particle_conserving)
    {
        std::stringstream ss;
        ss << "Cluster grand potential for Hilbert space truncation not implemented.";
        throw std::invalid_argument(ss.str());
    }
    
    if(H.dim <= 256)
    {
        MATTYPE H_matrix = MATTYPE::Zero(H.dim, H.dim);
        for(int i = 0; i < H.dim; ++i)
        {
            VECTYPE state = VECTYPE::Zero(H.dim);
            state[i] = 1.0;
            H_matrix.col(i) = H * state;
        }
        SelfAdjointEigenSolver<MATTYPE> es(H_matrix);
        return std::real(es.eigenvalues()[0]);
    }
    else
    {
        MATTYPE M = lanczos.get_tridiagonal(1234);
        FLTYPE genergy = 0.0;
        lanczos.get_groundstate(1234, M, genergy);
        return std::real(genergy);
    }
}

template double Potential::calcClusterGrandPotential<FLTYPE, FARRAY, MatrixXd>(Hamiltonian H);
template double Potential::calcClusterGrandPotential<std::complex<FLTYPE>, VectorXcd, MatrixXcd>(Hamiltonian H);

template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
double Potential::calcClusterFreeEnergy(Hamiltonian H)
{
    Lanczos<VALTYPE, VECTYPE, MATTYPE> lanczos(&H);

    if(H.is_parity_conserving || H.is_particle_conserving)
    {
        std::stringstream ss;
        ss << "Cluster free energy for Hilbert space truncation not implemented.";
        throw std::invalid_argument(ss.str());
    }
    
    if(H.dim <= 256)
    {
        MATTYPE H_matrix = MATTYPE::Zero(H.dim, H.dim);
        for(int i = 0; i < H.dim; ++i)
        {
            VECTYPE state = VECTYPE::Zero(H.dim);
            state[i] = 1.0;
            H_matrix.col(i) = H * state;
        }
        return -std::log(std::abs((-100.0*H_matrix).exp().trace()))/100.0;
    }
    else
    {
        std::stringstream ss;
        ss << "Cluster free energy not implemented for big hamiltonian.";
        throw std::invalid_argument(ss.str());
        return 0;
    }
}

template double Potential::calcClusterFreeEnergy<FLTYPE, FARRAY, MatrixXd>(Hamiltonian H);
template double Potential::calcClusterFreeEnergy<std::complex<FLTYPE>, VectorXcd, MatrixXcd>(Hamiltonian H);



