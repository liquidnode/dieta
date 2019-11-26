/*
This code is part of DIETA

Authored by Kirill Alpin
*/

#include "density.h"


double GkInt(double x, double y, void* data)
{
    InData* sdata = (InData*)data;
    CPT* cpt = sdata->cpt;
    std::complex<double> omega = ((InData*)data)->omega;
    MatrixXcd G = ((InData*)data)->G;
    
#ifndef PERIODIZATION
    if(sdata->nambu)
    {
        if(sdata->mean)
        {
            if(sdata->real)
                return std::real(cpt->cluster_nambu_G_precalc(omega, Vector3d(x,y,0), G).block(0,0,cpt->H->numsites,cpt->H->numsites).trace() / (double)(cpt->H->numsites));
            else
                return std::imag(cpt->cluster_nambu_G_precalc(omega, Vector3d(x,y,0), G).block(0,0,cpt->H->numsites,cpt->H->numsites).trace() / (double)(cpt->H->numsites));
        }
        else
        {
            if(sdata->real)
                return std::real(cpt->cluster_nambu_G_precalc(omega, Vector3d(x,y,0), G).block(0,0,cpt->H->numsites,cpt->H->numsites)(sdata->site,sdata->site));
            else
                return std::imag(cpt->cluster_nambu_G_precalc(omega, Vector3d(x,y,0), G).block(0,0,cpt->H->numsites,cpt->H->numsites)(sdata->site,sdata->site));
        }
    }
    else
    {
        if(sdata->mean)
        {
            if(sdata->real)
                return std::real(cpt->cluster_G_precalc(omega, Vector3d(x,y,0), G).trace() / (double)(cpt->H->numsites));
            else
                return std::imag(cpt->cluster_G_precalc(omega, Vector3d(x,y,0), G).trace() / (double)(cpt->H->numsites));
        }
        else
        {
            if(sdata->real)
                return std::real(cpt->cluster_G_precalc(omega, Vector3d(x,y,0), G)(sdata->site,sdata->site));
            else
                return std::imag(cpt->cluster_G_precalc(omega, Vector3d(x,y,0), G)(sdata->site,sdata->site));
        }
    }
#else
    if(!sdata->mean)
    {
        std::stringstream ss;
        ss << "Site density not supported with PERIODIZATION flag set.";
        throw std::invalid_argument(ss.str());
    }
    if(sdata->real)
        return std::real(cpt->G_with_self_energy_precalc(omega, Vector3d(x,y,0), G));
    else
        return std::imag(cpt->G_with_self_energy_precalc(omega, Vector3d(x,y,0), G));
#endif
}

double f(double x, void* data)
{
    const double pi = std::acos(-1);
    const std::complex<double> img(0, 1);
    InData* sdata = (InData*)data;
    std::complex<double> start = sdata->start;
    std::complex<double> end = sdata->end;
    std::complex<double> omega = start + (end - start) * x;
    CPT* cpt = sdata->cpt;
    sdata->omega = omega;
#ifdef PERIODIZATION
    if(sdata->nambu)
        sdata->G = cpt->greens->get_nambu_greens(omega).inverse();
    else
        sdata->G = cpt->greens->get_greens(omega).inverse();
#else
    if(sdata->nambu)
        sdata->G = cpt->greens->get_nambu_greens(omega);
    else
        sdata->G = cpt->greens->get_greens(omega);
#endif
    
    //integrate G over brillouin zone
    double init = 0.0;
    std::complex<double> val = 0.0;
    if(std::abs(std::imag(end - start)) > 1e-11)
    {
        sdata->real = true;
        val = BZIntegrator::bzIntegrate<double>(cpt, GkInt, sdata, init, false);
    }
    else
    {
        sdata->real = false;
        val = BZIntegrator::bzIntegrate<double>(cpt, GkInt, sdata, init, false);
    }
    val *= -img / (2.0*pi);
    val *= end - start; //jacobi
    if(std::real(val) != std::real(val))
    {
        std::cout << FRED("ALARM") << std::endl;
        return 0;
    }
    return std::real(val);
}


double senechal_integration_f(double omega, void* data)
{
    const double pi = std::acos(-1);
    const std::complex<double> img(0, 1);
	InData* sdata = new InData((InData*)data);
    CPT* cpt = sdata->cpt;
    if(sdata->inverse_omega)
        omega = 1.0 / omega;
    sdata->omega = img*omega;
    
    #ifdef PERIODIZATION
    if(sdata->nambu)
        sdata->G = cpt->greens->get_nambu_greens(sdata->omega).inverse();
    else
        sdata->G = cpt->greens->get_greens(sdata->omega).inverse();
    #else
    if(sdata->nambu)
        sdata->G = cpt->greens->get_nambu_greens(sdata->omega);
    else
        sdata->G = cpt->greens->get_greens(sdata->omega);
    #endif
    
    //integrate G over brillouin zone
    sdata->real = true;
    double init = 0.0;
    double val  = BZIntegrator::bzIntegrate<double>(cpt, GkInt, sdata, init, false);
    int b = cpt->H->numsites;
    if(sdata->mean)
        val -= std::real((1.0 / ((img*omega - 1.0)))*cpt->greens->omegaInfLimit().block(0,0,b,b).diagonal().mean());
    else
        val -= std::real((1.0 / ((img*omega - 1.0))))*std::real(cpt->greens->omegaInfLimit()(sdata->site,sdata->site));
    
    if(std::real(val) != std::real(val))
        return 0;
    if(sdata->inverse_omega)
        val *= omega * omega;
    
    delete sdata;
    return val;
}

/*
 * calculates the expectation value <n_(site)> or the mean particle density,
 * if (mean) is set to true
 */
double Density::calcDensity(CPT* cpt, int site, bool mean, bool nambu)
{
    const double pi = std::acos(-1);
    
    InData* data = new InData();
    data->cpt = cpt;
    data->mean = mean;
    data->site = site;
    data->onedim = cpt->H->graph.is_onedim;
    data->nambu = nambu;
    data->inverse_omega = false;
    data->is_quadri = cpt->greens->get_H()->graph.poly_bz.size() != 0;
    
#ifdef SENECHAL
    double freq = cpt->greens->get_high_energy()*4.0;
    double low_energy = gauss_legendre(128,senechal_integration_f,(void*)data,0.0,freq);
    data->inverse_omega = true;
    double hi_energy = gauss_legendre(128,senechal_integration_f,(void*)data,0.0,1.0/freq);
    delete data;
    
    if(cpt->greens->get_H()->graph.particleType == ParticleType::FERMION)
        return (low_energy + hi_energy) / pi;
    else
        return -(low_energy + hi_energy) / pi;
#else
    
    double width = 0.1;
    double eps_int = 1e-10;
    double freq = cpt->greens->get_high_energy()*1.5;
    
    data->start = std::complex<FLTYPE>(-eps_int,-width);
    data->end = std::complex<FLTYPE>(-eps_int,width);
    double A = gauss_legendre(128,f,(void*)data,0.0,1.0);
    
    data->start = std::complex<FLTYPE>(-eps_int,width);
    data->end = std::complex<FLTYPE>(-freq,width);
    double B = gauss_legendre(128,f,(void*)data,0.0,1.0);
    
    data->start = std::complex<FLTYPE>(-freq,width);
    data->end = std::complex<FLTYPE>(-freq,-width);
    double C = gauss_legendre(128,f,(void*)data,0.0,1.0);
    
    data->start = std::complex<FLTYPE>(-freq,-width);
    data->end = std::complex<FLTYPE>(-eps_int,-width);
    double D = gauss_legendre(128,f,(void*)data,0.0,1.0);
    
    return -(A + B + C + D);
#endif
}




//vector stuff







VectorXd GkIntV(double x, double y, void* data)
{
    InData* sdata = (InData*)data;
    CPT* cpt = sdata->cpt;
    std::complex<double> omega = ((InData*)data)->omega;
    MatrixXcd G = ((InData*)data)->G;
    
    int b = cpt->H->numsites;
    if(sdata->with_hole_dens)
        b += cpt->H->numsites;
    if(sdata->nambu || sdata->with_hole_dens)
    {
        if(sdata->real)
            return cpt->cluster_nambu_G_precalc(omega, Vector3d(x,y,0), G).block(0,0,b,b).diagonal().real();
        else
            return cpt->cluster_nambu_G_precalc(omega, Vector3d(x,y,0), G).block(0,0,b,b).diagonal().imag();
    }
    else
    {
        if(sdata->real)
            return cpt->cluster_G_precalc(omega, Vector3d(x,y,0), G).diagonal().real();
        else
            return cpt->cluster_G_precalc(omega, Vector3d(x,y,0), G).diagonal().imag();
    }
}

VectorXd senechal_integration_fV(double omega, void* data)
{
    const double pi = std::acos(-1);
    const std::complex<double> img(0, 1);
	InData* sdata = new InData((InData*)data);
    CPT* cpt = sdata->cpt;
    if(sdata->inverse_omega)
        omega = 1.0 / omega;
    sdata->omega = img*omega;
    
    if(sdata->nambu)
        sdata->G = cpt->greens->get_nambu_greens(sdata->omega);
    else
    {
        if(sdata->with_hole_dens)
        {
            sdata->G = MatrixXcd::Zero(cpt->H->numsites * 2, cpt->H->numsites * 2);
            sdata->G << cpt->greens->get_greens(sdata->omega), 
            MatrixXcd::Zero(cpt->H->numsites, cpt->H->numsites),
            MatrixXcd::Zero(cpt->H->numsites, cpt->H->numsites), 
            cpt->greens->bfh_sign(cpt->greens->get_greens(-sdata->omega));
        }
        else
            sdata->G = cpt->greens->get_greens(sdata->omega);
    }
    
    //integrate G over brillouin zone
    sdata->real = true;
    int b = cpt->H->numsites;
    if(sdata->with_hole_dens)
        b += cpt->H->numsites;
    VectorXd init = VectorXd::Zero(b);
    VectorXcd val = BZIntegrator::bzIntegrate<VectorXd>(cpt, GkIntV, sdata, init, false);
    val -= std::real(1.0 / ((sdata->omega - 1.0))) * cpt->greens->omegaInfLimit().block(0,0,b,b).diagonal();
    
    
    if(sdata->inverse_omega)
        val *= omega * omega;
    
    if(val.real() != val.real())
    {
        std::cout << FRED("ALARM") << std::endl;
        return VectorXd::Zero(b);
    }
    
    delete sdata;
    return val.real();
}

/*
 * calculates the expectation vector <n_i>, i.e. the particle density simultaneously for all sites
 */
VectorXd Density::calcAllDensity(CPT* cpt, bool nambu, bool with_hole_dens)
{
    const double pi = std::acos(-1);
    
    InData* data = new InData();
    data->cpt = cpt;
    data->onedim = cpt->H->graph.is_onedim;
    data->nambu = nambu;
    data->with_hole_dens = with_hole_dens;
    data->inverse_omega = false;
    data->is_quadri = cpt->H->graph.poly_bz.size() != 0;
    
    int b = cpt->H->numsites;
    if(data->with_hole_dens)
        b += cpt->H->numsites;
    
    VectorXd init = VectorXd::Zero(b);
    double freq = cpt->greens->get_high_energy()*4.0;
    VectorXd low_energy = gauss_legendreI<VectorXd>(256,senechal_integration_fV,(void*)data,0.0,freq, init, true);
    data->inverse_omega = true;
    VectorXd hi_energy = gauss_legendreI<VectorXd>(256,senechal_integration_fV,(void*)data,0.0,1.0/freq, init, true);
    delete data;
    
    if(cpt->greens->get_H()->graph.particleType == ParticleType::FERMION)
        return (low_energy + hi_energy) / pi;
    else
        return -(low_energy + hi_energy) / pi;
}





//matrix stuff









MatrixXcd GkIntM(double x, double y, void* data)
{
    InData* sdata = (InData*)data;
    CPT* cpt = sdata->cpt;
    std::complex<double> omega = ((InData*)data)->omega;
    MatrixXcd G = ((InData*)data)->G;
    
    int b = cpt->H->numsites;
    if(sdata->with_hole_dens)
        b += cpt->H->numsites;
    if(sdata->nambu || sdata->with_hole_dens)
        return cpt->cluster_nambu_G_precalc(omega, Vector3d(x,y,0), G).block(0,0,b,b);
    else
        return cpt->cluster_G_precalc(omega, Vector3d(x,y,0), G);
}

MatrixXcd senechal_integration_fM(double omega, void* data)
{
    const double pi = std::acos(-1);
    const std::complex<double> img(0, 1);
	InData* sdata = new InData((InData*)data);
    CPT* cpt = sdata->cpt;
    if(sdata->inverse_omega)
        omega = 1.0 / omega;
    sdata->omega = img*omega;
    
    if(sdata->nambu)
        sdata->G = cpt->greens->get_nambu_greens(sdata->omega);
    else
    {
        if(sdata->with_hole_dens)
        {
            sdata->G = MatrixXcd::Zero(cpt->H->numsites * 2, cpt->H->numsites * 2);
            sdata->G << cpt->greens->get_greens(sdata->omega), 
            MatrixXcd::Zero(cpt->H->numsites, cpt->H->numsites),
            MatrixXcd::Zero(cpt->H->numsites, cpt->H->numsites), 
            cpt->greens->bfh_sign(cpt->greens->get_greens(-sdata->omega));
        }
        else
            sdata->G = cpt->greens->get_greens(sdata->omega);
    }
    
    //integrate G over brillouin zone
    sdata->real = true;
    int b = cpt->H->numsites;
    if(sdata->with_hole_dens)
        b += cpt->H->numsites;
    
    MatrixXcd init = MatrixXcd::Zero(b, b);
    MatrixXcd val = BZIntegrator::bzIntegrate<MatrixXcd>(cpt, GkIntM, sdata, init, false);
    val -= std::real(1.0 / ((sdata->omega - 1.0))) * cpt->greens->omegaInfLimit().block(0,0,b,b);
    
    if(sdata->inverse_omega)
        val *= omega * omega;
    
    if(val.real() != val.real() || val.imag() != val.imag())
    {
        std::cout << FRED("ALARM") << std::endl;
        return MatrixXcd::Zero(b, b);
    }
    
    delete sdata;
    return val;
}

/*
 * calculates the expectation values <a^dag a>, with a being a vector a=(a_1,a_2,...)
 * if (nambu) is true, a is the vector a=(a_1,a_2,...,a^dag_1,a^dag_2,...)
 */
MatrixXcd Density::calcEverything(CPT* cpt, bool nambu, bool with_hole_dens)
{
    const double pi = std::acos(-1);
    
    InData* data = new InData();
    data->cpt = cpt;
    data->onedim = cpt->H->graph.is_onedim;
    data->nambu = nambu;
    data->with_hole_dens = with_hole_dens;
    data->is_quadri = cpt->H->graph.poly_bz.size() != 0;
    
    const std::complex<double> img(0, 1);
    
    int b = cpt->H->numsites;
    if(data->with_hole_dens)
        b += cpt->H->numsites;
    
    MatrixXcd init = MatrixXcd::Zero(b,b);
    
    double freq = cpt->greens->get_high_energy()*1.5;
    MatrixXcd low_energy = gauss_legendreI<MatrixXcd>(256,senechal_integration_fM,(void*)data,-freq,freq, init, true);
    data->inverse_omega = true;
    MatrixXcd hi_energy = gauss_legendreI<MatrixXcd>(128,senechal_integration_fM,(void*)data,0.0,1.0/freq, init, true);
    hi_energy += gauss_legendreI<MatrixXcd>(128,senechal_integration_fM,(void*)data,-1.0/freq,0.0, init, true);
    delete data;
    
    if(cpt->greens->get_H()->graph.particleType == ParticleType::FERMION)
        return (low_energy + hi_energy) / (2.0 * pi);
    else
        return -(low_energy + hi_energy) / (2.0 * pi);
}






//cluster stuff







/*
 * calculates the cluster particle density vector
 */
template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
VectorXd Density::calcAllClusterDensity(Hamiltonian* H, VECTYPE& groundstate)
{
    Lanczos<VALTYPE, VECTYPE, MATTYPE> lanczos(H);
    
    if(H->is_parity_conserving || H->is_particle_conserving)
    {
        std::stringstream ss;
        ss << "Cluster density for Hilbert space truncation not implemented.";
        throw std::invalid_argument(ss.str());
    }
    
    if(H->dim <= 256)
    {
        MATTYPE H_matrix = MATTYPE::Zero(H->dim, H->dim);
        for(int i = 0; i < H->dim; ++i)
        {
            VECTYPE state = VECTYPE::Zero(H->dim);
            state[i] = 1.0;
            H_matrix.col(i) = (*H) * state;
        }
        SelfAdjointEigenSolver<MATTYPE> es(H_matrix);
        MATTYPE evecs = es.eigenvectors();
        int sm = 0;
        FLTYPE genergy = std::real(es.eigenvalues()[sm]);
        for(int i = 0; i < es.eigenvalues().size(); ++i)
            if(std::real(es.eigenvalues()[i]) < genergy)
            {
                genergy = std::real(es.eigenvalues()[i]);
                sm = i;
            }
        groundstate = evecs.col(sm);
    }
    else
    {
        MATTYPE M = lanczos.get_tridiagonal(1234);
        FLTYPE genergy = 0.0;
        groundstate = lanczos.get_groundstate(1234, M, genergy);
    }
    
    VectorXd dens = VectorXd::Zero(H->numsites);
    for(int i = 0; i < H->numsites; ++i)
        dens[i] = std::real(groundstate.dot(H->nop_site(groundstate, i)) / std::sqrt(groundstate.dot(groundstate)));
    return dens;
}

template VectorXd Density::calcAllClusterDensity<FLTYPE, FARRAY, MatrixXd>(Hamiltonian* H, FARRAY& groundstate);
template VectorXd Density::calcAllClusterDensity<std::complex<FLTYPE>, VectorXcd, MatrixXcd>(Hamiltonian* H, VectorXcd& groundstate);


/*
 * calculates the expectation values <a^dag a>, with a being a vector a=(a_1,a_2,...) using the cluster groundstate
 */
template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
VECTYPE Density::calcEverythingCluster(Hamiltonian* H)
{
    Lanczos<VALTYPE, VECTYPE, MATTYPE> lanczos(H);
    
    if(H->is_parity_conserving || H->is_particle_conserving)
    {
        std::stringstream ss;
        ss << "Cluster density for Hilbert space truncation not implemented.";
        throw std::invalid_argument(ss.str());
    }
        
    
    VECTYPE groundstate;
    if(H->dim <= 256)
    {
        MATTYPE H_matrix = MATTYPE::Zero(H->dim, H->dim);
        for(int i = 0; i < H->dim; ++i)
        {
            VECTYPE state = VECTYPE::Zero(H->dim);
            state[i] = 1.0;
            H_matrix.col(i) = (*H) * state;
        }
        SelfAdjointEigenSolver<MATTYPE> es(H_matrix);
        MATTYPE evecs = es.eigenvectors();
        int sm = 0;
        FLTYPE genergy = std::real(es.eigenvalues()[sm]);
        for(int i = 0; i < es.eigenvalues().size(); ++i)
            if(std::real(es.eigenvalues()[i]) < genergy)
            {
                genergy = std::real(es.eigenvalues()[i]);
                sm = i;
            }
        groundstate = evecs.col(sm);
    }
    else
    {
        MATTYPE M = lanczos.get_tridiagonal(1234);
        FLTYPE genergy = 0.0;
        groundstate = lanczos.get_groundstate(1234, M, genergy);
    }
    
    //TODO noetig??
    groundstate /= std::sqrt(groundstate.dot(groundstate));
    
    for(int i = 0; i < H->numsites; ++i)
        H->densities[i] = std::real(groundstate.dot(H->nop_site(groundstate, i)));
    H->h_densities = H->densities;
    
    for(int n = 0; n < H->numsites; ++n)
        for(int m = 0; m < n; ++m)
        {
            VECTYPE state = VECTYPE::Zero(H->dim);
            for(unsigned int i = 0; i < state.size(); ++i)
            {
                unsigned int occup = H->index_to_occup(i);
                
                unsigned int siteA = m;
                unsigned int siteB = n;
                
                //forward
                unsigned int occup_work = occup;
                bool success = true;
                FLTYPE coeff = 1.0;
                
                success &= H->remove_particle(occup_work, siteA, coeff);
                success &= H->add_particle(occup_work, siteB, coeff);
                if(success)
                    state[H->occup_to_index(occup_work)] += coeff * groundstate[i];
            }
            H->current_densities(n,m) = groundstate.dot(state);
            H->current_densities(m,n) = std::conj(H->current_densities(n,m));
            H->h_current_densities(n,m) = H->current_densities(n,m);
            H->h_current_densities(m,n) = std::conj(H->current_densities(n,m));
            
            if(H->graph.is_nambu)
            {
                state = VECTYPE::Zero(H->dim);
                for(unsigned int i = 0; i < state.size(); ++i)
                {
                    unsigned int occup = H->index_to_occup(i);
                    
                    unsigned int siteA = m;
                    unsigned int siteB = n;
                    
                    //forward
                    unsigned int occup_work = occup;
                    bool success = true;
                    FLTYPE coeff = 1.0;
                    
                    success &= H->remove_particle(occup_work, siteA, coeff);
                    success &= H->remove_particle(occup_work, siteB, coeff);
                    if(success)
                        state[H->occup_to_index(occup_work)] += coeff * groundstate[i];
                }
                H->anom_current_densities(n,m) = groundstate.dot(state);
                H->anom_current_densities(m,n) = std::conj(H->anom_current_densities(n,m));
            }
        }
    return groundstate;
}

template FARRAY Density::calcEverythingCluster<FLTYPE, FARRAY, MatrixXd>(Hamiltonian* H);
template VectorXcd Density::calcEverythingCluster<std::complex<FLTYPE>, VectorXcd, MatrixXcd>(Hamiltonian* H);

/*
 * calculates the expectation values <a>, with a being a vector a=(a_1,a_2,...) using the cluster groundstate
 */
template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
void Density::calcClusterSingleOperatorDensity(CPT* cpt)
{
    cpt->H->single_operator_densities = VECTYPE::Zero(cpt->H->numsites);
    
    for(int j = 0; j < cpt->H->numsites; ++j)
    {
        VECTYPE state = VECTYPE::Zero(cpt->H->dim);
        for(unsigned int i = 0; i < state.size(); ++i)
        {
            unsigned int occup = cpt->H->index_to_occup(i);
            
            unsigned int occup_work = occup;
            bool success = true;
            FLTYPE coeff = 1.0;
            
            success &= cpt->H->remove_particle(occup_work, j, coeff);
            if(success)
            {
                state[cpt->H->occup_to_index(occup_work)] += coeff * ((Greens<VALTYPE, VECTYPE, MATTYPE>*)cpt->greens)->groundstate[i];
            }
        }
        cpt->H->single_operator_densities[j] = ((Greens<VALTYPE, VECTYPE, MATTYPE>*)cpt->greens)->groundstate.dot(state);
    }
    cpt->H->single_operator_densities_set = true;
}

template void Density::calcClusterSingleOperatorDensity<FLTYPE, FARRAY, MatrixXd>(CPT* cpt);
template void Density::calcClusterSingleOperatorDensity<std::complex<FLTYPE>, VectorXcd, MatrixXcd>(CPT* cpt);
