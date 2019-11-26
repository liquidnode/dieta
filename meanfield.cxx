/*
This code is part of DIETA

Authored by Kirill Alpin
*/

#include "meanfield.h"

MeanField::MeanField(CPT* _cpt)
{
    cpt = _cpt;
}

/*
 * performs a mean field decoupling via iteration and current mean field parameters
 */
void MeanField::findFixpointWithCurrent()
{
    cpt->greens->set_eps(0.0);
    cpt->H->current_set = true;
    
    //set starting point with d variables
    cpt->H->densities = VectorXd::Zero(cpt->H->numsites);
    cpt->H->single_operator_densities = VectorXcd::Zero(cpt->H->numsites);
    VectorXcd single_last = VectorXcd::Zero(cpt->H->numsites);
    VectorXcd single_tmp = VectorXcd::Zero(cpt->H->numsites);
    
    
    cpt->H->h_densities = cpt->H->densities;
    cpt->H->h_current_densities = cpt->H->current_densities;
    
    std::cout << "Start densities: " << std::endl << cpt->H->densities << std::endl << std::endl;
    
    if(!cpt->H->graph.is_nambu)
    {
        std::stringstream ss;
        ss << "Non nambu current mean field not implemented.";
        throw std::invalid_argument(ss.str());
    }
    
    
    bool converged = false;
    VectorXd densitieslast = cpt->H->densities;
    VectorXd h_densitieslast = cpt->H->h_densities;
    MatrixXcd currentlast = cpt->H->current_densities;
    MatrixXcd h_currentlast = cpt->H->h_current_densities;
    MatrixXcd anomcurrentlast = cpt->H->anom_current_densities;
    bool first = true;
    
    int i = 0;
    while(!converged)
    {
        if(cpt->H->graph.cluster_dens)
        {
            std::cout << cpt->H->h_densities << std::endl;
            if(cpt->H->is_complex)
                ((Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>*)cpt->greens)->groundstate = Density::calcEverythingCluster<std::complex<FLTYPE>, VectorXcd, MatrixXcd>(cpt->H);
            else
                ((Greens<FLTYPE, FARRAY, MatrixXd>*)cpt->greens)->groundstate = Density::calcEverythingCluster<FLTYPE, FARRAY, MatrixXd>(cpt->H);
            
            cpt->H->h_current_densities = cpt->H->current_densities;
            cpt->H->h_densities = cpt->H->densities;
        }
        else if(cpt->H->graph.better_densities)
        {
            if(cpt->H->graph.max_particle != 1)
            {
                std::stringstream ss;
                ss << "Corrected densities are only possible when max_particle is equal to 1.";
                throw std::invalid_argument(ss.str());
            }
            cpt->recalculate_greens();
            MatrixXcd everything = Density::calcEverything(cpt, cpt->H->graph.is_nambu, true);
            
            //densities
            cpt->H->densities = everything.diagonal().real().head(cpt->H->numsites);
            VectorXd ph_densities = everything.diagonal().real().tail(cpt->H->numsites);
            cpt->H->densities = (cpt->H->densities + (VectorXd::Ones(cpt->H->numsites) - ph_densities)) / 2.0;
            cpt->H->h_densities = cpt->H->densities;
            
            //current
            cpt->H->current_densities = (everything.block(0,0,cpt->H->numsites,cpt->H->numsites) + everything.block(cpt->H->numsites,cpt->H->numsites,cpt->H->numsites,cpt->H->numsites)) / 2.0;
            
            std::complex<FLTYPE> mcurrent = 0.0;
            int count = 0;
            for(int n = 0; n < cpt->H->numsites; ++n)
                for(int m = 0; m < cpt->H->numsites; ++m)
                {
                    if(n != m)
                    {
                        mcurrent += cpt->H->current_densities(n,m);
                        count++;
                    }
                }
            mcurrent /= count;
            std::cout << "Mean current:\t" << mcurrent << std::endl;
            
            cpt->H->h_current_densities = cpt->H->current_densities;
            
            //anom current
            cpt->H->anom_current_densities = (everything.block(cpt->H->numsites,0,cpt->H->numsites,cpt->H->numsites) + everything.block(cpt->H->numsites,0,cpt->H->numsites,cpt->H->numsites).transpose()) / 2.0; //<aa> terms
            
            std::complex<FLTYPE> macurrent = 0.0;
            count = 0;
            for(int n = 0; n < cpt->H->numsites; ++n)
                for(int m = 0; m < cpt->H->numsites; ++m)
                {
                    if(n != m)
                    {
                        macurrent += cpt->H->anom_current_densities(n,m);
                        count++;
                    }
                }
            macurrent /= count;
            std::cout << "Mean anom current:\t" << macurrent << std::endl;
        }
        else
        {
            MatrixXcd everything = Density::calcEverything(cpt, cpt->H->graph.is_nambu, true);
            
            if(cpt->H->graph.max_particle == 1)
            {
                cpt->H->densities = everything.diagonal().real().head(cpt->H->numsites);
                cpt->H->h_densities = VectorXd::Ones(cpt->H->numsites) - everything.diagonal().real().tail(cpt->H->numsites);
            }
            else
            {
                cpt->H->densities = everything.diagonal().real().head(cpt->H->numsites);
                cpt->H->h_densities = everything.diagonal().real().tail(cpt->H->numsites) - VectorXd::Ones(cpt->H->numsites);
            }
            
            cpt->H->current_densities = everything.block(0,0,cpt->H->numsites,cpt->H->numsites);
            cpt->H->h_current_densities = everything.block(cpt->H->numsites,cpt->H->numsites,cpt->H->numsites,cpt->H->numsites);
            
            std::complex<FLTYPE> mcurrent = 0.0;
            int count = 0;
            for(int n = 0; n < cpt->H->numsites; ++n)
                for(int m = 0; m < cpt->H->numsites; ++m)
                {
                    if(n != m)
                    {
                        mcurrent += cpt->H->current_densities(n,m);
                        count++;
                    }
                }
            mcurrent /= count;
            std::cout << "Mean current:\t" << mcurrent << std::endl;
            
            //anom current
            cpt->H->anom_current_densities = everything.block(cpt->H->numsites,0,cpt->H->numsites,cpt->H->numsites); //<aa> terms
            
            std::complex<FLTYPE> macurrent = 0.0;
            count = 0;
            for(int n = 0; n < cpt->H->numsites; ++n)
                for(int m = 0; m < cpt->H->numsites; ++m)
                {
                    if(n != m)
                    {
                        macurrent += cpt->H->anom_current_densities(n,m);
                        count++;
                    }
                }
            macurrent /= count;
            std::cout << "Mean anom current:\t" << macurrent << std::endl;
        }
            
        if(cpt->H->graph.is_half_filling)
        {
            //force half filling
            double mean = cpt->H->densities.mean();
            std::cout << "Mean " << mean << std::endl;
            cpt->H->densities += (0.5 - mean) * VectorXd::Ones(cpt->H->densities.size());
        }
        
        if(cpt->H->graph.better_densities)
            std::cout << std::endl  << "Corrected particle densities: " << std::endl << cpt->H->densities << std::endl << std::endl;
        else
            std::cout << std::endl  << "Particle densities: " << std::endl << cpt->H->densities << std::endl << std::endl;
        
        if(!first && (densitieslast - cpt->H->densities).cwiseAbs().maxCoeff() < 1e-7 && (currentlast - cpt->H->current_densities).cwiseAbs().maxCoeff() < 1e-7 && (anomcurrentlast - cpt->H->anom_current_densities).cwiseAbs().maxCoeff() < 1e-7)
        {
            converged = true;
            break;
        }
        densitieslast = cpt->H->densities;
        h_densitieslast = cpt->H->h_densities;
        currentlast = cpt->H->current_densities;
        anomcurrentlast = cpt->H->anom_current_densities;
        first = false;
        
        i++;
        
    }
}

/*
 * performs a mean field decoupling via iteration using only particle densities
 */
void MeanField::findFixpoint()
{
    cpt->greens->set_eps(0.0);
    
    
    bool converged = false;
    VectorXd densitieslast = cpt->H->densities;
    cpt->H->h_densities = cpt->H->densities;
    bool first = true;
    while(!converged)
    {
        if(cpt->H->graph.cluster_dens)
        {
            if(cpt->H->is_complex)
                cpt->H->densities = Density::calcAllClusterDensity<std::complex<FLTYPE>, VectorXcd, MatrixXcd>(cpt->H, ((Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>*)cpt->greens)->groundstate);
            else
                cpt->H->densities = Density::calcAllClusterDensity<FLTYPE, FARRAY, MatrixXd>(cpt->H, ((Greens<FLTYPE, FARRAY, MatrixXd>*)cpt->greens)->groundstate);
        }
        else
        {
            cpt->recalculate_greens();
            if(cpt->H->graph.better_densities)
            {
                if(cpt->H->graph.max_particle != 1)
                {
                    std::stringstream ss;
                    ss << "Corrected densities are only possible when max_particle is equal to 1.";
                    throw std::invalid_argument(ss.str());
                }
                VectorXd alldensities = Density::calcAllDensity(cpt, cpt->H->graph.is_nambu, true);
                cpt->H->densities = alldensities.head(cpt->H->numsites);
                VectorXd ph_densities = alldensities.tail(cpt->H->numsites);
                std::cout << "Particle densities: " << std::endl << cpt->H->densities << std::endl << std::endl;
                std::cout << "Hole densities: " << std::endl<< ph_densities << std::endl;
                cpt->H->densities = (cpt->H->densities + (VectorXd::Ones(cpt->H->numsites) - ph_densities)) / 2.0;
            }
            else
                cpt->H->densities = Density::calcAllDensity(cpt, cpt->H->graph.is_nambu);
        }
            
        if(cpt->H->graph.is_half_filling)
        {
            //force half filling
            double mean = cpt->H->densities.mean();
            std::cout << "Mean " << mean << std::endl;
            cpt->H->densities += (0.5 - mean) * VectorXd::Ones(cpt->H->densities.size());
        }
        
        VectorXd target = cpt->H->densities;
        double mean = cpt->H->densities.mean();
        std::cout << "Mean " << mean << std::endl;
        cpt->H->densities = cpt->H->densities;
        cpt->H->h_densities = cpt->H->densities;
            
        if(cpt->H->graph.better_densities)
            std::cout << std::endl  << "Corrected particle densities: " << std::endl << cpt->H->densities << std::endl << std::endl;
        else
            std::cout << std::endl  << "Particle densities: " << std::endl << cpt->H->densities << std::endl << std::endl;
        
        if(!first && (densitieslast - target).cwiseAbs().maxCoeff() < 1e-7)
        {
            converged = true;
            break;
        }
        densitieslast = cpt->H->densities;
        first = false;
    }
}
