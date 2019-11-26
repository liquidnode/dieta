/*
This code is part of DIETA

Authored by Kirill Alpin
*/

#include "cpt.h"

/*
 * initialize the CPT main class
 */
CPT::CPT(std::string graph_filename, std::string params_filename, std::string lanczos_params_filename, bool load, bool precalc)
{
    H = new Hamiltonian(graph_filename, params_filename);
    if(!load)
    {
        if(H->is_complex)
            greens = new Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>(H, precalc);
        else
            greens = new Greens<FLTYPE, FARRAY, MatrixXd>(H, precalc);
        if(lanczos_params_filename != "")
            greens->save_params(lanczos_params_filename);
    }
    else
    {
        if(H->is_complex)
            greens = new Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>(H, lanczos_params_filename);
        else
            greens = new Greens<FLTYPE, FARRAY, MatrixXd>(H, lanczos_params_filename);
    }
    vmatrix = new VMatrix(H);
}

CPT::~CPT()
{
    delete H;
    delete greens;
    delete vmatrix;
}

/*
 * calculate the Greens function <<a^dag,a>>(omega),(k)
 */
std::complex<FLTYPE> CPT::G(std::complex<FLTYPE> omega, Vector3d k)
{
    const std::complex<FLTYPE> i(0, 1);
    MatrixXcd G = greens->get_nambu_greens(omega);
    G = (G.inverse() - vmatrix->get_nambu_V(k)).inverse();
    
    //Fourier transform
    std::complex<FLTYPE> ret(0,0);
    for(int m = 0; m < H->numsites; ++m)
        for(int n = 0; n < H->numsites; ++n)
            ret += G(m,n) * std::exp(-i*k.dot(H->graph.site_pos[m%H->numsites].pos-H->graph.site_pos[n%H->numsites].pos));
    return ret/((FLTYPE)(H->numsites));
}

/*
 * calculate anomalous the Greens function <<a,a>>(omega),(k)
 */
std::complex<FLTYPE> CPT::F(std::complex<FLTYPE> omega, Vector3d k)
{
    const std::complex<FLTYPE> i(0, 1);
    MatrixXcd G = greens->get_nambu_greens(omega);
    G = (G.inverse() - vmatrix->get_nambu_V(k)).inverse();
    
    //Fourier transform
    std::complex<FLTYPE> ret(0,0);
    for(int m = H->numsites; m < H->numsites*2; ++m)
        for(int n = 0; n < H->numsites; ++n)
            ret += G(m,n) * std::exp(-i*k.dot(H->graph.site_pos[m%H->numsites].pos-H->graph.site_pos[n%H->numsites].pos));
    return ret/((FLTYPE)(H->numsites));
}

/*
 * calculate the full response function (<<a^dag,a>>+<<a,a>>+<<a,a^dag>>+<<a^dag,a^dag>>)(omega),(k)
 */
std::complex<FLTYPE> CPT::nambu_G(std::complex<FLTYPE> omega, Vector3d k)
{
    const std::complex<FLTYPE> i(0, 1);
    MatrixXcd G = greens->get_nambu_greens(omega);
    G = (G.inverse() - vmatrix->get_nambu_V(k)).inverse();
    
    //Fourier transform
    std::complex<FLTYPE> ret(0,0);
    for(int m = 0; m < H->numsites*2; ++m)
        for(int n = 0; n < H->numsites*2; ++n)
            ret += G(m,n) * std::exp(-i*k.dot(H->graph.site_pos[m%H->numsites].pos-H->graph.site_pos[n%H->numsites].pos));
    return ret/((FLTYPE)(H->numsites*2));
}

/*
 * calculate the CPT Greens function without periodization
 */
MatrixXcd CPT::cluster_nambu_G(std::complex<FLTYPE> omega, Vector3d k)
{
    MatrixXcd G = greens->get_nambu_greens(omega);
    G = G*((MatrixXcd::Identity(H->numsites*2,H->numsites*2) - vmatrix->get_nambu_V(k) * G).inverse());
    return G;
}

/*
 * calculate the nambu CPT Greens function without periodization and with a precalculated cluster Greens function at omega
 */
MatrixXcd CPT::cluster_nambu_G_precalc(std::complex<FLTYPE> omega, Vector3d k, MatrixXcd G)
{
    return G*((MatrixXcd::Identity(H->numsites*2,H->numsites*2) - vmatrix->get_nambu_V(k) * G).inverse());
}

/*
 * calculate the CPT Greens function without periodization
 */
MatrixXcd CPT::cluster_G(std::complex<FLTYPE> omega, Vector3d k)
{
    MatrixXcd G = greens->get_greens(omega);
    G = G*((MatrixXcd::Identity(H->numsites,H->numsites) - vmatrix->getV(k) * G).inverse());
    return G;
}

/*
 * calculate the CPT Greens function without periodization and with a precalculated cluster Greens function
 */
MatrixXcd CPT::cluster_G_precalc(std::complex<FLTYPE> omega, Vector3d k, MatrixXcd G)
{
    return G*((MatrixXcd::Identity(H->numsites,H->numsites) - vmatrix->getV(k) * G).inverse());
}

/*
 * perform a nambu Greens function sweep from (k_start) to (k_end) with (k_num_steps) number of steps. 
 * resulting CPT Greens function values are stored in an array and then returned.
 * self energy periodization can be applied by (with_self_energy_periodization)
 */
std::vector<std::complex<FLTYPE>> CPT::nambu_G_k_sweep(std::complex<FLTYPE> omega, Vector3d k_start, Vector3d k_end, int k_num_steps, bool with_self_energy_periodization)
{
    if(!with_self_energy_periodization)
    {
        const std::complex<FLTYPE> i(0, 1);
        std::vector<std::complex<FLTYPE>> ret;
        MatrixXcd G = greens->get_nambu_greens(omega);
        for(int j = 0; j < k_num_steps; ++j)
        {
            Vector3d k = (k_end * (FLTYPE)j + k_start*((FLTYPE)k_num_steps - (FLTYPE)j))/((FLTYPE)k_num_steps);
            
            MatrixXcd tG = G*((MatrixXcd::Identity(H->numsites*2,H->numsites*2) - vmatrix->get_nambu_V(k)*G).inverse());
            //fourier transform
            std::complex<FLTYPE> lret(0,0);
            for(int m = 0; m < H->numsites*2; ++m)
                for(int n = 0; n < H->numsites*2; ++n)
                {
                    lret += tG(m,n) * std::exp(-i*k.dot(H->graph.site_pos[m%H->numsites].pos-H->graph.site_pos[n%H->numsites].pos));
                }
            lret /= H->numsites*2.0;
            ret.push_back(lret);
        }
        return ret;
    }
    else
    {
        std::vector<std::complex<FLTYPE>> ret;
        MatrixXcd G = greens->get_nambu_greens(omega);
        G = G.inverse();
        for(int j = 0; j < k_num_steps; ++j)
        {
            Vector3d k = (k_end * (FLTYPE)j + k_start*((FLTYPE)k_num_steps - (FLTYPE)j))/((FLTYPE)k_num_steps);
            //calculate full response <<a^dag,a>>+<<a,a>>+<<a,a^dag>>+<<a^dag,a^dag>> with self energy periodization
            ret.push_back(partial_ft_G_with_self_energy_precalc(omega, k, G)(0,0)+
                partial_ft_G_with_self_energy_precalc(omega, k, G)(1,0)+
                partial_ft_G_with_self_energy_precalc(omega, k, G)(0,1)+
                partial_ft_G_with_self_energy_precalc(omega, k, G)(1,1));
        }
        return ret;
    }
}

/*
 * perform a nambu Greens function sweep from (k_start) to (k_end) with (k_num_steps) number of steps. 
 * resulting CPT Greens function values are stored in an array and then returned.
 * self energy periodization can be applied by (with_self_energy_periodization)
 */
std::vector<std::complex<FLTYPE>> CPT::G_k_sweep(std::complex<FLTYPE> omega, Vector3d k_start, Vector3d k_end, int k_num_steps, bool with_self_energy_periodization)
{
    if(!with_self_energy_periodization)
    {
        const std::complex<FLTYPE> i(0, 1);
        std::vector<std::complex<FLTYPE>> ret;
        MatrixXcd G = greens->get_greens(omega);
        for(int j = 0; j < k_num_steps; ++j)
        {
            Vector3d k = (k_end * (FLTYPE)j + k_start*((FLTYPE)k_num_steps - (FLTYPE)j))/((FLTYPE)k_num_steps);
            
            MatrixXcd tG = G*((MatrixXcd::Identity(H->numsites,H->numsites) - vmatrix->getV(k)*G).inverse());
            
            //fourier transform
            std::complex<FLTYPE> lret(0,0);
            for(int m = 0; m < H->numsites; ++m)
                for(int n = 0; n < H->numsites; ++n)
                    lret += tG(m,n) * std::exp(-i*k.dot(H->graph.site_pos[m%H->numsites].pos-H->graph.site_pos[n%H->numsites].pos));
            lret /= H->numsites;
            ret.push_back(lret);
        }
        return ret;
    }
    else
    {
        std::vector<std::complex<FLTYPE>> ret;
        MatrixXcd G = greens->get_nambu_greens(omega);
        G = G.inverse();
        for(int j = 0; j < k_num_steps; ++j)
        {
            Vector3d k = (k_end * (FLTYPE)j + k_start*((FLTYPE)k_num_steps - (FLTYPE)j))/((FLTYPE)k_num_steps);
            ret.push_back(partial_ft_G_with_self_energy_precalc(omega, k, G)(0,0));
        }
        return ret;
    }
}

/*
 * forwards a recalulation of the Greens function parameters
 */
void CPT::recalculate_greens(bool no_interaction)
{
    greens->prepare_greens(no_interaction);
}

/*
 * returns inverse of the noninteracting Greens function of the entire lattice
 */
MatrixXcd CPT::cluster_noninteracting_inverse_G(std::complex<FLTYPE> omega, Vector3d k)
{
    MatrixXcd Ginv = greens->get_noninteracting_inverse_greens(omega, greens->get_H()->graph.is_nambu, false);
    if(greens->get_H()->graph.is_nambu)
        return Ginv - vmatrix->get_nambu_V(k);
    else
        return Ginv - vmatrix->getV(k);
}

/*
 * returns the noninteracting Greens function of the entire lattice
 */
MatrixXcd CPT::cluster_noninteracting_G(std::complex<FLTYPE> omega, Vector3d k)
{
    MatrixXcd Ginv = greens->get_noninteracting_inverse_greens(omega, greens->get_H()->graph.is_nambu, false);
    if(greens->get_H()->graph.is_nambu)
        return (Ginv - vmatrix->get_nambu_V(k)).inverse();
    else
        return (Ginv - vmatrix->getV(k)).inverse();
}

/*
 * performs a noninteracting Greens function sweep from (k_start) to (k_end) with (k_num_steps) number of steps. 
 * resulting CPT Greens function values are stored in an array and then returned.
 */
std::vector<std::complex<FLTYPE>> CPT::noninteracting_G_k_sweep(std::complex<FLTYPE> omega, Vector3d k_start, Vector3d k_end, int k_num_steps)
{
    const std::complex<FLTYPE> i(0, 1);
    std::vector<std::complex<FLTYPE>> ret;
    for(int j = 0; j < k_num_steps; ++j)
    {
        Vector3d k = (k_end * (FLTYPE)j + k_start*((FLTYPE)k_num_steps - (FLTYPE)j))/((FLTYPE)k_num_steps);
        
        MatrixXcd tG = cluster_noninteracting_G(omega, k);
        
        //fourier transform
        std::complex<FLTYPE> lret(0,0);
        for(int m = 0; m < H->numsites; ++m)
            for(int n = 0; n < H->numsites; ++n)
            {
                lret += tG(m,n) * std::exp(-i*k.dot(H->graph.site_pos[m%H->numsites].pos-H->graph.site_pos[n%H->numsites].pos));
            }
        lret /= H->numsites;
        ret.push_back(lret);
    }
    return ret;
}

/*
 * returns the periodized cluster self energy
 */
std::complex<FLTYPE> CPT::self_energy(std::complex<FLTYPE> omega, Vector3d k)
{
    const std::complex<FLTYPE> i(0, 1);
    MatrixXcd Sigma;
    //the falses are for the exclusion of one body terms caused by mean field decoupling
    if(greens->get_H()->graph.is_nambu)
        Sigma = greens->get_noninteracting_inverse_greens(omega, greens->get_H()->graph.is_nambu, false) - greens->get_nambu_greens(omega).inverse();
    else
        Sigma = greens->get_noninteracting_inverse_greens(omega, greens->get_H()->graph.is_nambu, false) - greens->get_greens(omega).inverse();
    std::complex<FLTYPE> lret(0,0);
    for(int m = 0; m < (greens->get_H()->graph.is_nambu ? H->numsites * 2 : H->numsites); ++m)
        for(int n = 0; n < (greens->get_H()->graph.is_nambu ? H->numsites * 2 : H->numsites); ++n)
            lret += Sigma(m,n) * std::exp(-i*k.dot(H->graph.site_pos[m%H->numsites].pos-H->graph.site_pos[n%H->numsites].pos));
    lret /= (greens->get_H()->graph.is_nambu ? H->numsites * 2 : H->numsites);
    return lret;
}

/*
 * returns the CPT Greens function with full self energy periodization
 */
std::complex<FLTYPE> CPT::G_with_self_energy(std::complex<FLTYPE> omega, Vector3d k)
{
    const std::complex<FLTYPE> i(0, 1);
    const double pi = std::acos(-1);
    MatrixXcd G0 = cluster_noninteracting_G(omega, k);
    std::complex<FLTYPE> lret(0,0);
    for(int m = 0; m < (greens->get_H()->graph.is_nambu ? H->numsites * 2 : H->numsites); ++m)
        for(int n = 0; n < (greens->get_H()->graph.is_nambu ? H->numsites * 2 : H->numsites); ++n)
        {
            lret += G0(m,n) * std::exp(-i*k.dot(H->graph.site_pos[m%H->numsites].pos-H->graph.site_pos[n%H->numsites].pos));
        }
    lret /= (greens->get_H()->graph.is_nambu ? H->numsites * 2 : H->numsites);
    return 1.0/((1.0/lret) - self_energy(omega, k));
}

/*
 * returns the CPT Greens function with full self energy periodization using a precalculated cluster Greens function
 */
std::complex<FLTYPE> CPT::self_energy_precalc(std::complex<FLTYPE> omega, Vector3d k, MatrixXcd Ginv)
{
    const std::complex<FLTYPE> i(0, 1);
    MatrixXcd Sigma;
    Sigma = greens->get_noninteracting_inverse_greens(omega, greens->get_H()->graph.is_nambu, false) - Ginv;
    std::complex<FLTYPE> lret(0,0);
    for(int m = 0; m < (greens->get_H()->graph.is_nambu ? H->numsites * 2 : H->numsites); ++m)
        for(int n = 0; n < (greens->get_H()->graph.is_nambu ? H->numsites * 2 : H->numsites); ++n)
        {
            lret += Sigma(m,n) * std::exp(-i*k.dot(H->graph.site_pos[m%H->numsites].pos-H->graph.site_pos[n%H->numsites].pos));
        }
    lret /= (greens->get_H()->graph.is_nambu ? H->numsites * 2 : H->numsites);
    return lret;
}

/*
 * returns the CPT Greens function with full self energy periodization using a precalculated cluster Greens function
 */
std::complex<FLTYPE> CPT::G_with_self_energy_precalc(std::complex<FLTYPE> omega, Vector3d k, MatrixXcd Ginv)
{
    const std::complex<FLTYPE> i(0, 1);
    MatrixXcd G0 = cluster_noninteracting_G(omega, k);
    std::complex<FLTYPE> lret(0,0);
    for(int m = 0; m < (greens->get_H()->graph.is_nambu ? H->numsites * 2 : H->numsites); ++m)
        for(int n = 0; n < (greens->get_H()->graph.is_nambu ? H->numsites * 2 : H->numsites); ++n)
        {
            lret += G0(m,n) * std::exp(-i*k.dot(H->graph.site_pos[m%H->numsites].pos-H->graph.site_pos[n%H->numsites].pos));
        }
    lret /= (greens->get_H()->graph.is_nambu ? H->numsites * 2 : H->numsites);
    return 1.0/((1.0/lret) - self_energy_precalc(omega, k, Ginv));
}

/*
 * returns the cluster self energy with partial self energy periodization using a precalculated cluster Greens function
 */
MatrixXcd CPT::partial_ft_self_energy_precalc(std::complex<FLTYPE> omega, Vector3d k, MatrixXcd Ginv)
{
    if(!greens->get_H()->graph.is_nambu)
    {
        MatrixXcd lret = MatrixXcd::Zero(1, 1);
        lret(0,0) = self_energy_precalc(omega, k, Ginv);
        return lret;
    }
    const std::complex<FLTYPE> i(0, 1);
    MatrixXcd Sigma;
    //the false is for the exclusion of one body terms caused by mean field decoupling
    Sigma = greens->get_noninteracting_inverse_greens(omega, greens->get_H()->graph.is_nambu, false) - Ginv;
    const double pi = std::acos(-1);
    MatrixXcd lret = MatrixXcd::Zero(2, 2);
    for(int m = 0; m < H->numsites * 2; ++m)
        for(int n = 0; n < H->numsites * 2; ++n)
        {
            int phaseA = 1;//(m / H->numsites) * 2 - 1;
            int phaseB = 1;//(n / H->numsites) * 2 - 1;
            lret(m/H->numsites,n/H->numsites) += Sigma(m,n) * std::exp(-i*k.dot(phaseA*H->graph.site_pos[m%H->numsites].pos-phaseB*H->graph.site_pos[n%H->numsites].pos));
        }
    lret /= H->numsites;
    return lret;
}

/*
 * returns the CPT Greens function with partial self energy periodization using a precalculated cluster Greens function
 */
MatrixXcd CPT::partial_ft_G_with_self_energy_precalc(std::complex<FLTYPE> omega, Vector3d k, MatrixXcd Ginv)
{
    if(!greens->get_H()->graph.is_nambu)
    {
        MatrixXcd lret = MatrixXcd::Zero(1, 1);
        lret(0,0) = G_with_self_energy_precalc(omega, k, Ginv);
        return lret;
    }
    const std::complex<FLTYPE> i(0, 1);
    MatrixXcd G0 = cluster_noninteracting_G(omega, k);
    MatrixXcd lret = MatrixXcd::Zero(2,2);
    for(int m = 0; m < H->numsites * 2; ++m)
        for(int n = 0; n < H->numsites * 2; ++n)
        {
            int phaseA = 1;//(m / H->numsites) * 2 - 1;
            int phaseB = 1;//(n / H->numsites) * 2 - 1;
            lret(m/H->numsites,n/H->numsites) += G0(m,n) * std::exp(-i*k.dot(phaseA*H->graph.site_pos[m%H->numsites].pos-phaseB*H->graph.site_pos[n%H->numsites].pos));
        }
    lret /= H->numsites;
    return (lret.inverse() - partial_ft_self_energy_precalc(omega, k, Ginv)).inverse();
}
