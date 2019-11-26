/*
This code is part of DIETA

Authored by Kirill Alpin
*/

#include <iostream>
#include <fstream>
#include <bitset>
#include "hamiltonian.h"


/*
 * precalculate binomial coefficients
 */
int binomialCoeff_c(int n, int k)
{
    if(k > n)
        return 0;
    int res = 1;
 
    // Since C(n, k) = C(n, n-k)
    if ( k > n - k )
        k = n - k;
 
    // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
    for (int i = 0; i < k; ++i)
    {
        res *= (n - i);
        res /= (i + 1);
    }
 
    return res;
}

/*
 * used in calculating occupation numbers
 */
int hibit(unsigned int n) {
    n |= (n >>  1);
    n |= (n >>  2);
    n |= (n >>  4);
    n |= (n >>  8);
    n |= (n >> 16);
    return ((n - (n >> 1)) << 1) - 1;
}

/*
 * used in calculating occupation numbers
 */
static inline uint32_t log2(const int xp) {
    uint32_t x = (uint32_t)xp;
    uint32_t y;
  asm ( "\tbsr %1, %0\n"
      : "=r"(y)
      : "r" (x)
  );
    return y+1;
}

/*
 * get precalculated binomial coefficients
 */
unsigned int binomialCoeff(int m, int n, unsigned int* precalc, int particle_number)
{
    return precalc[n + m * (particle_number + 2)];
}

/*
 * initialize Hamiltonian
 */
Hamiltonian::Hamiltonian(std::string graph_filename, std::string params_filename) : graph(graph_filename, params_filename)
{
    //initialize variables
    binom_coeff = 0;
    occup_list = 0;
    occup_list_part = 0;
    occup_list_hole = 0;
    occup_lists = 0;
    dims = 0;
    ph_transform = graph.ph_transform;
    //get number of sites
    numsites = graph.site_pos.size();
    //set vectors to zero
    densities = VectorXd::Zero(numsites);
    h_densities = VectorXd::Zero(numsites);
    current_densities = MatrixXcd::Zero(numsites, numsites);
    h_current_densities = MatrixXcd::Zero(numsites, numsites);
    anom_current_densities = MatrixXcd::Zero(numsites, numsites);
    single_operator_densities = VectorXcd::Zero(numsites);
    single_operator_densities_set = false;
    std::cout << "Number of sites:\t" << numsites << std::endl;
    space_dims = graph.translations.size();
    //get Hilbert space trunctation parameters
    particle_number = graph.particle_number;
    is_particle_conserving = graph.is_particle_conserving;
    is_parity_conserving = graph.is_parity_conserving;
    hardcore = graph.hardcore_current;
    //bath site stuff
    if(graph.num_bath > 0)
        std::cout << "Number of bath sites:\t" << graph.num_bath << std::endl;
    bath_mu = VectorXd::Zero(graph.num_bath);
    bath_hopping = MatrixXcd::Zero(numsites, graph.num_bath);
    
    //set meanfields using the parameters in the params file
    for(int i = 0; i < graph.param_names.size(); ++i)
    {
        //set particle density mean field parameters <n>
        if(graph.param_names[i][0] == 'd')
        {
            std::string sub = graph.param_names[i].substr(1);
            int site = std::stoi(sub);
            if(site < numsites)
            {
                densities[site] = std::real(graph.params[i]);
            }
        }
        //set single operator expectation values
        if(graph.param_names[i][0] == 'f')
        {
            std::string sub = graph.param_names[i].substr(1);
            int site = std::stoi(sub);
            if(site < numsites)
            {
                single_operator_densities[site] = graph.params[i];
                single_operator_densities_set = true;
            }
        }
        if(graph.param_names[i][0] == 'm')
        {
            if(graph.param_names[i][1] == 'a' && std::isdigit(graph.param_names[i][2]))
            {
                //set <a_i a_j> mean field parameters
                //anom current densities
                std::string sub = graph.param_names[i].substr(2);
                std::string delimiter = "_";
                int siteA = std::stoi(sub.substr(0, sub.find(delimiter)));
                sub.erase(0, sub.find(delimiter) + delimiter.length());
                int siteB = std::stoi(sub.substr(0, sub.find(delimiter)));
                if(siteA < numsites && siteB < numsites)
                {
                    anom_current_densities(siteA, siteB) = graph.params[i];
                    anom_current_densities(siteB, siteA) = std::conj(graph.params[i]);
                    current_set = true;
                }
            }
            else if(std::isdigit(graph.param_names[i][1]))
            {
                //set <a^dag_i a_j> mean field parameters
                //current densities
                std::string sub = graph.param_names[i].substr(1);
                std::string delimiter = "_";
                int siteA = std::stoi(sub.substr(0, sub.find(delimiter)));
                sub.erase(0, sub.find(delimiter) + delimiter.length());
                int siteB = std::stoi(sub.substr(0, sub.find(delimiter)));
                if(siteA < numsites && siteB < numsites)
                {
                    current_densities(siteA, siteB) = graph.params[i];
                    current_densities(siteB, siteA) = std::conj(graph.params[i]);
                    current_set = true;
                }
            }
        }
        if(graph.param_names[i][0] == 'b')
        {
            //set bath site parameters
            if(graph.param_names[i][1] == 'm' && std::isdigit(graph.param_names[i][2]))
            {
                //set bath site chemical potential
                std::string sub = graph.param_names[i].substr(2);
                int site = std::stoi(sub);
                if(site < graph.num_bath)
                    bath_mu[site] = std::real(graph.params[i]);
            }
            else if(std::isdigit(graph.param_names[i][1]))
            {
                //set bath site hybridization
                std::string sub = graph.param_names[i].substr(1);
                std::string delimiter = "_";
                int siteA = std::stoi(sub.substr(0, sub.find(delimiter)));
                sub.erase(0, sub.find(delimiter) + delimiter.length());
                int siteB = std::stoi(sub.substr(0, sub.find(delimiter)));
                if(siteA < numsites && siteB < graph.num_bath)
                    bath_hopping(siteA, siteB) = graph.params[i];
            }
        }
    }
    //make parameters particle hole symmetric
    h_densities = densities;
    h_current_densities = current_densities;
    
    if(numsites*(graph.max_particle == 1 ? 1 : (graph.max_particle < 4 ? 2 : log2(graph.max_particle))) > sizeof(unsigned int)*8)
    {
        std::stringstream ss;
        ss << "Number of sites exceeds the occupation number representation size.";
        throw std::invalid_argument(ss.str());
    }
    if(is_parity_conserving)
    {
        //set up Hilbert space trunctation based on parity conservation
        
        //due to parity symmetry, the groundstate particle number defaults to numsites/2
        particle_number = numsites / 2;
        if(graph.max_particle != 1)
        {
            std::stringstream ss;
            ss << "Particle number conservation not implemented for " << graph.max_particle << " particles per site.";
            throw std::invalid_argument(ss.str());
        }
        
        //precalculate binomial coefficients
        binom_coeff = new unsigned int[(numsites+1)*(numsites+1)];
        for(int m = 0; m < numsites + 1; ++m)
            for(int n = 0; n < numsites + 1; ++n)
                binom_coeff[n + m * (numsites + 1)] = binomialCoeff_c(m, n);
        
        //precalculate the occupation list
        occup_lists = new unsigned int*[numsites+1];
        dims = new unsigned int[numsites+1];
        for(int i = 0; i < numsites+1; ++i)
        {
            dims[i] = binomialCoeff_c(numsites, i);
            compute_occup_lookup_table(occup_lists[i], i, dims[i]);
        }
        //calculate the Hilbert dimension of the subsectors
        unsigned int dim_even = 0;
        unsigned int dim_odd = 0;
        for(int i = 0; i < numsites+1; i+=2)
            dim_even += dims[i];
        for(int i = 1; i < numsites+1; i+=2)
            dim_odd += dims[i];
        
        if(numsites % 2 == 0)
        {
            dim = dim_even;
            dim_part = dim_odd;
            dim_hole = dim_odd;
        }
        else
        {
            dim = dim_odd;
            dim_part = dim_even;
            dim_hole = dim_even;
        }
    }
    else if(is_particle_conserving)
    {
        //set up Hilbert space trunctation based on particle number conservation
        if(graph.max_particle != 1)
        {
            std::stringstream ss;
            ss << "Particle number conservation not implemented for " << graph.max_particle << " particles per site.";
            throw std::invalid_argument(ss.str());
        }
        //precompute binomial coefficients
        binom_coeff = new unsigned int[(numsites+1)*(particle_number+2)];
        for(int m = 0; m < numsites + 1; ++m)
            for(int n = 0; n < particle_number + 2; ++n)
                binom_coeff[n + m * (particle_number + 2)] = binomialCoeff_c(m, n);
        std::cout << "Number of particles:\t" << particle_number << std::endl;
        //set dimensions
        dim = binomialCoeff_c(numsites, particle_number);
        dim_part = particle_number + 1 > numsites ? 1 : binomialCoeff_c(numsites, particle_number + 1);
        dim_hole = particle_number - 1 < 0 ? 1 : binomialCoeff_c(numsites, particle_number - 1);
        
        //precompute occupation numbers for N, N+1 and N-1 particle sector
        compute_occup_lookup_table(occup_list, particle_number, dim);
        compute_occup_lookup_table(occup_list_part, particle_number + 1 > numsites ? numsites : particle_number + 1, dim_part);
        compute_occup_lookup_table(occup_list_hole, particle_number - 1 < 0 ? 0 : particle_number - 1, dim_hole);
    }
    else
    {
        //set up dimensions for the case of no Hilbert space trunctation
        dim = get_dimension();
        dim_part = dim;
        dim_hole = dim;
    }
    //check if coefficients are complex
    is_complex = graph.is_complex;
    //precompute sqrt coefficients for boson operators
    sqrt_coeff = new FLTYPE[graph.max_particle+1];
    for(int i = 0; i < graph.max_particle+1; ++i)
        sqrt_coeff[i] = std::sqrt((FLTYPE)i);
}

Hamiltonian::~Hamiltonian()
{
    //deallocate all stuff
    if(binom_coeff != 0)
    {
        delete[] binom_coeff;
        binom_coeff = 0;
    }
    if(occup_list != 0)
    {
        delete[] occup_list;
        occup_list = 0;
    }
    if(occup_list_hole != 0)
    {
        delete[] occup_list_hole;
        occup_list_hole = 0;
    }
    if(occup_list_part != 0)
    {
        delete[] occup_list_part;
        occup_list_part = 0;
    }
    if(occup_lists != 0)
    {
        for(int i = 0; i < numsites; ++i)
            delete[] occup_lists[i];
        delete[] occup_lists;
        occup_lists = 0;
    }
    if(dims != 0)
    {
        delete[] dims;
        dims = 0;
    }
    if(sqrt_coeff != 0)
    {
        delete[] sqrt_coeff;
        sqrt_coeff = 0;
    }
}

/*
 * precomputes the occupation lookup table in case of Hilbert space trunctation
 */
void Hamiltonian::compute_occup_lookup_table(unsigned int* &occup_lookup, unsigned int particle_number_now, unsigned int dim_now)
{
    std::vector<unsigned int> now_occup;
    for(int m = 0; m < numsites - particle_number_now; ++m)
        now_occup.push_back(0);
    for(int n = 0; n < particle_number_now; ++n)
        now_occup.push_back(1);
    
    std::sort(now_occup.begin(), now_occup.end());
    occup_lookup = new unsigned int[dim_now];
    for(int i = 0; i < dim_now; ++i)
    {
        occup_lookup[i] = 0;
        for(int j = 0; j < now_occup.size(); ++j)
            occup_lookup[i] |= now_occup[j] << j;
        std::next_permutation(now_occup.begin(), now_occup.end());
    }
}

/*
 * returns the total Hibert space dimension
 */
unsigned int Hamiltonian::get_dimension()
{
    if(graph.max_particle == 1)
        return 1 << numsites;
    else
    {
        unsigned int ret = 1;
        for(int i = 0; i < numsites; ++i)
            ret *= graph.max_particle + 1;
        return ret;
    }
}

/*
 * converts an occupation number vector (occup) to the Fock state index in the state vector
 */
unsigned int Hamiltonian::occup_to_index(unsigned int occup)
{
    if(is_parity_conserving)
    {
        int ret = 0;
        std::bitset<sizeof(unsigned int)*8> bito(occup);
        int remaining_particles = bito.count();
        int remaining_sites = numsites;
        for(int i = remaining_particles % 2; i < remaining_particles; i+=2)
            ret += dims[i];
        for(int i = 0; i < numsites - 1; ++i)
        {
            int action_i = (occup >> i) & 1;
            if(action_i == 1)
            {
                remaining_sites--;
                ret += binomialCoeff(remaining_sites, remaining_particles, binom_coeff, numsites - 1);
                remaining_particles--;
            }
            else
                remaining_sites--;
        }
        return ret;
    }
    else if(is_particle_conserving)
    {
        int ret = 0;
        int remaining_particles = current_particle_space;
        int remaining_sites = numsites;
        for(int i = 0; i < numsites - 1; ++i)
        {
            int action_i = (occup >> i) & 1;
            if(action_i == 1)
            {
                remaining_sites--;
                ret += binomialCoeff(remaining_sites, remaining_particles, binom_coeff, particle_number);
                remaining_particles--;
            }
            else
                remaining_sites--;
        }
        return ret;
    }
    else
    {
        if(graph.max_particle == 1)
            return occup;
        else
        {
            unsigned int ret = 0;
            unsigned int base = 1;
            for(int i = 0; i < numsites; ++i)
            {
                if(graph.max_particle < 4)
                {
                    unsigned int soccup = (occup >> (i*2)) & 3;
                    ret += soccup * base;
                    base *= graph.max_particle + 1;
                }
                else
                {
                    unsigned int soccup = (occup >> (i*log2(graph.max_particle))) & hibit(graph.max_particle);
                    ret += soccup * base;
                    base *= graph.max_particle + 1;
                }
            }
            return ret;
        }
    }
}

/*
 * converts a Fock state index to a occupation number vector (stored in one unsigned int)
 */
unsigned int Hamiltonian::index_to_occup(unsigned int index)
{
    unsigned int occup;
    if(is_parity_conserving)
    {
        unsigned int d = 0;
        for(int i = current_particle_space % 2; i < numsites + 1; i+=2)
        {
            if(d <= index && d + dims[i] > index)
            {
                occup = occup_lists[i][index - d];
                break;
            }
            else
                d += dims[i];
        }
    }
    else if(is_particle_conserving)
    {
        if(current_particle_space == particle_number)
            occup = occup_list[index];
        else if(current_particle_space == particle_number + 1)
            occup = occup_list_part[index];
        else
            occup = occup_list_hole[index];
    }
    else
    {
        if(graph.max_particle == 1)
            occup = index;
        else
        {
            unsigned int ret = 0;
            unsigned int work = index;
            for(int i = 0; i < numsites; ++i)
            {
                if(graph.max_particle < 4)
                {
                    ret |= (work % (graph.max_particle + 1)) << (i*2);
                    work = work / (graph.max_particle + 1);
                }
                else
                {
                    ret |= (work % (graph.max_particle + 1)) << (i*log2(graph.max_particle));
                    work = work / (graph.max_particle + 1);
                }
            }
            occup = ret;
        }
    }
    return occup;
}

/*
 * returns the occupation number at site index (site) given a occupation number vector (occup)
 */
unsigned int Hamiltonian::get_occup(unsigned int occup, unsigned int site)
{
    if(graph.max_particle == 1)
        return (occup >> site) & 1;
    else
        if(graph.max_particle < 4)
            return (occup >> (site*2)) & 3;
        else
            return (occup >> (site*log2(graph.max_particle))) & hibit(graph.max_particle);
}

/*
 * computes a^dag_site a_site|occup>. it is possible to use an automatic Holstein Primakoff expansion for bosons
 */
FLTYPE Hamiltonian::n_operator(unsigned int occup, unsigned int site)
{
    if(ph_transform)
        occup = ~occup;
    if(graph.max_particle == 1)
        return (FLTYPE)((occup >> site) & 1);
    else
    {
        if(graph.max_particle < 4)
        {
            FLTYPE sq_n = holprim_expansion(((occup >> (site*2)) & 3) - 1);
            return (FLTYPE)((occup >> (site*2)) & 3) * sq_n * sq_n;
        }
        else
        {
            FLTYPE sq_n = holprim_expansion(((occup >> (site*log2(graph.max_particle))) & hibit(graph.max_particle)) - 1);
            return (FLTYPE)((occup >> (site*log2(graph.max_particle))) & hibit(graph.max_particle)) * sq_n * sq_n;
        }
    }
}

/*
 * get Holstein Primakoff expansion of n. essentially equivalent to a taylor series of sqrt(1-n)
 */
FLTYPE Hamiltonian::holprim_expansion(int n)
{
    FLTYPE ret = 1.0;
    int np = n;
    for(int i = 0; i < MAX_HOLPRIM_ORDER; ++i)
    {
        ret -= (FLTYPE)np * HOLPRIM_COEFF[i];
        np *= n;
    }
    return ret;
}

/*
 * adds a particle to the occupation number vector (occup) at site (site). if successful, returns true. 
 * if not, this means the site is already fully occupied. 
 * (coeff) returns a possible coefficient applied onto the state, like the phase due to normal ordering of fermions
 * or sqrt coefficient for bosons
 */
bool Hamiltonian::add_particle(unsigned int &occup, unsigned int site, FLTYPE& coeff)
{
    if(graph.particleType == ParticleType::FERMION)
    {
        //compute normal ordering phase factor of fermions
        int rr = 0;
        for(int i = 0; i < site; ++i)
            rr += get_occup(occup, i);
        coeff *= (1 - ((rr & 1) * 2));
    }
    
    //apply possible particle hole transfromation
    if(ph_transform)
        occup = ~occup;
    //check if site is already fully occupied
    if(get_occup(occup, site) == graph.max_particle)
        return false;
    //check if maximal allowed particle number is 1. indicating either fermions or hardcore bosons
    if(graph.max_particle == 1)
    {
        //add particle at site (site)
        occup = occup | (1 << site);
        //reverse particle hole transfromation
        if(ph_transform)
            occup = ~occup;
        return true;
    }
    else
    {
        if(graph.max_particle < 4)
        {
            //calculate sqrt coefficient for bosons. automatic Holstein Primakoff expansion is also possible
            coeff *= sqrt_coeff[((occup >> (site*2))&3) + 1] * holprim_expansion((occup >> (site*2))&3);
            //apply bit magic to add particle
            occup = (((occup)&(~(3<<(site*2))))) | ((((occup >> (site*2))&3) + 1) << (site*2));
            //reverse particle hole transfromation
            if(ph_transform)
                occup = ~occup;
            return true;
        }
        else
        {
            //calculate sqrt coefficient for bosons. automatic Holstein Primakoff expansion is also possible
            coeff *= sqrt_coeff[((occup >> (site*log2(graph.max_particle)))&hibit(graph.max_particle)) + 1] * holprim_expansion((occup >> (site*log2(graph.max_particle)))&hibit(graph.max_particle));
            //apply more bit magic to add particle
            occup = (((occup)&(~(hibit(graph.max_particle)<<(site*log2(graph.max_particle)))))) | ((((occup >> (site*log2(graph.max_particle)))&hibit(graph.max_particle)) + 1) << (site*log2(graph.max_particle)));
            //reverse particle hole transfromation
            if(ph_transform)
                occup = ~occup;
            return true;
        }
    }
}

/*
 * equivalent to add_particle, only that here a particle is removed from (occup)
 */
bool Hamiltonian::remove_particle(unsigned int &occup, unsigned int site, FLTYPE& coeff)
{
    if(graph.particleType == ParticleType::FERMION)
    {
        //compute normal ordering phase factor of fermions
        int rr = 0;
        for(int i = 0; i < site; ++i)
            rr += get_occup(occup, i);
        coeff *= (1 - ((rr & 1) * 2));
    }
    
    //apply possible particle hole transfromation
    if(ph_transform)
        occup = ~occup;
    //check if site is already at zero occupation
    if(get_occup(occup, site) == 0)
        return false;
    //check if maximal allowed particle number is 1. indicating either fermions or hardcore bosons
    if(graph.max_particle == 1)
    {
        //remove particle at site (site)
        occup = occup & (~(1 << site));
        //reverse particle hole transfromation
        if(ph_transform)
            occup = ~occup;
        return true;
    }
    else
    {
        if(graph.max_particle < 4)
        {
            //calculate sqrt coefficient for bosons. automatic Holstein Primakoff expansion is also possible
            coeff *= sqrt_coeff[((occup >> (site*2))&3)] * holprim_expansion(((occup >> (site*2))&3) - 1);
            //apply bit magic to remove particle
            occup = (((occup)&(~(3<<(site*2))))) | ((((occup >> (site*2))&3) - 1) << (site*2));
            //reverse particle hole transfromation
            if(ph_transform)
                occup = ~occup;
            return true;
        }
        else
        {
            //calculate sqrt coefficient for bosons. automatic Holstein Primakoff expansion is also possible
            coeff *= sqrt_coeff[((occup >> (site*log2(graph.max_particle)))&hibit(graph.max_particle))] * holprim_expansion(((occup >> (site*log2(graph.max_particle)))&hibit(graph.max_particle)) - 1);
            //apply more bit magic to remove particle
            occup = (((occup)&(~(hibit(graph.max_particle)<<(site*log2(graph.max_particle)))))) | ((((occup >> (site*log2(graph.max_particle)))&hibit(graph.max_particle)) - 1) << (site*log2(graph.max_particle)));
            //reverse particle hole transfromation
            if(ph_transform)
                occup = ~occup;
            return true;
        }
    }
}

FARRAY Hamiltonian::single_particle_state(unsigned int site)
{
    //TODO legacy code. not compatible with complex hamiltonian
    FARRAY ret = FARRAY::Zero(dim);
    unsigned int occup = 0;
    FLTYPE coeff = 1.0;
    add_particle(occup, site, coeff);
    ret[occup_to_index(occup)] = 1.0;
    return ret;
}

/*
 * check sanity of Hibert space dimension
 */
void Hamiltonian::check_hilbert_dim(unsigned int dim_now)
{
    if(is_particle_conserving)
    {
        if(dim_now != dim_hole && dim_now != dim_part && dim_now != dim)
        {
            std::stringstream ss;
            ss << "Invalid hilbert space dimension.";
            throw std::invalid_argument(ss.str());
        }
        if(current_particle_space != particle_number 
            && current_particle_space != particle_number + 1
            && current_particle_space != particle_number - 1)
        {
            std::stringstream ss;
            ss << "Invalid subspace.";
            throw std::invalid_argument(ss.str());
        }
    }
    else
    {
        if(dim_now != dim)
        {
            std::stringstream ss;
            ss << "Invalid hilbert space dimension.";
            throw std::invalid_argument(ss.str());
        }
    }
}

/*
 * set current particle space on which the code operates
 */
void Hamiltonian::set_particle_space(unsigned int p)
{
    if(is_particle_conserving || is_parity_conserving)
        current_particle_space = p;
}

/*
 * returns a^dag_(site_num) a_(site_num)|s>
 */
template<typename T>
T Hamiltonian::nop_site(T s, int site_num)
{
    T ret = T::Zero(s.size());
    check_hilbert_dim(s.size());
    if(is_particle_conserving || is_parity_conserving)
    {
        std::stringstream ss;
        ss << "Not implemented.";
        throw std::invalid_argument(ss.str());
    }
    for(unsigned int i = 0; i < s.size(); ++i)
    {
        unsigned int occup = index_to_occup(i);
        auto const site = graph.site_pos[site_num];
        switch(site.type)
        {
            case 'A': //chemical potential
                {
                    FLTYPE site_occup = n_operator(occup, site.num);
                    if(ph_transform)
                        site_occup = 1.0 - site_occup;
                    ret[i] += site_occup * s[i];
                }
                break;
            case 'F': //on site term a+adag
                {
                }
                break;
            case 'N': //noop
                break;
            default:
                std::stringstream ss;
                ss << "Site type " << site.type << " not implemented.";
                throw std::invalid_argument(ss.str());
                break;
        }
    }
    return ret;
}

template VectorXd Hamiltonian::nop_site(VectorXd s, int site_num);
template VectorXcd Hamiltonian::nop_site(VectorXcd s, int site_num);

/*
 * used to incorporate complex valued Hamiltonians
 */
std::complex<FLTYPE> Hamiltonian::K(std::complex<FLTYPE> val, Tag<VectorXcd>)
{
    return val;
}

FLTYPE Hamiltonian::K(std::complex<FLTYPE> val, Tag<VectorXd>)
{
    return val.real();
}

/*
 * set fields F_field, mu_field and meanF_field given the parameters of the Hamiltonian
 */
void Hamiltonian::setFields()
{
    F_field = VectorXcd::Zero(numsites);
    mu_field = VectorXd::Zero(numsites);
    meanF_field = VectorXcd::Zero(numsites);
    for(auto const& site: graph.site_pos)
    {
        switch(site.type)
        {
            case 'A': //chemical potential
                {
                    mu_field[site.num] += std::real(P(site.ind, Tag<VectorXcd>()));
                }
                break;
            default:
                break;
        }
    }
    for(auto const& edge: graph.edges)
    {
        switch(edge.type)
        {
            case 'F': //a+a^dag term. edge.b is unused
                {
                    unsigned int siteA = edge.a;
                    F_field[siteA] += P(edge.ind, Tag<VectorXcd>());
                }
                break;
            default:
                break;
        }
    }
    for(auto const& edge: graph.mean_edges)
    {
        switch(edge.type)
        {
            case 'C':
                {
                    unsigned int siteA = edge.a;
                    unsigned int siteB = edge.b;
                    mu_field[siteA] += std::real(P(edge.ind, Tag<VectorXcd>())) * densities[siteB];
                    mu_field[siteB] += std::real(P(edge.ind, Tag<VectorXcd>())) * densities[siteA];
                }
                break;
            case 'D':
                {
                    unsigned int siteA = edge.a;
                    unsigned int siteB = edge.b;
                    meanF_field[siteB] += P(edge.ind, Tag<VectorXcd>()) * densities[siteA];
                }
                break;
            default:
                break;
        }
    }
    F_field += meanF_field;
    mu_field *= -1;
}

/*
 * applies the Hamiltonian onto a state by multiplication
 */
template<typename T>
T Hamiltonian::operator* (const T &s)
{
    //set the returning state vector
    T ret = T::Zero(s.size());
    //check if vector size is compatible with the Hamiltonian
    check_hilbert_dim(s.size());
    if(is_particle_conserving)
    {
        //only a subset of operation is supported for particle number conserving Hamiltonians
        //iterate over all Fock state indices to construct the returning vector
        #pragma omp parallel for
        for(unsigned int i = 0; i < s.size(); ++i)
        {
            unsigned int occup = index_to_occup(i);
            //on site stuff
            //iterate over all sites
            for(auto const& site: graph.site_pos)
            {
                switch(site.type)
                {
                    case 'A': //chemical potential
                        {
                            FLTYPE site_occup = n_operator(occup, site.num);
                            ret[i] += PA(site.ind) * site_occup * s[i];
                        }
                        break;
                    case 'N': //noop
                        break;
                    default:
                        std::stringstream ss;
                        ss << "Site type " << site.type << " not implemented or not particle number conserving.";
                        throw std::invalid_argument(ss.str());
                        break;
                }
            }
            //inter site stuff
            //iterate over all edges connecting sites inside the graph
            for(auto const& edge: graph.edges)
            {
                switch(edge.type)
                {
                    case 'A': //hopping term c^dag*c+c*c^dag
                        {
                            unsigned int siteA = edge.a;
                            unsigned int siteB = edge.b;
                            
                            //forward
                            unsigned int occup_work = occup;
                            bool success = true;
                            FLTYPE coeff = 1.0;
                            
                            success &= remove_particle(occup_work, siteA, coeff);
                            success &= add_particle(occup_work, siteB, coeff);
                            if(success)
                                ret[i] += (coeff * std::conj(PA(edge.ind))) * s[occup_to_index(occup_work)];
                            
                            //backward
                            occup_work = occup;
                            success = true;
                            coeff = 1.0;
                            success &= add_particle(occup_work, siteA, coeff);
                            success &= remove_particle(occup_work, siteB, coeff);
                            if(success)
                                ret[i] += (coeff * PA(edge.ind)) * s[occup_to_index(occup_work)];
                        }
                        break;
                    case 'C': //interaction n_a*n_b
                        if(!no_interaction)
                        {
                            unsigned int siteA = edge.a;
                            unsigned int siteB = edge.b;
                            
                            ret[i] += PA(edge.ind) * n_operator(occup, siteA) * n_operator(occup, siteB) * s[i];
                        }
                        break;
                    default:
                        std::stringstream ss;
                        ss << "Edge type " << edge.type << " not implemented or not particle number conserving.";
                        throw std::invalid_argument(ss.str());
                        break;
                }
            }
        }
    }
    else
    {
        //iterate over all Fock state indices to construct the returning vector
        #pragma omp parallel for
        for(unsigned int i = 0; i < s.size(); ++i)
        {
            unsigned int occup = index_to_occup(i);
            //on site stuff
            //iterate over all sites
            for(auto const& site: graph.site_pos)
            {
                switch(site.type)
                {
                    case 'A': //chemical potential
                        {
                            FLTYPE site_occup = n_operator(occup, site.num);
                            ret[i] += PA(site.ind) * (FLTYPE)site_occup * s[i];
                        }
                        break;
                    case 'N': //noop
                        break;
                    default:
                        std::stringstream ss;
                        ss << "Site type " << site.type << " not implemented.";
                        throw std::invalid_argument(ss.str());
                        break;
                }
            }
            //inter site stuff
            //iterate over all edges connecting sites inside the graph
            for(auto const& edge: graph.edges)
            {
                switch(edge.type)
                {
                    case 'A': //hopping term c^dag*c+c*c^dag
                        {
                            unsigned int siteA = edge.a;
                            unsigned int siteB = edge.b;
                            
                            //forward
                            unsigned int occup_work = occup;
                            bool success = true;
                            FLTYPE coeff = 1.0;
                            success &= remove_particle(occup_work, siteA, coeff);
                            success &= add_particle(occup_work, siteB, coeff);
                            if(success)
                                ret[i] += (coeff * std::conj(PA(edge.ind))) * s[occup_to_index(occup_work)];
                            
                            //backward
                            occup_work = occup;
                            success = true;
                            coeff = 1.0;
                            success &= remove_particle(occup_work, siteB, coeff);
                            success &= add_particle(occup_work, siteA, coeff);
                            if(success)
                                ret[i] += (coeff * PA(edge.ind)) * s[occup_to_index(occup_work)];
                        }
                        break;
                    case 'B': //pair term c*c+c^dag*c^dag
                        {
                            unsigned int siteA = edge.a;
                            unsigned int siteB = edge.b;
                            
                            //forward
                            unsigned int occup_work = occup;
                            bool success = true;
                            FLTYPE coeff = 1.0;
                            success &= add_particle(occup_work, siteA, coeff);
                            success &= add_particle(occup_work, siteB, coeff);
                            if(success)
                                ret[i] += (coeff * std::conj(PA(edge.ind))) * s[occup_to_index(occup_work)];
                             
                            //backward
                            occup_work = occup;
                            success = true;
                            coeff = 1.0;
                            success &= remove_particle(occup_work, siteB, coeff);
                            success &= remove_particle(occup_work, siteA, coeff);
                            if(success)
                                ret[i] += (coeff * PA(edge.ind)) * s[occup_to_index(occup_work)];
                        }
                        break;
                    case 'C': //interaction n_a*n_b or n_a*(n_a-1)
                        if(!no_interaction)
                        {
                            unsigned int siteA = edge.a;
                            unsigned int siteB = edge.b;
                            if(siteA == siteB)
                                ret[i] += PA(edge.ind) * n_operator(occup, siteA) * (n_operator(occup, siteB) - 1) * s[i];
                            else
                                ret[i] += PA(edge.ind) * n_operator(occup, siteA) * n_operator(occup, siteB) * s[i];
                        }
                        break;
                    case 'D': //interaction n_a*c_b+n_a*c^dag_b
                        if(!no_interaction)
                        {
                            unsigned int siteA = edge.a;
                            unsigned int siteB = edge.b;
                            
                            unsigned int occup_work = occup;
                            FLTYPE coeff = 1.0;
                            bool success = add_particle(occup_work, siteB, coeff);
                            if(success)
                                ret[i] += (coeff * std::conj(PA(edge.ind))) * n_operator(occup, siteA) * s[occup_to_index(occup_work)];
                            
                            occup_work = occup;
                            coeff = 1.0;
                            success = remove_particle(occup_work, siteB, coeff);
                            if(success)
                                ret[i] += (coeff * PA(edge.ind)) * n_operator(occup, siteA) * s[occup_to_index(occup_work)];
                        }
                        break;
                    case 'F': //a+a^dag term. edge.b is unused
                        {
                            unsigned int siteA = edge.a;
                            
                            unsigned int occup_work = occup;
                            FLTYPE coeff = 1.0;
                            bool success = add_particle(occup_work, siteA, coeff);
                            if(success)
                                ret[i] += (coeff * PA(edge.ind)) * s[occup_to_index(occup_work)];
                            
                            occup_work = occup;
                            coeff = 1.0;
                            success = remove_particle(occup_work, siteA, coeff);
                            if(success)
                                ret[i] += (coeff * std::conj(PA(edge.ind))) * s[occup_to_index(occup_work)];
                        }
                        break;
                    default:
                        std::stringstream ss;
                        ss << "Edge type " << edge.type << " not implemented.";
                        throw std::invalid_argument(ss.str());
                        break;
                }
            }
            
            //additional terms in the Hamiltonian coming from mean field decoupling
            //iterate over all mean field edges in the graph
            for(auto const& edge: graph.mean_edges)
            {
                switch(edge.type)
                {
                    case 'C':
                        if(!no_interaction)
                        {
                            unsigned int siteA = edge.a;
                            unsigned int siteB = edge.b;
                            
                             //interaction decoupling n_a*n_b -> n_a*<n_b>+<n_a>*n_b
                            ret[i] += PA(edge.ind) * n_operator(occup, siteA) * densities[siteB] * s[i]
                                + PA(edge.ind) * densities[siteA] * n_operator(occup, siteB) * s[i]
                                - PA(edge.ind) * densities[siteA] * densities[siteB] * s[i];
                             
                            //if current mean field parameters are set, use them to decouple the interaction
                            if(current_set)
                            {
                                //interaction decoupling a^dag_a*a_a*a^dag_b*a_b -> a^dag_a*a_b*<a^dag_b*a_a>+<a^dag_a*a_b>*a^dag_b*a_a
                                
                                //forward
                                unsigned int occup_work = occup;
                                bool success = true;
                                FLTYPE coeff = 1.0; 
                                if(hardcore)
                                    coeff = -1.0; // due to commutation
                                success &= remove_particle(occup_work, siteA, coeff);
                                success &= add_particle(occup_work, siteB, coeff);
                                if(success) //<a^dag_a*a_b> as current
                                    ret[i] += KA(current_densities(siteB,siteA) * (coeff * std::conj(PA(edge.ind))) * s[occup_to_index(occup_work)]);
                                
                                //backward
                                occup_work = occup;
                                success = true;
                                coeff = 1.0;
                                if(hardcore)
                                    coeff = -1.0;
                                success &= add_particle(occup_work, siteA, coeff);
                                success &= remove_particle(occup_work, siteB, coeff);
                                if(success) //<a^dag_b*a_a> as current
                                    ret[i] += KA(current_densities(siteA,siteB)*(coeff * PA(edge.ind)) * s[occup_to_index(occup_work)]);
                                coeff = 1.0;
                                if(hardcore)
                                    coeff = -1.0;
                                ret[i] -= KA(coeff*current_densities(siteA,siteB)*(std::conj(current_densities(siteA,siteB)) * std::abs(PA(edge.ind))) * s[i]);
                                //interaction decoupling a^dag_a*a_a*a^dag_b*a_b -> a^dag_a*a^dag_b*<a_a*a_b> + <a^dag_a*a^dag_b>*a_a*a_b
                                //forward
                                occup_work = occup;
                                success = true;
                                coeff = 1.0;
                                success &= add_particle(occup_work, siteA, coeff);
                                success &= add_particle(occup_work, siteB, coeff);
                                if(success) // <a_a*a_b> as anom current
                                    ret[i] += KA(std::conj(anom_current_densities(siteA,siteB)) *(coeff * std::conj(PA(edge.ind))) * s[occup_to_index(occup_work)]);
                                
                                //backward
                                occup_work = occup;
                                success = true;
                                coeff = 1.0;
                                success &= remove_particle(occup_work, siteB, coeff);
                                success &= remove_particle(occup_work, siteA, coeff);
                                if(success) // <a^dag_a*a^dag_b> as anom current
                                    ret[i] += KA(anom_current_densities(siteA,siteB) * (coeff * PA(edge.ind)) * s[occup_to_index(occup_work)]);
                                
                                ret[i] -= KA(anom_current_densities(siteA,siteB)*(std::conj(anom_current_densities(siteA,siteB)) * std::abs(PA(edge.ind))) * s[i]);
                            }
                        }
                        break;
                    case 'D': 
                        if(!no_interaction)
                        {
                            //interaction n_a*c_b+n_a*c^dag_b -> <n_a>*c_b+<n_a>*c^dag_b
                            unsigned int siteA = edge.a;
                            unsigned int siteB = edge.b;
                            
                            unsigned int occup_work = occup;
                            FLTYPE coeff = 1.0;
                            bool success = false;
                            success = add_particle(occup_work, siteB, coeff);
                            if(success)
                                ret[i] += (coeff * std::conj(PA(edge.ind))) * densities[siteA] * s[occup_to_index(occup_work)];
                            
                            occup_work = occup;
                            coeff = 1.0;
                            success = remove_particle(occup_work, siteB, coeff);
                            if(success)
                                ret[i] += (coeff * PA(edge.ind)) * densities[siteA] * s[occup_to_index(occup_work)];
                            
                            //if current mean field parameters are set, use them to decouple the interaction
                            if(current_set)
                            {
                                std::complex<FLTYPE> currents = current_densities(siteB,siteA) - anom_current_densities(siteA,siteB);
                                
                                occup_work = occup;
                                coeff = 1.0;
                                success = remove_particle(occup_work, siteA, coeff);
                                if(success)
                                    ret[i] += KA(std::conj(currents) * (coeff * PA(edge.ind)) * s[occup_to_index(occup_work)]);
                                
                                occup_work = occup;
                                coeff = 1.0;
                                success = add_particle(occup_work, siteA, coeff);
                                if(success)
                                    ret[i] += KA(currents * (coeff * std::conj(PA(edge.ind))) * s[occup_to_index(occup_work)]);
                            }
                        }
                        break;
                    default:
                        std::stringstream ss;
                        ss << "Mean-edge type " << edge.type << " not implemented.";
                        throw std::invalid_argument(ss.str());
                        break;
                }
            }
        }
    }
    return ret;
}

template VectorXd Hamiltonian::operator* (const VectorXd &s);
template VectorXcd Hamiltonian::operator* (const VectorXcd &s);

/*
 * used to incorporate complex Hamiltonians
 */
std::complex<FLTYPE> Hamiltonian::P(int i, Tag<VectorXcd>)
{
    return graph.params[i];
}

FLTYPE Hamiltonian::P(int i, Tag<VectorXd>)
{
    return graph.params[i].real();
}

/*
 * adds a particle to a given state (ts) at site (site) with a prefactor (factor)
 */
template<typename T, typename U, typename G>
void Hamiltonian::add_particle_state(G s, int site, T &ts, U factor)
{
    if(is_particle_conserving)
    {
        if(particle_number + 1 > numsites * graph.max_particle)
        {
            ts[0] = 1.0;
            return;
        }
    }
    for(unsigned int i = 0; i < s.size(); ++i)
    {
        current_particle_space = particle_number;
        unsigned int occup = index_to_occup(i);
        FLTYPE coeff = 1.0;
        current_particle_space = particle_number + 1;
        if(add_particle(occup, site, coeff))
            ts[occup_to_index(occup)] += (coeff * factor) * s[i];
    }
}

template void Hamiltonian::add_particle_state(FARRAY s, int site, FARRAY &ts, FLTYPE factor);
template void Hamiltonian::add_particle_state(FARRAY s, int site, VectorXcd &ts, FLTYPE factor);
template void Hamiltonian::add_particle_state(FARRAY s, int site, VectorXcd &ts, std::complex<FLTYPE> factor);
template void Hamiltonian::add_particle_state(VectorXcd s, int site, VectorXcd &ts, FLTYPE factor);
template void Hamiltonian::add_particle_state(VectorXcd s, int site, VectorXcd &ts, std::complex<FLTYPE> factor);

/*
 * removes a particle to a given state (ts) at site (site) with a prefactor (factor)
 */
template<typename T, typename U, typename G>
void Hamiltonian::remove_particle_state(G s, int site, T &ts, U factor)
{
    if(is_particle_conserving)
    {
        if(current_particle_space < 0)
        {
            ts[0] = 1.0;
            return;
        }
    }
    for(unsigned int i = 0; i < s.size(); ++i)
    {
        current_particle_space = particle_number;
        unsigned int occup = index_to_occup(i);
        FLTYPE coeff = 1.0;
        current_particle_space = particle_number - 1;
        if(remove_particle(occup, site, coeff))
            ts[occup_to_index(occup)] += (coeff * factor) * s[i];
    }
}

template void Hamiltonian::remove_particle_state(FARRAY s, int site, FARRAY &ts, FLTYPE factor);
template void Hamiltonian::remove_particle_state(FARRAY s, int site, VectorXcd &ts, FLTYPE factor);
template void Hamiltonian::remove_particle_state(FARRAY s, int site, VectorXcd &ts, std::complex<FLTYPE> factor);
template void Hamiltonian::remove_particle_state(VectorXcd s, int site, VectorXcd &ts, FLTYPE factor);
template void Hamiltonian::remove_particle_state(VectorXcd s, int site, VectorXcd &ts, std::complex<FLTYPE> factor);
