/*
This code is part of DIETA

Authored by Kirill Alpin
*/

#include <fstream>
#include <iostream>
#include <sstream>
#include "parser.h"

/*
 * create an instance of LatticeGraph
 */
LatticeGraph::LatticeGraph(std::string graph_filename, std::string params_filename)
{
    //default are non-Nambu Greens functions
    is_nambu = false;
    //default is a Hamiltonian without particle number conservation
    is_particle_conserving = false;
    is_half_filling = false;
    //default is to use ordinary CPT Greens function densities
    cluster_dens = false;
    better_densities = false;
    //default is no application of no particle hole transformation
    ph_transform = false;
    //default is a non parity conserving Hamiltonian
    is_parity_conserving = false;
    //default is to not use currents as mean field parameter
    use_current_decoupling = false;
    //default is no sign hardcore boson sgin in the current mean field decoupling
    hardcore_current = false;
    //default is that the original F field is equal to the applied one
    F_field_zero_in_continuum_limit = false;
    real_F_field_is_mean_field = false;
    //default is a one dimensional system
    is_onedim = true;
    //default is an automatic grand potential integration mode (see potential.cxx)
    intMode = GPIntegrationMode::AUTOMATIC;
    //default is a system of hardcore bosons => maximum number of particles per site = 1 and bosons
    max_particle = 1;
    particleType = ParticleType::BOSON;
    //default are zero bath sites
    num_bath = 0;
    //default values for omega sweep parameters for CPT spectra
    omega_start = 0.0;
    omega_end = 3.0;
    omega_steps = 120;
    //first load the Hamiltonian and starting mean field parameters 
    //(see also the function Hamiltonian::Hamiltonian in hamiltonian.cxx)
    std::cout << "Load " << params_filename << std::endl;
    std::ifstream infile(params_filename);
    std::string line;
    //iterate over the lines in params_filename
    while(std::getline(infile, line))
    {
        // a / indicates a comment
        if(line[0] != '/')
        {
            std::string pname;
            std::string para;
            std::istringstream iss(line);
            //read the parameter name (pname) and its value (para)
            iss >> pname >> para;
            int namid = param_names.size();
            //add parameter name to dictionary
            param_names.push_back(pname);
            //parse para to obtain its value
            params.push_back(parse_float(para));
            if(iss.rdbuf()->in_avail() != 0)
            {
                char c;
                iss.get(c);
                if(c != '\n')
                {
                    //if this line contains a second value, 
                    //its value is the complex part of the current parameter
                    
                    //indicate that the Hamiltonian has complex parameters
                    is_complex = true;
                    //parse the imaginary value of the current parameters
                    std::string cpara;
                    iss >> cpara;
                    params[namid] = std::complex<FLTYPE>(params[namid].real(), parse_float(cpara));
                }
            }
            std::cout << pname << "\t" << params[namid] << std::endl;
        }
    }
    infile.close();
    
    
    //load the graph file, which contains the structure of the systems Hamiltonian
    //and extra flags
    std::cout << "Load " << graph_filename << std::endl;
    std::ifstream infile2(graph_filename);
    //default values for the artificial peak broadening and number of k steps for CPT spectra
    eps = 0.05;
    k_steps = 60;
    
    int phase = 0;
    std::string line2;
    while(std::getline(infile2, line2))
    {
        if(line2[0]=='/')
        {
            // a / indicates a comment
            continue;
        }
        else if(line2[0]=='#')
        {
            //# indicates a change to the next phase
            phase++;
        }
        else
        {
            std::istringstream iss(line2);
            switch(phase) 
            {
                case 1:
                    {
                        //phase 1 describes the sites of the Hamiltonian
                        Vec pos;
                        pos.num = site_pos.size();
                        pos.is_edge = false;
                        std::string param_name;
                        //get position of the site
                        for(int i = 0; i < 3; ++i)
                        {
                            std::string expression_string;
                            iss >> expression_string;
                            pos.pos[i] = parse_float(expression_string);
                        }
                        //read type of the site
                        //and the parameter name of the on-site chemical potential
                        iss >> pos.type >> param_name;
                        //retrieve the parameter index of the on-site chemical potential
                        for(int i = 0; i < param_names.size(); ++i)
                                if(param_names[i] == param_name)
                                    pos.ind = i;
                        //add to the list of sites
                        site_pos.push_back(pos);
                    }
                    break;
                case 2:
                    {
                        //phase 2 describes the interactions of the different sites
                        Edge edge;
                        //read the indices of the sites involved in the current interaction
                        iss >> edge.a >> edge.b >> edge.type;
                        //the interaction is decoupled via mean field, if the interaction type is M
                        if(edge.type == 'M')
                        {
                            std::string param_name;
                            //get the translation vector index to a possible different cluster
                            //and type of interaction (see hamiltonian.cxx for a description of types)
                            iss >> edge.trans_index >> edge.type;
                            //read the parameter name and index of for the current interaction strength
                            iss >> param_name;
                            for(int i = 0; i < param_names.size(); ++i)
                                if(param_names[i] == param_name)
                                    edge.ind = i;
                            //add to the list of mean field decoupling edges
                            mean_edges.push_back(edge);
                        }
                        else if(edge.type == 'V') //the interaction is written in the V matrix, if its type is V
                        {
                            //get the translation vector index to a possible different cluster
                            //and type of interaction (see vmatrix.cxx for a description of types)
                            iss >> edge.trans_index >> edge.type;
                            std::string param_name;
                            //read the parameter name and index of for the current interaction strength
                            iss >> param_name;
                            for(int i = 0; i < param_names.size(); ++i)
                                if(param_names[i] == param_name)
                                    edge.ind = i;
                            //add to the list of V matrix edges
                            inter_edges.push_back(edge);
                        }
                        else //if the interaction type is none of the above, it is an ordinary interaction in the Hamiltonian
                        {
                            std::string param_name;
                            //read the parameter name and index of for the current interaction strength
                            iss >> param_name;
                            for(int i = 0; i < param_names.size(); ++i)
                                if(param_names[i] == param_name)
                                    edge.ind = i;
                            //add to the list of Hamiltonian interaction edges
                            edges.push_back(edge);
                        }
                    }
                    break;
                case 3:
                    {
                        //phase 3 reads the possible inter cluster translation vectors
                        Vec trans;
                        for(int i = 0; i < 3; ++i)
                        {
                            std::string expression_string;
                            iss >> expression_string;
                            trans.pos[i] = parse_float(expression_string);
                        }
                        translations.push_back(trans);
                    }
                    break;
                case 4:
                    {
                        //phase 4 reads the k points the CPT spectrum visits
                        Vector3d sympoint;
                        for(int i = 0; i < 3; ++i)
                        {
                            std::string expression_string;
                            iss >> expression_string;
                            sympoint[i] = parse_float(expression_string);
                        }
                        sym_points.push_back(sympoint);
                    }
                    break;
                case 5:
                    {
                        //phase 5 reads the artificial peak broadening value for CPT spectra
                        std::string expression_string;
                        iss >> expression_string;
                        eps = parse_float(expression_string);
                    }
                    break;
                case 6:
                    {
                        //phase 6 reads the number of k steps between the ones in sym_points
                        std::string expression_string;
                        iss >> expression_string;
                        k_steps = parse_float(expression_string);
                    }
                    break;
                case 7:
                    {
                        //phase 7 reads the parametrization of the Brillouin zone for integration
                        //if none are given, a 1D integration from -pi to pi for a 1D system
                        //and a 2D integration from -pi to pi in both coordinates for a 2D system is assumed
                        
                        //else parametrize the Brillouin zone using distorted squares using four k vectors
                        if(poly_bz.size() % 4 == 0 && poly_bz.size() != 0 && poly_bz.size() / 4 > poly_weight.size())
                        {
                            //read the weigth/prefactor of the current square in the integration
                            std::string expression_string;
                            iss >> expression_string;
                            poly_weight.push_back(parse_float(expression_string));
                        }
                        else
                        {
                            //read a k vector of the edge of the current square
                            Vector3d s;
                            for(int i = 0; i < 3; ++i)
                            {
                                std::string expression_string;
                                iss >> expression_string;
                                s[i] = parse_float(expression_string);
                            }
                            poly_bz.push_back(s);
                        }
                    }
                    break;
                case 8:
                    {
                        //phase 8 reads extra flags
                        std::string parse_string;
                        iss >> parse_string;
                        if(parse_string == "half_filling")
                        {
                            //use Hilbert space truncation
                            //Hamiltonian must be particle number conserving
                            //and groundstate is at half filling
                            std::cout << "Half-filling" << std::endl;
                            if(site_pos.size() % 2 != 0)
                            {
                                std::cout << FYEL("WARNING!") << " Exact half-filling not possible, because number of sites is odd." << std::endl;
                            }
                            if(site_pos.size() > 9)
                                is_particle_conserving = true;
                            is_half_filling = true;
                            particle_number = site_pos.size() / 2;
                        }
                        
                        if(parse_string == "cluster_densities")
                        {
                            //use cluster densities for mean field decoupling
                            cluster_dens = true;
                        }
                        if(parse_string == "vacuum")
                        {
                            //use Hilbert space truncation
                            //Hamiltonian must be particle number conserving
                            //and groundstate is the particle vacuum
                            is_particle_conserving = true;
                            particle_number = 0;
                        }
                        if(parse_string == "nambu")
                        {
                            //use the Nambu Greens function formalism
                            is_nambu = true;
                        }
                        if(parse_string == "better_densities")
                        {
                            //use better densities for mean field decoupling
                            better_densities = true;
                        }
                        if(parse_string == "ph_transform")
                        {
                            //apply a particle hole transformation on the states
                            ph_transform = true;
                        }
                        if(parse_string == "parity_conserving")
                        {
                            //use Hilbert space truncation
                            //Hamiltonian must be parity conserving
                            is_parity_conserving = true;
                        }
                        if(parse_string == "current_decoupling")
                        {
                            //use current mean field decoupling
                            //experimental feature
                            use_current_decoupling = true;
                            hardcore_current = true;
                        }
                        if(parse_string == "nonhardcore_current_decoupling")
                        {
                            //do not apply hardcore boson sign in current mean field decoupling
                            //experimental feature
                            if(max_particle == 1)
                                hardcore_current = false;
                            use_current_decoupling = true;
                        }
                        if(parse_string == "F_field_zero_in_continuum_limit")
                        {
                            //used in the boson condensation term of the grand potential
                            //describes a spontanious broken U(1) symmetry
                            //since F field is zero in the original system
                            F_field_zero_in_continuum_limit = true;
                        }
                        if(parse_string == "real_F_field_is_mean_field")
                        {
                            //used in the boson condensation term of the grand potential
                            //describes a system where the F field stems only from mean field decoupling
                            real_F_field_is_mean_field = true;
                        }
                        if(parse_string == "bath")
                        {
                            //number of bath sites
                            iss >> num_bath;
                        }
                        if(parse_string == "omega_start")
                        {
                            //starting value of omega in CPT spectra
                            iss >> omega_start;
                        }
                        if(parse_string == "omega_end")
                        {
                            //end value of omega in CPT spectra
                            iss >> omega_end;
                        }
                        if(parse_string == "omega_steps")
                        {
                            //omega resolution
                            iss >> omega_steps;
                        }
                        if(parse_string == "gp_automatic")
                        {
                            //automatically choose a grand potential integration mode
                            intMode = GPIntegrationMode::AUTOMATIC;
                        }
                        if(parse_string == "gp_exact")
                        {
                            //use the exact integration of eigenvalues over the BZ 
                            //to determine the grand potential
                            intMode = GPIntegrationMode::EXACT;
                        }
                        if(parse_string == "gp_contour")
                        {
                            //use a contour integration to determine the grand potential
                            intMode = GPIntegrationMode::CONTOUR;
                        }
                        if(parse_string == "max_particle")
                        {
                            //choose maximal amount of particles per site_pos
                            //only relevant for bosons
                            //when max_particle==1 and boson flag is set => system consists of hardcore bosons
                            iss >> max_particle;
                            if(max_particle > 1)
                                hardcore_current = false;
                        }
                        if(parse_string == "boson")
                        {
                            //system consists of bosons
                            particleType = ParticleType::BOSON;
                        }
                        if(parse_string == "fermion")
                        {
                            //system consists of fermions
                            //mix of bosons and fermions not supported yet
                            //to add support:
                            //1. add type (boson/fermion) specification to site_pos vector, to identify type of site
                            //2. change fermi sign calculation in Hamiltonian::add_particle and Hamiltonian::remove_particle depending on the site type
                            //3. change fermion sign in Greens function Greens::bfh_sign depending on site type
                            particleType = ParticleType::FERMION;
                        }
                    }
                    break;
            }
        }
    }
    
    //legacy code??
    //F indicates edge of the cluster?
    for(auto const& edge : mean_edges)
    {
        if(edge.type == 'F')
        {
            site_pos[edge.a].is_edge = true;
            if(edge.b != -1)
                site_pos[edge.b].is_edge = true;
        }
    }
    
    //check if system is one dimensional or not using the cluster translation vectors
    Vector3d main_dim = Vector3d(0,0,0);
    for(auto const& trans : translations)
    {
        if(main_dim.dot(main_dim) > 1e-10 && trans.pos.dot(trans.pos) > 1e-10 && trans.pos.dot(main_dim) < 1e-10)
        {
            std::cout << "System is not one dimensional." << std::endl;
            is_onedim = false;
        }
        if(main_dim.dot(main_dim) < 1e-10 && trans.pos.dot(trans.pos) > 1e-10)
        {
            main_dim = trans.pos;
        }
    }
    
    if(is_onedim)
        std::cout << "System is one dimensional." << std::endl;
    
    infile2.close();
}

/*
 * used to parse floats form files
 * the string can contain mathematical expressions like cos,sin,sqrt...
 */
FLTYPE LatticeGraph::parse_float(std::string expression_string)
{
    typedef exprtk::symbol_table<FLTYPE> symbol_table_t;
    typedef exprtk::expression<FLTYPE>     expression_t;
    typedef exprtk::parser<FLTYPE>             parser_t;
    
    symbol_table_t symbol_table;
    symbol_table.add_constants();
    
    expression_t expression;
    expression.register_symbol_table(symbol_table);

    parser_t parser;
    parser.compile(expression_string,expression);
    
    return expression.value();
}
