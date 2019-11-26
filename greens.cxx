/*
This code is part of DIETA

Authored by Kirill Alpin
*/

#include <bitset>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>
#include "greens.h"


/*
 * initialize Greens. if (precalc) is true, the Greens function is calculated directly
 */
template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
Greens<VALTYPE, VECTYPE, MATTYPE>::Greens(Hamiltonian* hamiltonian, bool precalc) : lanczos(hamiltonian)
{
    highenergy = 0.0;
    H = hamiltonian;
    eps = hamiltonian->graph.eps;
    std::cout << "Hilbert dimension: " << H->dim << std::endl;
    //apply full diagonalization when Hilbert space is smaller than 256
    full_diag = H->dim <= 256;
    if(full_diag)
        std::cout << "Full diagonalization mode" << std::endl;
    if(H->is_complex)
        std::cout << "Complex mode"<< std::endl;
    if(precalc)
        prepare_greens();
}

template Greens<FLTYPE, FARRAY, MatrixXd>::Greens(Hamiltonian* hamiltonian, bool precalc);
template Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::Greens(Hamiltonian* hamiltonian, bool precalc);

/*
 * initialize Greens from (file)
 */
template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
Greens<VALTYPE, VECTYPE, MATTYPE>::Greens(Hamiltonian* hamiltonian, std::string file) : lanczos(hamiltonian)
{
    highenergy = 0.0;
    H = hamiltonian;
    eps = hamiltonian->graph.eps;
    std::cout << "Hilbert dimension: " << H->dim << std::endl;
    full_diag = H->dim <= 256;
    if(H->is_complex)
        std::cout << "Complex mode"<< std::endl;
    load_params(file);
}

template Greens<FLTYPE, FARRAY, MatrixXd>::Greens(Hamiltonian* hamiltonian, std::string file);
template Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::Greens(Hamiltonian* hamiltonian, std::string file);

/*
 * prepare greens function for use. the groundstate, Q matrices and energies are calculated
 */
template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
void Greens<VALTYPE, VECTYPE, MATTYPE>::prepare_greens(bool no_interaction)
{
    H->no_interaction = no_interaction;
    if(full_diag)
    {
        //prepare Greens function using full diagonalisation
        if(H->is_particle_conserving || H->is_parity_conserving)
        {
            std::stringstream ss;
            ss << "Particle or parity conservation and full diagonalisation is not supported.";
            throw std::invalid_argument(ss.str());
        }
        //build Hamiltonian matrix by apply all orthogonal state vectors
        MATTYPE H_matrix = MATTYPE::Zero(H->dim, H->dim);
        H->set_particle_space(H->particle_number);
        for(int i = 0; i < H->dim; ++i)
        {
            VECTYPE state = VECTYPE::Zero(H->dim);
            state[i] = 1.0;
            H_matrix.col(i) = (*H) * state;
        }
        
        //check if Hamiltonian is hermitian
        if((H_matrix-H_matrix.adjoint()).cwiseAbs().maxCoeff() > 1.0e-10)
        {
            std::stringstream ss;
            ss << FRED("ERROR!") << "Matrix is not hermitian! Dein Gehirn ist zu klein!";
            throw std::invalid_argument(ss.str());
        }
        
        //get eigenvalues and eigenvectors
        SelfAdjointEigenSolver<MATTYPE> es(H_matrix);
        lanczos.energies = es.eigenvalues();
        MATTYPE evecs = es.eigenvectors();
        int sm = 0;
        //get groundstate energy
        FLTYPE genergy = std::real(lanczos.energies[sm]);
        highenergy = std::abs(lanczos.energies[lanczos.energies.size() - 1] - genergy);
        for(int i = 0; i < lanczos.energies.size(); ++i)
            if(std::real(lanczos.energies[i]) < genergy)
            {
                genergy = std::real(lanczos.energies[i]);
                sm = i;
            }
        
        if(lanczos.energies.size() == 2)
        {
            std::cout << "Cluster energies: " << lanczos.energies[0] << "\t" << lanczos.energies[1] << std::endl;
        }
        
        if(lanczos.energies.size() > 2)
        {
            std::cout << "Cluster energies: " << lanczos.energies[0] << "\t" << lanczos.energies[1] << "\t" << lanczos.energies[2] << std::endl;
        }
        
        //check degeneracy
        deg = 1;
        for(int i = 1; i < H->dim; ++i)
        {
            if(std::abs(lanczos.energies[i] - genergy) < 1.0e-10)
                deg++;
            else
                break;
        }
        if(deg != 1)
        {
            std::cout << FYEL("WARNING!") << " Groundstate is " << deg << " degenerate." << std::endl;
        }
        //apply epsilon broadening to the groundstate energy
        groundenergy = std::complex<FLTYPE>(genergy, eps);
        //get groundstate
        groundstate = VECTYPE::Zero(H->dim);
        for(int i = 0; i < deg; ++i)
            groundstate += evecs.col(i);
        groundstate /= deg;
        //(highest) is the highest relevant state index
        int highest = 0;
        
        int num_states = H->dim;
        
        if(deg == 1)
        {
            //build Q matrices when groundstate is non degenerate
            Q_matrix_elec = MATTYPE::Zero(H->numsites, num_states);
            Q_matrix_hole = MATTYPE::Zero(H->numsites, num_states);
            VECTYPE state = VECTYPE::Zero(H->dim);
            for(int i = 0; i < num_states; ++i)
            {
                double m = 0.0;
                for(int j = 0; j < H->numsites; ++j)
                {
                    state = VECTYPE::Zero(H->dim);
                    H->add_particle_state(groundstate, j, state);
                    Q_matrix_elec(j,i) = state.dot(evecs.col(i));
                    if(m < std::abs(Q_matrix_elec(j,i)))
                        m = std::abs(Q_matrix_elec(j,i));
                    state = VECTYPE::Zero(H->dim);
                    H->remove_particle_state(groundstate, j, state);
                    Q_matrix_hole(j,i) = state.dot(evecs.col(i));
                    if(m < std::abs(Q_matrix_hole(j,i)))
                        m = std::abs(Q_matrix_hole(j,i));
                }
                if(m > 0.1)
                    highest = i;
            }
        }
        else
        {
            //build Q matrices when groundstate is degenerate
            Q_matrix_elec = MATTYPE::Zero(H->numsites * deg, num_states);
            Q_matrix_hole = MATTYPE::Zero(H->numsites * deg, num_states);
            VECTYPE state = VECTYPE::Zero(H->dim);
            
            for(int d = 0; d < deg; ++d)
            {
                for(int i = 0; i < num_states; ++i)
                {
                    double m = 0.0;
                    for(int j = 0; j < H->numsites; ++j)
                    {
                        state = VECTYPE::Zero(H->dim);
                        H->add_particle_state((VECTYPE)evecs.col(d), j, state);
                        Q_matrix_elec(j+d*H->numsites,i) = state.dot(evecs.col(i));
                        if(m < std::abs(Q_matrix_elec(j+d*H->numsites,i)))
                            m = std::abs(Q_matrix_elec(j+d*H->numsites,i));
                    }
                    if(m > 0.1)
                        highest = i;
                }
                for(int i = 0; i < num_states; ++i)
                {
                    double m = 0.0;
                    for(int j = 0; j < H->numsites; ++j)
                    {       
                        state = VECTYPE::Zero(H->dim);
                        H->remove_particle_state((VECTYPE)evecs.col(d), j, state);
                        Q_matrix_hole(j+d*H->numsites,i) = state.dot(evecs.col(i));
                        if(m < std::abs(Q_matrix_hole(j+d*H->numsites,i)))
                            m = std::abs(Q_matrix_hole(j+d*H->numsites,i));
                    }
                    if(m > 0.1 && i > highest)
                        highest = i;
                }
            }
        }
        
        //get highest energy from highest relevant state. used for integration over omega
        highenergy = std::max(std::abs(lanczos.energies[highest] - genergy), highenergy);
        std::cout << "Goundstate energy: " << genergy << std::endl;
        std::cout << "Max energy: " << highenergy << std::endl;
        
    }
    else
    {
        //prepare Greens function using band lanczos
        //12345 is the seed of the lanczos runs
        //make first run to get the tridiagonal Matrix representation
        MATTYPE M = lanczos.get_tridiagonal(12345);
        FLTYPE genergy = 0.0;
        //make groundstate run
        groundstate = lanczos.get_groundstate(12345, M, genergy);
        groundenergy = std::complex<FLTYPE>(genergy, eps);
        
        //check groundstate sanity
        VECTYPE tmp = (*H) * groundstate;
        FLTYPE ch_energy = std::real(tmp.dot(groundstate));
        if(std::abs(ch_energy - genergy) < 1.0e-10)
            std::cout << "INFO: " << FGRN("GROUNDSTATE CHECK SUCCESSFUL") << "\tError: " <<  std::abs(ch_energy - genergy) << std::endl;
        else
            std::cout << FYEL("WARNING") << ": " << FRED("GROUNDSTATE CHECK FAILED") << "\tError: " <<  std::abs(ch_energy - genergy) << std::endl;
        
        //(phi_elec) and (phi_hole) hold the vectors a^dag_i|groundstate> and a_i|groundstate> 
        // used as a starting basis for the next band lanczos iterations
        MATTYPE phi_elec = MATTYPE::Zero(H->numsites, H->dim_part);
        MATTYPE phi_hole = MATTYPE::Zero(H->numsites, H->dim_hole);
        for(int i = 0; i < H->numsites; ++i)
        {
            VECTYPE tmp = VECTYPE::Zero(H->dim_hole);
            H->remove_particle_state(groundstate, i, tmp);
            phi_hole.row(i) = tmp.adjoint();
            tmp = VECTYPE::Zero(H->dim_part);
            H->add_particle_state(groundstate, i, tmp);
            phi_elec.row(i) = tmp.adjoint();
        }
        
        //band lanczos run to determine (T_elec), (Q_matrix_elec) and (Q_matrix_elec_h)
        lanczos.get_Qmatrix(groundstate, T_elec, Q_matrix_elec, Q_matrix_elec_h, false, phi_elec, phi_hole);
        //band lanczos run to determine (T_hole), (Q_matrix_hole_e) and (Q_matrix_hole)
        lanczos.get_Qmatrix(groundstate, T_hole, Q_matrix_hole_e, Q_matrix_hole, true, phi_elec, phi_hole);
        
        //use T_elec eigenbasis for the particle Q matrices
        SelfAdjointEigenSolver<MATTYPE> es(T_elec);
        T_elec_energies = es.eigenvalues();
        Q_matrix_elec = Q_matrix_elec * es.eigenvectors();
        Q_matrix_elec_h = Q_matrix_elec_h * es.eigenvectors();
        
        //use T_hole eigenbasis for the hole Q matrices
        SelfAdjointEigenSolver<MATTYPE> es2(T_hole);
        T_hole_energies = es2.eigenvalues();
        Q_matrix_hole_e = Q_matrix_hole_e * es2.eigenvectors();
        Q_matrix_hole = Q_matrix_hole * es2.eigenvectors();
        
        //get highes relevant energy
        int highest_elec_energy = 0;
        int highest_elec = 0;
        for(int i = 0; i < Q_matrix_elec.cols(); ++i)
        {
            double m = 0.0;
            for(int j = 0; j < H->numsites; ++j)
            {
                if(m < std::abs(Q_matrix_elec(j,i)))
                    m = std::abs(Q_matrix_elec(j,i));
            }
            if(m > 0.1)
                highest_elec = i;
            if(m > 0.1)
                highest_elec_energy = i;
        }
        int highest_hole_energy = 0;
        int highest_hole = 0;
        for(int i = 0; i < Q_matrix_hole.cols(); ++i)
        {
            double m = 0.0;
            for(int j = 0; j < H->numsites; ++j)
            {
                if(m < std::abs(Q_matrix_hole(j,i)))
                    m = std::abs(Q_matrix_hole(j,i));
            }
            if(m > 0.1)
                highest_hole = i;
            if(m > 0.1)
                highest_hole_energy = i;
        }
        
        std::cout << "Goundstate energy: " << genergy << std::endl;
        highenergy = std::max(std::abs(std::max(T_hole_energies[highest_hole_energy],
            T_elec_energies[highest_elec_energy]) - genergy), highenergy); 
        
        
        std::cout << "Max energy: " << highenergy << std::endl;
    }
}


template void Greens<FLTYPE, FARRAY, MatrixXd>::prepare_greens(bool no_interaction);
template void Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::prepare_greens(bool no_interaction);

/*
 * alternative way of calculating cluster densities
 * unused for now
 */
template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
FARRAY Greens<VALTYPE, VECTYPE, MATTYPE>::clusterDensities()
{
    VectorXd clusterDens = ((FARRAY::Ones(Q_matrix_elec.rows()) - ((Q_matrix_elec * Q_matrix_elec.adjoint()).diagonal().real() - bfh_sign(Q_matrix_hole * Q_matrix_hole.adjoint()).diagonal().real())) / 2.0);
    
    if(clusterDens.size() > H->numsites)
    {
        //deal with degeneracy
        VectorXd tmp = clusterDens;
        clusterDens = VectorXd::Zero(H->numsites);
        for(int j = 0; j < (tmp.size() / H->numsites); ++j)
            clusterDens += tmp.segment(j * H->numsites, H->numsites);
        clusterDens /= tmp.size() / H->numsites;
    }
    return clusterDens;
}

template FARRAY Greens<FLTYPE, FARRAY, MatrixXd>::clusterDensities();
template FARRAY Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::clusterDensities();

/*
 * calculate CPT Greensfunction for omega to infinity
 */
template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
MatrixXcd Greens<VALTYPE, VECTYPE, MATTYPE>::omegaInfLimit()
{
    MatrixXcd G;
    VectorXcd diag_elec = VectorXcd::Ones(Q_matrix_elec.cols()); //omega / (omega + O(1)) -> 1 for omega -> inf
    VectorXcd diag_hole = VectorXcd::Ones(Q_matrix_hole.cols());
    if(deg == 1)
    {
        G = (Q_matrix_elec * diag_elec.asDiagonal() * Q_matrix_elec.adjoint()) - bfh_sign((Q_matrix_hole * diag_hole.asDiagonal() * Q_matrix_hole.adjoint()));
    }
    else
    {
        G = MatrixXcd::Zero(H->numsites, H->numsites);
        int numstates = Q_matrix_elec.cols();
        for(int d = 0; d < deg; ++d)
        {
            G += (Q_matrix_elec.block(H->numsites*d, 0, H->numsites, numstates) * diag_elec.asDiagonal() * Q_matrix_elec.block(H->numsites*d, 0, H->numsites, numstates).adjoint()) - bfh_sign((Q_matrix_hole.block(H->numsites*d, 0, H->numsites, numstates) * diag_hole.asDiagonal() * Q_matrix_hole.block(H->numsites*d, 0, H->numsites, numstates).adjoint()));
        }
        G = G / deg;
    }
    if(!H->graph.is_nambu)
        return G;
    else
    {
        MatrixXcd F;
        MatrixXcd F_conj;
        if(deg == 1)
        {
            if(full_diag)
            {
                F_conj = (Q_matrix_hole * diag_elec.asDiagonal() * Q_matrix_elec.adjoint()) - bfh_sign((Q_matrix_elec * diag_hole.asDiagonal() * Q_matrix_hole.adjoint()));
                F = (Q_matrix_elec * diag_elec.asDiagonal() * Q_matrix_hole.adjoint()) - bfh_sign((Q_matrix_hole * diag_hole.asDiagonal() * Q_matrix_elec.adjoint()));
            }
            else
            {
                F_conj = (Q_matrix_elec_h * diag_elec.asDiagonal() * Q_matrix_elec.adjoint()) - bfh_sign((Q_matrix_hole_e * diag_hole.asDiagonal() * Q_matrix_hole.adjoint()));
                F = (Q_matrix_elec * diag_elec.asDiagonal() * Q_matrix_elec_h.adjoint()) - bfh_sign((Q_matrix_hole * diag_hole.asDiagonal() * Q_matrix_hole_e.adjoint()));
            }
        }
        else
        {
            F = MatrixXcd::Zero(H->numsites, H->numsites);
            F_conj = MatrixXcd::Zero(H->numsites, H->numsites);
            int numstates = Q_matrix_elec.cols();
            for(int d = 0; d < deg; ++d)
            {
                if(full_diag)
                {
                    F_conj += (Q_matrix_hole.block(H->numsites*d, 0, H->numsites, numstates) * diag_elec.asDiagonal() * Q_matrix_elec.block(H->numsites*d, 0, H->numsites, numstates).adjoint()) - bfh_sign((Q_matrix_elec.block(H->numsites*d, 0, H->numsites, numstates) * diag_hole.asDiagonal() * Q_matrix_hole.block(H->numsites*d, 0, H->numsites, numstates).adjoint()));
                
                    F += (Q_matrix_elec.block(H->numsites*d, 0, H->numsites, numstates) * diag_elec.asDiagonal() * Q_matrix_hole.block(H->numsites*d, 0, H->numsites, numstates).adjoint()) - bfh_sign((Q_matrix_hole.block(H->numsites*d, 0, H->numsites, numstates) * diag_hole.asDiagonal() * Q_matrix_elec.block(H->numsites*d, 0, H->numsites, numstates).adjoint()));
                }
                else
                {
                    std::stringstream ss;
                    ss << FRED("ERROR!") << "You know what you did!";
                    throw std::invalid_argument(ss.str());
                }
            }
            F = F / deg;
            F_conj = F_conj / deg;
        }
        
        MatrixXcd ret(H->numsites * 2, H->numsites * 2);
        ret << G, F,
            F_conj, -bfh_sign(G); //minus da omega/(-omega + O(1)) -> -1 fuer omega -> inf
        return ret;
    }
}

template MatrixXcd Greens<FLTYPE, FARRAY, MatrixXd>::omegaInfLimit();
template MatrixXcd Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::omegaInfLimit();


/* =========================================================================
 * code for saving and loading Greens function data
 */

template void Greens<FLTYPE, FARRAY, MatrixXd>::parse_vector(std::ofstream& ofile, FARRAY vec);
template void Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::parse_vector(std::ofstream& ofile, FARRAY vec);

template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
void Greens<VALTYPE, VECTYPE, MATTYPE>::parse_vector(std::ofstream& ofile, FARRAY vec)
{
    typename FARRAY::Index size=vec.size();
    ofile.write((char*) (&size), sizeof(typename FARRAY::Index));
    ofile.write((char*) vec.data(), size*sizeof(FLTYPE) );
}

template void Greens<FLTYPE, FARRAY, MatrixXd>::parse_vector(std::ofstream& ofile, VectorXcd vec);
template void Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::parse_vector(std::ofstream& ofile, VectorXcd vec);

template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
void Greens<VALTYPE, VECTYPE, MATTYPE>::parse_vector(std::ofstream& ofile, VectorXcd vec)
{
    typename VectorXcd::Index size=vec.size();
    ofile.write((char*) (&size), sizeof(typename VectorXcd::Index));
    ofile.write((char*) vec.data(), size*sizeof(FLTYPE)*2 );
}

template void Greens<FLTYPE, FARRAY, MatrixXd>::parse_matrix(std::ofstream& ofile, MatrixXd matrix);
template void Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::parse_matrix(std::ofstream& ofile, MatrixXcd matrix);

template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
void Greens<VALTYPE, VECTYPE, MATTYPE>::parse_matrix(std::ofstream& ofile, MATTYPE matrix)
{
    typename MATTYPE::Index rows=matrix.rows(), cols=matrix.cols();
    ofile.write((char*) (&rows), sizeof(typename MATTYPE::Index));
    ofile.write((char*) (&cols), sizeof(typename MATTYPE::Index));
    ofile.write((char*) matrix.data(), rows*cols*sizeof(VALTYPE) );
}

template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
void Greens<VALTYPE, VECTYPE, MATTYPE>::save_params(std::string file)
{
    std::cout << "Save " << file << std::endl;
    std::ofstream ofile(file, std::ios::out | std::ios::binary | std::ios::trunc);
    
    FLTYPE packed_groundenergy[] = {std::real(groundenergy), std::imag(groundenergy)};
    ofile.write(reinterpret_cast<char*>(packed_groundenergy), 2*sizeof(FLTYPE));
    
    parse_vector(ofile, groundstate);
    parse_vector(ofile, T_elec_energies);
    parse_vector(ofile, T_hole_energies);
    parse_matrix(ofile, Q_matrix_elec);
    parse_matrix(ofile, Q_matrix_elec_h);
    parse_matrix(ofile, Q_matrix_hole);
    parse_matrix(ofile, Q_matrix_hole_e);
    ofile.close();
}

template void Greens<FLTYPE, FARRAY, MatrixXd>::save_params(std::string file);
template void Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::save_params(std::string file);

template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
std::string Greens<VALTYPE, VECTYPE, MATTYPE>::next_line(std::ifstream& infile)
{
    std::string line;
    std::getline(infile, line);
    return line;
}

template std::string Greens<FLTYPE, FARRAY, MatrixXd>::next_line(std::ifstream& infile);
template std::string Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::next_line(std::ifstream& infile);

template void Greens<FLTYPE, FARRAY, MatrixXd>::load_vector(std::ifstream& infile, FARRAY& vec);
template void Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::load_vector(std::ifstream& infile, FARRAY& vec);

template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
void Greens<VALTYPE, VECTYPE, MATTYPE>::load_vector(std::ifstream& infile, FARRAY& vec)
{
    typename FARRAY::Index size=0;
    infile.read((char*) (&size),sizeof(typename FARRAY::Index));
    vec.resize(size);
    infile.read( (char *) vec.data() , size*sizeof(FLTYPE) );
}

template void Greens<FLTYPE, FARRAY, MatrixXd>::load_vector(std::ifstream& infile, VectorXcd& vec);
template void Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::load_vector(std::ifstream& infile, VectorXcd& vec);

template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
void Greens<VALTYPE, VECTYPE, MATTYPE>::load_vector(std::ifstream& infile, VectorXcd& vec)
{
    typename VectorXcd::Index size=0;
    infile.read((char*) (&size),sizeof(typename VectorXcd::Index));
    vec.resize(size);
    infile.read( (char *) vec.data() , size*2*sizeof(FLTYPE) );
}

template void Greens<FLTYPE, FARRAY, MatrixXd>::load_matrix(std::ifstream& infile, MatrixXd& matrix);
template void Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::load_matrix(std::ifstream& infile, MatrixXcd& matrix);

template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
void Greens<VALTYPE, VECTYPE, MATTYPE>::load_matrix(std::ifstream& infile, MATTYPE& matrix)
{
    typename MATTYPE::Index rows=0, cols=0;
    infile.read((char*) (&rows),sizeof(typename MATTYPE::Index));
    infile.read((char*) (&cols),sizeof(typename MATTYPE::Index));
    matrix.resize(rows, cols);
    infile.read( (char *) matrix.data() , rows*cols*sizeof(VALTYPE) );
}

template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
void Greens<VALTYPE, VECTYPE, MATTYPE>::load_params(std::string file)
{
    std::cout << "Load " << file << std::endl;
    std::ifstream infile(file, std::ios::binary | std::ios::in);
    
    std::vector<std::complex<FLTYPE>> packed_groundenergy(1);
    infile.read(reinterpret_cast<char*>(packed_groundenergy.data()), 2*sizeof(FLTYPE));
    groundenergy = packed_groundenergy[0];
    
    load_vector(infile, groundstate);
    load_vector(infile, T_elec_energies);
    load_vector(infile, T_hole_energies);
    load_matrix(infile, Q_matrix_elec);
    load_matrix(infile, Q_matrix_elec_h);
    load_matrix(infile, Q_matrix_hole);
    load_matrix(infile, Q_matrix_hole_e);
    infile.close();
}

/*
 * =========================================================================
 */



/*
 * unsused
 * legacy code for calculating Greens functions using continued fractions
 */
template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
template <typename T> 
std::complex<FLTYPE> Greens<VALTYPE, VECTYPE, MATTYPE>::calc_fraction(std::complex<FLTYPE> omega, LanczosParams<T> param, bool hole)
{
    std::complex<FLTYPE> ret = 0.0;
    std::complex<FLTYPE> z = groundenergy;
    if(hole)
        z *= -1.0;
    for(int i = param.a.size()-1; i > -1; --i)
        ret = param.bsqr[i] / (omega + z - param.a[i] - ret);
    return ret;
}

template std::complex<FLTYPE> Greens<FLTYPE, FARRAY, MatrixXd>::calc_fraction(std::complex<FLTYPE> omega, LanczosParams<FLTYPE> param, bool hole);
template std::complex<FLTYPE> Greens<FLTYPE, FARRAY, MatrixXd>::calc_fraction(std::complex<FLTYPE> omega, LanczosParams<std::complex<FLTYPE>> param, bool hole);
template std::complex<FLTYPE> Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::calc_fraction(std::complex<FLTYPE> omega, LanczosParams<std::complex<FLTYPE>> param, bool hole);


/*
 * handles the fermion/boson sign in the Greens function
 */
template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
MatrixXcd Greens<VALTYPE, VECTYPE, MATTYPE>::bfh_sign(MatrixXcd g)
{
    if(H->graph.particleType == ParticleType::BOSON)
        return g.transpose();
    else if(H->graph.particleType == ParticleType::FERMION)
        return -g.transpose();
    else if(H->graph.particleType == ParticleType::HARDCORE)
    {
        MatrixXcd ret = g;
        for(int i = 0; i < g.cols(); ++i)
            ret(i,i) = -g(i,i);
        return ret.transpose();
    }
    return MatrixXcd();
}

template MatrixXcd Greens<FLTYPE, FARRAY, MatrixXd>::bfh_sign(MatrixXcd g);
template MatrixXcd Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::bfh_sign(MatrixXcd g);

/*
 * calculates the cluster Greens function <<a^dag_i, a_j>>
 */
template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
MatrixXcd Greens<VALTYPE, VECTYPE, MATTYPE>::get_greens(std::complex<FLTYPE> omega)
{
    if(full_diag)
    {
        //if full diagonalization is used
        VectorXcd diag_elec = VectorXcd::Zero(Q_matrix_elec.cols());
        VectorXcd diag_hole = VectorXcd::Zero(Q_matrix_elec.cols());
        for(int i = 0; i < Q_matrix_elec.cols(); ++i)
        {
            if(std::abs(omega - lanczos.energies[i] + groundenergy) > 1e-11)
            {
                diag_elec[i] = 1.0 / (omega - lanczos.energies[i] + groundenergy);
                diag_hole[i] = 1.0 / (omega + lanczos.energies[i] - groundenergy);
            }
        }
        if(deg == 1)
        {
            //no groundstate degeneracy code
            if(H->graph.num_bath > 0) //additional bath sites. see senechal paper "An introduction to quantum cluster methods"
                return (((Q_matrix_elec * diag_elec.asDiagonal() * Q_matrix_elec.adjoint()) - bfh_sign((Q_matrix_hole * diag_hole.asDiagonal() * Q_matrix_hole.adjoint()))).inverse() - get_bath_gamma(omega)).inverse();
            else
                return (Q_matrix_elec * diag_elec.asDiagonal() * Q_matrix_elec.adjoint()) - bfh_sign((Q_matrix_hole * diag_hole.asDiagonal() * Q_matrix_hole.adjoint()));
        }
        else
        {
            //groundstate is degenerate. T->0 limit of the Lehmann representation is used
            MatrixXcd G = MatrixXcd::Zero(H->numsites, H->numsites);
            int numstates = Q_matrix_elec.cols();
            for(int d = 0; d < deg; ++d)
            {
                G += (Q_matrix_elec.block(H->numsites*d, 0, H->numsites, numstates) * diag_elec.asDiagonal() * Q_matrix_elec.block(H->numsites*d, 0, H->numsites, numstates).adjoint()) - bfh_sign((Q_matrix_hole.block(H->numsites*d, 0, H->numsites, numstates) * diag_hole.asDiagonal() * Q_matrix_hole.block(H->numsites*d, 0, H->numsites, numstates).adjoint()));
            }
            if(H->graph.num_bath > 0) //additional bath sites. see senechal paper "An introduction to quantum cluster methods"
                return ((G / deg).inverse() - get_bath_gamma(omega)).inverse();
            else
                return G / deg;
        }
    }
    else
    {
        //determine cluster Greens function when lanczos was used
        VectorXcd diag_elec = VectorXcd::Zero(T_elec_energies.size());
        VectorXcd diag_hole = VectorXcd::Zero(T_hole_energies.size());
        for(int i = 0; i < T_elec_energies.size(); ++i)
            diag_elec[i] = 1.0 / (omega - T_elec_energies[i] + groundenergy);
        for(int i = 0; i < T_hole_energies.size(); ++i)
            diag_hole[i] = 1.0 / (omega + T_hole_energies[i] - groundenergy);
        MatrixXcd G = (Q_matrix_elec * diag_elec.asDiagonal() * Q_matrix_elec.adjoint()) - bfh_sign((Q_matrix_hole * diag_hole.asDiagonal() * Q_matrix_hole.adjoint()));
        
        if(H->graph.num_bath > 0) //additional bath sites. see senechal paper "An introduction to quantum cluster methods"
            return (G.inverse() - get_bath_gamma(omega)).inverse();
        else
            return G;
    }
}

template MatrixXcd Greens<FLTYPE, FARRAY, MatrixXd>::get_greens(std::complex<FLTYPE> omega);
template MatrixXcd Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::get_greens(std::complex<FLTYPE> omega);

/*
 * calculates the anomalous cluster Greens function <<a_i, a_j>>. when conjugate is true <<a^dag_i, a^dag_j>>
 */
template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
MatrixXcd Greens<VALTYPE, VECTYPE, MATTYPE>::get_anom_greens(std::complex<FLTYPE> omega, bool conjugate)
{
    const std::complex<FLTYPE> img(0, 1);
    if(full_diag)
    {
        //if full diagonalization is used
        VectorXcd diag_elec = VectorXcd::Zero(Q_matrix_elec.cols());
        VectorXcd diag_hole = VectorXcd::Zero(Q_matrix_elec.cols());
        for(int i = 0; i < Q_matrix_elec.cols(); ++i)
        {
            if(std::abs(omega - lanczos.energies[i] + groundenergy) > 1e-11)
            {
                diag_elec[i] = 1.0 / (omega - lanczos.energies[i] + groundenergy);
                diag_hole[i] = 1.0 / (omega + lanczos.energies[i] - groundenergy);
            }
        }
        if(deg == 1)
        {
            //no groundstate degeneracy code
            if(conjugate)
                return (Q_matrix_hole * diag_elec.asDiagonal() * Q_matrix_elec.adjoint()) - bfh_sign((Q_matrix_elec * diag_hole.asDiagonal() * Q_matrix_hole.adjoint()));
            return (Q_matrix_elec * diag_elec.asDiagonal() * Q_matrix_hole.adjoint()) - bfh_sign((Q_matrix_hole * diag_hole.asDiagonal() * Q_matrix_elec.adjoint()));
        }
        else
        {
            //groundstate is degenerate. T->0 limit of the Lehmann representation is used
            MatrixXcd F = MatrixXcd::Zero(H->numsites, H->numsites);
            int numstates = Q_matrix_elec.cols();
            for(int d = 0; d < deg; ++d)
            {
                if(conjugate)
                    F += (Q_matrix_hole.block(H->numsites*d, 0, H->numsites, numstates) * diag_elec.asDiagonal() * Q_matrix_elec.block(H->numsites*d, 0, H->numsites, numstates).adjoint()) - bfh_sign((Q_matrix_elec.block(H->numsites*d, 0, H->numsites, numstates) * diag_hole.asDiagonal() * Q_matrix_hole.block(H->numsites*d, 0, H->numsites, numstates).adjoint()));
                else
                    F += (Q_matrix_elec.block(H->numsites*d, 0, H->numsites, numstates) * diag_elec.asDiagonal() * Q_matrix_hole.block(H->numsites*d, 0, H->numsites, numstates).adjoint()) - bfh_sign((Q_matrix_hole.block(H->numsites*d, 0, H->numsites, numstates) * diag_hole.asDiagonal() * Q_matrix_elec.block(H->numsites*d, 0, H->numsites, numstates).adjoint()));
            }
            return F / deg;
        }
    }
    else 
    {
        //determine anomalous cluster Greens function when lanczos was used
        VectorXcd diag_elec = VectorXcd::Zero(T_elec_energies.size());
        VectorXcd diag_hole = VectorXcd::Zero(T_hole_energies.size());
        for(int i = 0; i < T_elec_energies.size(); ++i)
            diag_elec[i] = 1.0 / (omega - T_elec_energies[i] + groundenergy);
        for(int i = 0; i < T_hole_energies.size(); ++i)
            diag_hole[i] = 1.0 / (omega + T_hole_energies[i] - groundenergy);
        if(conjugate)
            return (Q_matrix_elec_h * diag_elec.asDiagonal() * Q_matrix_elec.adjoint()) - bfh_sign((Q_matrix_hole_e * diag_hole.asDiagonal() * Q_matrix_hole.adjoint()));
        return (Q_matrix_elec * diag_elec.asDiagonal() * Q_matrix_elec_h.adjoint()) - bfh_sign((Q_matrix_hole * diag_hole.asDiagonal() * Q_matrix_hole_e.adjoint()));
    }
}

template MatrixXcd Greens<FLTYPE, FARRAY, MatrixXd>::get_anom_greens(std::complex<FLTYPE> omega, bool conjugate=false);
template MatrixXcd Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::get_anom_greens(std::complex<FLTYPE> omega, bool conjugate=false);

/*
 * calculate full nambu cluster Greens function <<(a^dag_i,a_i)^dag,(a^dag_j,a_j)>>
 */
template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
MatrixXcd Greens<VALTYPE, VECTYPE, MATTYPE>::get_nambu_greens(std::complex<FLTYPE> omega)
{
    MatrixXcd ret(H->numsites * 2, H->numsites * 2);
    ret << get_greens(omega), get_anom_greens(omega),
            get_anom_greens(omega, true), bfh_sign(get_greens(-omega));
    return ret;
}

template MatrixXcd Greens<FLTYPE, FARRAY, MatrixXd>::get_nambu_greens(std::complex<FLTYPE> omega);
template MatrixXcd Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::get_nambu_greens(std::complex<FLTYPE> omega);

/*
 * calculate the noninteracting cluster Greens function. 
 * if (with_mean_field) is true, mean field one body terms are also acounted for in this calculation
 */
template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
MatrixXcd Greens<VALTYPE, VECTYPE, MATTYPE>::get_noninteracting_inverse_greens(std::complex<FLTYPE> omega, bool nambu, bool with_mean_field)
{
    MatrixXcd ret = MatrixXcd::Zero(nambu ? H->numsites * 2 : H->numsites, nambu ? H->numsites * 2 : H->numsites);
    
    for(auto const& site: H->graph.site_pos)
    {
        switch(site.type)
        {
            case 'A': //chemical potential
                {
                    ret(site.num,site.num) += H->graph.params[site.ind];
                    if(nambu)
                        ret(site.num+H->numsites,site.num+H->numsites) += H->graph.params[site.ind];
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
    
    for(auto const& edge: H->graph.edges)
    {
        const std::complex<FLTYPE> i(0, 1);
        switch(edge.type)
        {
            case 'A': //hopping term c^dag*c+c*c^dag
                {
                    ret(edge.b,edge.a) += H->graph.params[edge.ind];
                    if(edge.a != edge.b)
                        ret(edge.a,edge.b) += std::conj(H->graph.params[edge.ind]);
                    if(nambu)
                    {
                        ret(edge.b+H->numsites,edge.a+H->numsites) += H->graph.params[edge.ind];
                        if(edge.a != edge.b)
                            ret(edge.a+H->numsites,edge.b+H->numsites) += std::conj(H->graph.params[edge.ind]);
                    }
                }
                break;
            case 'B': //hopping term c*c+c^dag*c^dag
                {
                    if(nambu)
                    {
                        ret(edge.b,edge.a+H->numsites) += H->graph.params[edge.ind];
                        if(edge.a != edge.b)
                            ret(edge.a,edge.b+H->numsites) += std::conj(H->graph.params[edge.ind]);
                        ret(edge.b+H->numsites,edge.a) += H->graph.params[edge.ind];
                        if(edge.a != edge.b)
                            ret(edge.a+H->numsites,edge.b) += std::conj(H->graph.params[edge.ind]);
                    }
                }
                break;
            default:
                break;
        }
    }
    
    if(with_mean_field)
        for(auto const& edge: H->graph.mean_edges)
        {
            const std::complex<FLTYPE> i(0, 1);
            Vector3d translation = H->graph.translations[edge.trans_index].pos;
            switch(edge.type)
            {
                case 'C': //interaction decoupling n_a*n_b -> n_a*<n_b>+<n_a>*n_b
                    {
                        unsigned int siteA = edge.a;
                        unsigned int siteB = edge.b;
                        
                        ret(edge.a,edge.a) += H->graph.params[edge.ind] * H->densities[siteB];
                        ret(edge.b,edge.b) += H->graph.params[edge.ind] * H->densities[siteA];
                        if(nambu)
                        {
                            ret(edge.a+H->numsites,edge.a+H->numsites) += H->graph.params[edge.ind] * H->densities[siteB];
                            ret(edge.b+H->numsites,edge.b+H->numsites) += H->graph.params[edge.ind] * H->densities[siteA];
                        }
                    }
                    break;
                case 'F':
                    {
                        
                    }
                    break;
                case 'D': //interaction n_a*c_b+n_a*c^dag_b -> <n_a>*c_b+<n_a>*c^dag_b
                    {
                        
                    }
                    break;
                default:
                    std::stringstream ss;
                    ss << "Mean-edge type " << edge.type << " not implemented.";
                    throw std::invalid_argument(ss.str());
                    break;
            }
        }
    
    //bath site stuff
    MatrixXcd gamma;
    if(nambu)
    {
        gamma = MatrixXcd::Zero(H->numsites * 2, H->numsites * 2);
        gamma << get_bath_gamma(omega), MatrixXcd::Zero(H->numsites, H->numsites),
            MatrixXcd::Zero(H->numsites, H->numsites), bfh_sign(get_bath_gamma(-omega));
    }
    else
    {
        gamma = get_bath_gamma(omega);
    }
            
    ret *= -1;
    
    for(int i = 0; i < H->numsites; ++i)
        ret(i,i) += omega;
    if(nambu)
        for(int i = 0; i < H->numsites; ++i)
            ret(i+H->numsites,i+H->numsites) -= omega;
        
    if(H->graph.num_bath > 0)
        return ret - gamma;
    else
        return ret;
}

template MatrixXcd Greens<FLTYPE, FARRAY, MatrixXd>::get_noninteracting_inverse_greens(std::complex<FLTYPE> omega, bool nambu, bool with_mean_field);
template MatrixXcd Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::get_noninteracting_inverse_greens(std::complex<FLTYPE> omega, bool nambu, bool with_mean_field);

/*
 * get bath site gamma. see senechal paper "An introduction to quantum cluster methods"
 */
template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
MatrixXcd Greens<VALTYPE, VECTYPE, MATTYPE>::get_bath_gamma(std::complex<FLTYPE> omega)
{
    VectorXcd diag = VectorXcd::Zero(H->graph.num_bath);
    for(int i = 0; i < H->graph.num_bath; ++i)
        diag[i] = 1.0 / (omega - H->bath_mu[i]);
    return H->bath_hopping*diag.asDiagonal()*H->bath_hopping.adjoint();
}

template MatrixXcd Greens<FLTYPE, FARRAY, MatrixXd>::get_bath_gamma(std::complex<FLTYPE> omega);
template MatrixXcd Greens<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::get_bath_gamma(std::complex<FLTYPE> omega);
