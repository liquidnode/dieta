/*
This code is part of DIETA

Authored by Kirill Alpin
*/

#include "lanczos.h"

template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
Lanczos<VALTYPE, VECTYPE, MATTYPE>::Lanczos(Hamiltonian* hamiltonian)
{
    H = hamiltonian;
}

template Lanczos<FLTYPE, FARRAY, MatrixXd>::Lanczos(Hamiltonian* hamiltonian);
template Lanczos<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::Lanczos(Hamiltonian* hamiltonian);

#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))
std::vector<FLTYPE> tqli2(std::vector<FLTYPE> d, std::vector<FLTYPE> e, int n)
/***
 * from https://github.com/rgmelko
 * original from April 2005, Roger Melko, modified from Numerical Recipies in C v.2
 * modified from www.df.unipi.it/~moruzzi/
 ***/
{
    
    int m,l,iter,i,k;
    FLTYPE s,r,p,g,f,dd,c,b;

    for ( l=0; l < n; l++) 
    {
        iter = 0;
        do 
        { 
            for (m = l; m < n-1; m++) 
            { 
                dd = fabs( d[m] ) + fabs( d[m+1] );
                if (fabs( e[m] ) + dd == dd) break;
            }
            if (m != l) { 
                if (iter++ == 30) { 
                    std::stringstream ss;
                    ss << "Too many iterations in tqli().";
                    throw std::invalid_argument(ss.str());
                    return std::vector<FLTYPE>();
                }
                g = ( d[l+1] - d[l] )/( 2.0 * e[l] );
                r = sqrt( (g * g) + 1.0 );
                g = d[m] - d[l] + e[l] / ( g + SIGN(r,g) );
                s = c = 1.0;
                p = 0.0;
                for ( i = m - 1; i >= l; i--) 
                { 
                    f = s * e[i];
                    b = c * e[i];
                    if ( fabs( f ) >= fabs( g ) ) 
                    { 
                        c = g / f;
                        r = sqrt( (c*c) + 1.0 );
                        e[i + 1] = f * r;
                        c *= ( s = 1.0 / r );
                    }
                    else 
                    { 
                        s = f / g;
                        r = sqrt( (s*s) + 1.0 );
                        e[i + 1] = g * r;
                        s *= ( c = 1.0 / r );
                    }
                    g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    p = s * r;
                    d[i + 1] = g + p;
                    g = c * r - b;
                    /*EVECTS*/
                    /*if (Evects == 1) 
                    {
                        for ( k = 0; k < n; k++ ) 
                        { 
                            f = z[k][i + 1];
                            z[k][i + 1] = s * z[k][i] + c * f;
                            z[k][i] = c * z[k][i] - s * f;
                        }      
                    }*///Evects
                }
                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        } while (m != l);
    }
    return d;
}

/*
 * constructs a tridiagonal Matrix representation of the Hamiltonian using a Lanczos iteration
 */
template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
MATTYPE Lanczos<VALTYPE, VECTYPE, MATTYPE>::get_tridiagonal(int seed)
{
    std::cout << "Start get tridiagonal. ";
    H->set_particle_space(H->particle_number);
    
    bsqr = std::vector<FLTYPE>();
    a = std::vector<VALTYPE>();
    int numIter = 500;
    VECTYPE states[3];
    for(int i = 0; i < 3; ++i)
        states[i] = VECTYPE::Zero(H->dim);
    
    //random init state
    std::mt19937 gen(seed);
    std::normal_distribution<> dis(0.0, 1.0);
    for(int i = 0; i < H->dim; ++i)
        states[0][i] = dis(gen);
    states[0] /= std::sqrt(states[0].dot(states[0]));
    
    states[1] = (*H) * states[0];
    a.push_back(states[0].dot(states[1])/states[0].dot(states[0]));
    
    states[1] = states[1] - a[0]*states[0];
    
    //iterate lanczos
    bool nextlvl = false;
    int iterated = numIter;
    FLTYPE eig3 = 1e10;
    for(int i = 0; i < numIter; ++i)
    {
        int nm = i % 3;
        int n = (i + 1) % 3;
        int np = (i + 2) % 3;
        
        FLTYPE norm_nm = std::real(states[nm].dot(states[nm]));
        FLTYPE norm_n = std::real(states[n].dot(states[n]));
        
        bsqr.push_back(norm_n / norm_nm);
        
        states[nm] /= std::sqrt(norm_n);
        states[n] /= std::sqrt(norm_n);
        
        states[np] = (*H) * states[n];
        a.push_back(states[n].dot(states[np]));
        
        if(i != 0 && i > 10)
        {
            //get d and e
            std::vector<FLTYPE> d;
            std::vector<FLTYPE> e;
            int n = a.size();
            e.push_back(0);
            for(int i = 0; i < bsqr.size(); ++i)
                e.push_back(std::sqrt(bsqr[i]));
            if(typeid(VALTYPE) == typeid(FLTYPE))
                for(int i = 0; i < a.size(); ++i)
                    d.push_back(std::real(a[i]));
            else
                for(int i = 0; i < a.size(); ++i)
                    d.push_back(std::abs(a[i]));
            
            std::vector<FLTYPE> eigs = tqli2(d, e, n);
            std::sort(eigs.begin(), eigs.end());
            
            if(std::abs(eigs[0] - eig3) < 1e-10)
            {
                iterated = i;
                break;
            }
            eig3 = eigs[0];
        }
        
        states[np] = states[np] - a[i+1]*states[n] - bsqr[i]*states[nm];
    }
    
    std::cout << "Stop " << iterated << std::endl;
    
    //build tridiagonal matrix
    MATTYPE M = MATTYPE::Zero(a.size(), a.size());
    
    for(int i = 0; i < a.size(); ++i)
        M(i, i) = a[i];
    for(int i = 0; i < bsqr.size(); ++i)
    {
        M(i+1,i) = std::sqrt(bsqr[i]);
        M(i,i+1) = std::sqrt(bsqr[i]);
    }
    
    return M;
}

template MatrixXd Lanczos<FLTYPE, FARRAY, MatrixXd>::get_tridiagonal(int seed);
template MatrixXcd Lanczos<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::get_tridiagonal(int seed);

/*
 * returns the groundstate given a tridiagonal representation of the Hamiltonian
 */
template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
VECTYPE Lanczos<VALTYPE, VECTYPE, MATTYPE>::get_groundstate(int seed, MATTYPE tridiagonal, FLTYPE &genergy)
{
    std::cout << "Start get groundstate" << std::endl;
    H->set_particle_space(H->particle_number);
    
    SelfAdjointEigenSolver<MATTYPE> es(tridiagonal);
    energies = es.eigenvalues();
    evecs = es.eigenvectors();
    
    int sm = 0;
    genergy = std::real(energies[sm]);
    if(energies.size() > 2)
        std::cout << "Cluster energies: " << energies[0] << "\t" << energies[1] << "\t" << energies[2] << std::endl;
    for(int i = 0; i < energies.size(); ++i)
        if(std::real(energies[i]) < genergy)
        {
            genergy = std::real(energies[i]);
            sm = i;
        }
        
    VECTYPE gtridiag = evecs.col(sm);
    VECTYPE gstate = VECTYPE::Zero(H->dim);
    
    //reinit lanczos
    int numIter = bsqr.size();
    VECTYPE states[3];
    for(int i = 0; i < 3; ++i)
        states[i] = VECTYPE::Zero(H->dim);
    
    //random init state
    std::mt19937 gen(seed);
    std::normal_distribution<> dis(0.0, 1.0);
    for(int i = 0; i < H->dim; ++i)
        states[0][i] = dis(gen);
    states[0] /= std::sqrt(std::real(states[0].dot(states[0])));
    
    gstate += gtridiag[0] * states[0];
    
    states[1] = (*H) * states[0];
    
    states[1] = states[1] - a[0]*states[0];
    
    //iterate lanczos
    for(int i = 0; i < numIter; ++i)
    {
        int nm = i % 3;
        int n = (i + 1) % 3;
        int np = (i + 2) % 3;
        
        FLTYPE norm_n = std::real(states[n].dot(states[n]));
        
        states[nm] /= std::sqrt(norm_n);
        states[n] /= std::sqrt(norm_n);
        gstate += gtridiag[i+1] * states[n];
        
        if(i != numIter-1)
        {
            states[np] = (*H) * states[n] - a[i+1]*states[n] - bsqr[i]*states[nm];
        }
    }
    
    //TODO noetig??
    //gstate /= std::sqrt(gstate.dot(gstate));
    
    return gstate;
}

template FARRAY Lanczos<FLTYPE, FARRAY, MatrixXd>::get_groundstate(int seed, MatrixXd tridiagonal, FLTYPE &genergy);
template VectorXcd Lanczos<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::get_groundstate(int seed, MatrixXcd tridiagonal, FLTYPE &genergy);

/*
 * legacy code
 * used to construct continued fraction Greens functions
 */
template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
template <typename T, typename U> 
LanczosParams<T> Lanczos<VALTYPE, VECTYPE, MATTYPE>::get_params(U init_state, bool hole)
{
    std::cout << "Start get params.";
    FLTYPE sign = 1.0;
    if(hole)
        sign = -1.0;
    LanczosParams<T> ret;
    
    int numIter = 500;
    U states[3];
    for(int i = 1; i < 3; ++i)
        states[i] = U::Zero(H->dim);
    
    //random init state
    states[0] = init_state;
    
    states[1] = sign * ((*H) * states[0]);
    ret.bsqr.push_back(std::real(states[0].dot(states[0])));
    ret.a.push_back(states[0].dot(states[1])/states[0].dot(states[0]));
    
    states[1] = states[1] - ret.a[0]*states[0];
    
    int iterated = numIter;
    //iterate lanczos
    FLTYPE eig3 = 1e10f;
    for(int i = 0; i < numIter; ++i)
    {
        int nm = i % 3;
        int n = (i + 1) % 3;
        int np = (i + 2) % 3;
        
        FLTYPE norm_nm = std::real(states[nm].dot(states[nm]));
        FLTYPE norm_n = std::real(states[n].dot(states[n]));
        
        ret.bsqr.push_back(norm_n / norm_nm);
        
        states[nm] /= std::sqrt(norm_n);
        states[n] /= std::sqrt(norm_n);
        
        states[np] = sign * ((*H) * states[n]);
        ret.a.push_back(states[n].dot(states[np]));
        
        if(i % 5 == 0 && i > 10)
        {
            //get d and e
            std::vector<FLTYPE> d;
            std::vector<FLTYPE> e;
            int n = ret.a.size();
            for(int i = 0; i < ret.bsqr.size(); ++i)
                e.push_back(std::sqrt(ret.bsqr[i]));
            if(typeid(VALTYPE) == typeid(FLTYPE))
                for(int i = 0; i < ret.a.size(); ++i)
                    d.push_back(std::real(ret.a[i]));
            else
                for(int i = 0; i < ret.a.size(); ++i)
                    d.push_back(std::abs(ret.a[i]));
            
            std::vector<FLTYPE> eigs = tqli2(d, e, n);
            std::sort(eigs.begin(), eigs.end());
            
            if(std::abs(eigs[10] - eig3) < 1e-15)
            {
                iterated = i;
                break;
            }
            eig3 = eigs[10];
        }
        
        states[np] = states[np] - ret.a[i+1]*states[n] - ret.bsqr[i+1]*states[nm];
    }
    std::cout << " Stop " << iterated << std::endl;
    return ret;
}

template LanczosParams<FLTYPE> Lanczos<FLTYPE, FARRAY, MatrixXd>::get_params(FARRAY init_state, bool hole=false);
template LanczosParams<std::complex<FLTYPE>> Lanczos<FLTYPE, FARRAY, MatrixXd>::get_params(VectorXcd init_state, bool hole=false);
template LanczosParams<std::complex<FLTYPE>> Lanczos<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::get_params(VectorXcd init_state, bool hole=false);

/*
 * orthogonolizes the input vectors using Gram Schmidt
 */
template <typename VALTYPE, typename VECTYPE>
void gramSchmidt(std::vector<VECTYPE>& vecs)
{
    for(int j = 0; j < vecs.size(); ++j)
    {
        for(int i = 0; i < j; ++i)
        {
            VALTYPE rij = vecs[i].dot(vecs[j]);
            vecs[j] = vecs[j] - rij*vecs[i];
        }
        FLTYPE rjj = std::sqrt(std::real(vecs[j].dot(vecs[j])));
        if(rjj < 1e-10)
        {
            std::stringstream ss;
            ss << "Basis is not linear independent.";
            throw std::invalid_argument(ss.str());
        }
        else
        {
            vecs[j] /= rjj;
        }
    }
}

template void gramSchmidt<FLTYPE, FARRAY>(std::vector<FARRAY>& vecs);
template void gramSchmidt<std::complex<FLTYPE>, VectorXcd>(std::vector<VectorXcd>& vecs);

template void Lanczos<FLTYPE, FARRAY, MatrixXd>::get_Qmatrix(FARRAY groundstate, MatrixXd& T, MatrixXd& Q_elec, MatrixXd& Q_hole, bool hole, MatrixXd phi_elec, MatrixXd phi_hole);
template void Lanczos<std::complex<FLTYPE>, VectorXcd, MatrixXcd>::get_Qmatrix(VectorXcd groundstate, MatrixXcd& T, MatrixXcd& Q_elec, MatrixXcd& Q_hole, bool hole, MatrixXcd phi_elec, MatrixXcd phi_hole);

/*
 * constructs the Q matrices using the Band Lanczos method
 */
template <typename VALTYPE, typename VECTYPE, typename MATTYPE>
void Lanczos<VALTYPE, VECTYPE, MATTYPE>::get_Qmatrix(VECTYPE groundstate, MATTYPE& T, MATTYPE& Q_elec, MATTYPE& Q_hole, bool hole, MATTYPE phi_elec, MATTYPE phi_hole)
{
    bool reortho = false;
    
    std::cout << "Start get Q and T." << std::endl;
    if(hole)
        H->set_particle_space(H->particle_number - 1);
    else
        H->set_particle_space(H->particle_number + 1);
    
    // particle part
    std::vector<VECTYPE> q;
    int numIter = 600;
    Q_elec = MATTYPE::Zero(H->numsites, numIter + H->numsites + 1);
    Q_hole = MATTYPE::Zero(H->numsites, numIter + H->numsites + 1);
    int ndel = 0;
    for(int i = 0; i < H->numsites; ++i)
    {
        if(hole)
            q.push_back(phi_hole.row(i).adjoint());
        else
            q.push_back(phi_elec.row(i).adjoint());
        //check for linear independence
        VECTYPE testq = VECTYPE::Zero(q[0].size());
        for(int j = ndel; j < i; ++j)
            testq = testq + q[j - ndel].dot(q[i - ndel]) * q[j - ndel];
        testq = testq - q[i - ndel];
        if(std::abs(testq.dot(testq)) < 1e-14)
        {
            q.pop_back();
            ndel++;
        }
        else
        {
            FLTYPE norm = std::sqrt(std::real(q[i-ndel].dot(q[i-ndel])));
            if(norm < 1e-14)
            {
                q.pop_back();
                ndel++;
            }
            else
                q[i-ndel] = q[i-ndel] / norm;
        }
    }
    
    if(q.size() == 0)
    {
        T = MATTYPE::Zero(1, 1);
        Q_elec.conservativeResize(H->numsites, 1);
        Q_hole.conservativeResize(H->numsites, 1);
        std::cout << " Stop " << 0 << std::endl;
        return;
    }
    //q.push_back(groundstate); //braucht man das?
    
    gramSchmidt<VALTYPE, VECTYPE>(q);
    
    int L = q.size();
    int p_c = L + 1;
    
    std::vector<VECTYPE> allv;
    
    
    q.push_back(VECTYPE::Zero(q[0].size()));
    std::vector<VECTYPE> v;
    for(int i = 0; i < L + 1; ++i)
        v.push_back(VECTYPE::Zero(q[0].size()));
    
    //make lanczos run
    T = MATTYPE::Zero(numIter + L, numIter + L);
    
    std::vector<int> I;
    std::vector<VECTYPE> Iv;
    int shift = 0; 
    int realIter = numIter;
    int iterated = 0;
    int eigconv_index = 0;
    FLTYPE lastsm = 1e10;
    bool converged = false;
    int j = 0; 
    while(!converged)
    {
        FLTYPE beta_n = std::sqrt(std::real(q[(j+shift)%(L+1)].dot(q[(j+shift)%(L+1)])));
        if(beta_n <= 1.0e-8)
        {
            std::cout << "DEFLATE. beta: " << beta_n << std::endl;
            if(j + 1 - p_c > 0)
            {
                I.push_back(j-p_c);
                Iv.push_back(v[(j-p_c)%(L+1)]);
            }
            shift++;
            p_c--;
            if(p_c == 0)
            {
                realIter = j;
                break;
            }
            continue;
        }
        if(j - p_c >= 0)
            T(j,j-p_c) = beta_n;
        v[j%(L+1)] = q[(j+shift)%(L+1)] / beta_n;
        
        if(reortho)
        {
            for(int k = 0; k < (int)allv.size() - L - 1; ++k)
                v[j%(L+1)] -= allv[k].dot(v[j%(L+1)]) * allv[k];
            allv.push_back(v[j%(L+1)]);
        }
        
        if(!H->is_particle_conserving)
        {
            Q_elec.col(j) = phi_elec*v[j%(L+1)];
            Q_hole.col(j) = phi_hole*v[j%(L+1)];
        }
        else
        {
            if(hole)
                Q_hole.col(j) = phi_hole*v[j%(L+1)];
            else
                Q_elec.col(j) = phi_elec*v[j%(L+1)];
        }
        
        for(int k = j+1; k < j+p_c; ++k)
        {
            if(k - p_c >= 0)
            {
                T(j,k-p_c) = v[j%(L+1)].dot(q[(k+shift)%(L+1)]);
                q[(k+shift)%(L+1)] = q[(k+shift)%(L+1)] - v[j%(L+1)] * T(j,k-p_c);
            }
        }
        q[(j+p_c+shift)%(L+1)] = (*H) * v[j%(L+1)];
        for(int k = std::max(0,j-p_c); k < j; ++k)
        {
            T(k,j)=std::conj(T(j,k));
            q[(j+p_c+shift)%(L+1)] = q[(j+p_c+shift)%(L+1)] - v[k%(L+1)] * T(k,j);
        }
        std::vector<int> Ip = I;
        Ip.push_back(j);
        std::sort(Ip.begin(), Ip.end());
        for(std::vector<int>::iterator it = Ip.begin(); it != Ip.end(); ++it)
        {
            int k = *it;
            std::vector<int>::iterator ind = std::find(I.begin(), I.end(), k);
            if(ind != I.end())
            {
                int kind = std::distance(I.begin(), ind);
                T(k,j)=Iv[kind].dot(q[(j+p_c+shift)%(L+1)]);
                q[(j+p_c+shift)%(L+1)] = q[(j+p_c+shift)%(L+1)]-Iv[kind]*T(k,j);
            }
            else
            {
                T(k,j)=v[k%(L+1)].dot(q[(j+p_c+shift)%(L+1)]);
                q[(j+p_c+shift)%(L+1)] = q[(j+p_c+shift)%(L+1)]-v[k%(L+1)]*T(k,j);
            }
        }
        for(std::vector<int>::iterator it = I.begin(); it != I.end(); ++it)
            T(j,*it) = std::conj(T(*it,j));
        if(j > eigconv_index && j > L)
        {
            SelfAdjointEigenSolver<MATTYPE> es(T.block(0,0,j+1,j+1));
            auto eigs = es.eigenvalues();
            std::cout << std::abs(eigs[eigconv_index]-lastsm) << "\t" << eigs[eigconv_index] << std::endl;
            if(std::abs(eigs[eigconv_index]-lastsm) < 1.0e-9 && j % 4 == 0)
            {
                realIter = j;
                iterated = j;
                break;
            }
            if(j % 4 == 0)
                lastsm = eigs[eigconv_index];
        }
        if(j == numIter - 1)
        {
            realIter = j;
            iterated = j;
            break;
        }
        j++;
    }
    
    T.conservativeResize(realIter, realIter);
    Q_elec.conservativeResize(H->numsites, realIter);
    Q_hole.conservativeResize(H->numsites, realIter);
    
    std::cout << " Stop " << iterated << std::endl;
}
