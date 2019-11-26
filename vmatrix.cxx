/*
This code is part of DIETA

Authored by Kirill Alpin
*/

#include <iostream>
#include "vmatrix.h"

VMatrix::VMatrix(Hamiltonian* _H)
{
    H = _H;
}

/*
 * construct the V matrix
 */
MatrixXcd VMatrix::getV(Vector3d k)
{
    MatrixXcd ret = MatrixXcd::Zero(H->numsites, H->numsites);
    for(auto const& edge: H->graph.inter_edges)
    {
        const std::complex<FLTYPE> i(0, 1);
        if(edge.trans_index > H->graph.translations.size())
        {
            std::stringstream ss;
            ss << "Invalid translation index.";
            throw std::invalid_argument(ss.str());
        }
        Vector3d translation = H->graph.translations[edge.trans_index].pos;
        switch(edge.type)
        {
            case 'A': //hopping term c^dag*c+c*c^dag
                {
                    ret(edge.b,edge.a) += H->graph.params[edge.ind] * std::exp(-i*(k.dot(translation)));// * std::sqrt(1.0 - H->densities[edge.a]) * std::sqrt(1.0 - H->densities[edge.b]);
                    if(edge.a != edge.b || translation.cwiseAbs().maxCoeff() > 1e-10)
                        ret(edge.a,edge.b) += std::conj(H->graph.params[edge.ind]) * std::exp(i*(k.dot(translation)));//  * std::sqrt(1.0 - H->densities[edge.a]) * std::sqrt(1.0 - H->densities[edge.b]);
                }
                break;
            case 'B': //hopping term c*c+c^dag*c^dag
                {
                    std::stringstream ss;
                    ss << "Edge type B not possible with non nambu G.";
                    throw std::invalid_argument(ss.str());
                }
                break;
            case 'F':
                {
                    ret(edge.a,edge.a) += H->graph.params[edge.ind];
                }
                break;
            case 'G':
                {
                    ret(edge.a,edge.a) += H->graph.params[edge.ind];
                }
                break;
            default:
                std::stringstream ss;
                ss << "Edge type " << edge.type << " not implemented.";
                throw std::invalid_argument(ss.str());
                break;
        }
    }
    
    return ret;
}

/*
 * construct the V matrix using the Nambu formalism
 */
MatrixXcd VMatrix::get_nambu_V(Vector3d k)
{
    MatrixXcd ret = MatrixXcd::Zero(H->numsites * 2, H->numsites * 2);
    for(auto const& edge: H->graph.inter_edges)
    {
        const std::complex<FLTYPE> i(0, 1);
        if(edge.trans_index > H->graph.translations.size())
        {
            std::stringstream ss;
            ss << "Invalid translation index.";
            throw std::invalid_argument(ss.str());
        }
        Vector3d translation = H->graph.translations[edge.trans_index].pos;
        switch(edge.type)
        {
            case 'A': //hopping term c^dag*c+c*c^dag. if edge.a==edge.b c^dag*c i.e. chemical potential
                {
                    ret(edge.b,edge.a) += H->graph.params[edge.ind] * std::exp(-i*(k.dot(translation)));
                    if(edge.a != edge.b || translation.cwiseAbs().maxCoeff() > 1e-10)
                    {
                        ret(edge.a,edge.b) += std::conj(H->graph.params[edge.ind]) * std::exp(i*(k.dot(translation)));
                    }
                    ret(edge.b+H->numsites,edge.a+H->numsites) += H->graph.params[edge.ind] * std::exp(-i*(k.dot(translation)));
                    if(edge.a != edge.b || translation.cwiseAbs().maxCoeff() > 1e-10)
                        ret(edge.a+H->numsites,edge.b+H->numsites) += std::conj(H->graph.params[edge.ind]) * std::exp(i*(k.dot(translation)));
                }
                break;
            case 'B': //anomalous term c*c+c^dag*c^dag
                {
                    ret(edge.b,edge.a+H->numsites) += H->graph.params[edge.ind] * std::exp(-i*(k.dot(translation)));
                    if(edge.a != edge.b || translation.cwiseAbs().maxCoeff() > 1e-10)
                        ret(edge.a,edge.b+H->numsites) += std::conj(H->graph.params[edge.ind]) * std::exp(i*(k.dot(translation)));
                    ret(edge.b+H->numsites,edge.a) += H->graph.params[edge.ind] * std::exp(-i*(k.dot(translation)));
                    if(edge.a != edge.b || translation.cwiseAbs().maxCoeff() > 1e-10)
                        ret(edge.a+H->numsites,edge.b) += std::conj(H->graph.params[edge.ind]) * std::exp(i*(k.dot(translation)));
                }
                break;
            case 'F':
                {
                    ret(edge.a,edge.a) += H->graph.params[edge.ind];
                    ret(edge.a+H->numsites,edge.a+H->numsites) += H->graph.params[edge.ind];
                }
                break;
            case 'G':
                {
                    ret(edge.a,edge.a) += H->graph.params[edge.ind];
                }
                break;
            case 'H':
                {
                    ret(edge.a+H->numsites,edge.a+H->numsites) += H->graph.params[edge.ind];
                }
                break;
            default:
                std::stringstream ss;
                ss << "Edge type " << edge.type << " not implemented.";
                throw std::invalid_argument(ss.str());
                break;
        }
    }
    
    return ret;
}
