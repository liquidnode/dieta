/*
This code is part of DIETA

Authored by Kirill Alpin
*/

#pragma once
#include <Eigen/Dense>
#include "hamiltonian.h"

using Eigen::MatrixXcd;
using Eigen::Vector3d;

/*
 * this class holds the V matrix of the system
 */
class VMatrix
{
public:
    Hamiltonian* H;
    VMatrix(Hamiltonian* _H);
    MatrixXcd getV(Vector3d k);
    MatrixXcd get_nambu_V(Vector3d k);
};
