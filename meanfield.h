/*
This code is part of DIETA

Authored by Kirill Alpin
*/

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <Eigen/Dense>
#include "parser.h"
#include "hamiltonian.h"
#include "lanczos.h"
#include "greens.h"
#include "cpt.h"
#include "density.h"

/*
 * this class is used to perform mean field decoupling
 */
class MeanField
{
public:
    CPT* cpt;
    MeanField(CPT* _cpt);
    void findFixpoint();
    void findFixpointWithCurrent();
};
