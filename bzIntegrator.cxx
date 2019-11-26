/*
This code is part of DIETA

Authored by Kirill Alpin
*/

#include "bzIntegrator.h" 

/*
 * used to integrate an arbitray function over the Brillouin Zone
 */
template<class ITYPE>
ITYPE bzIntFunction(double eta, double nu, void* data)
{
    BZData<ITYPE>* sdata = (BZData<ITYPE>*)data;
    double x = 0.0;
    double y = 0.0;
    for(int i = 0; i < 4; ++i)
    {
        x += sdata->cpt->greens->get_H()->graph.poly_bz[i + sdata->poly_bz_start][0]*BZIntegrator::nodal(eta, nu, i);
        y += sdata->cpt->greens->get_H()->graph.poly_bz[i + sdata->poly_bz_start][1]*BZIntegrator::nodal(eta, nu, i);
    }
    return sdata->f(x, y, sdata->data) * 
        BZIntegrator::jacobian(eta, nu, sdata->cpt->greens->get_H()->graph.poly_bz, sdata->poly_bz_start);
}


/*
 * integrates an arbitray function (f) over the Brillouin Zone
 * (init) is the initial value of the summation, which is usually set to zero/a zero vector
 * if (parallel) is true, the integration is carried out in parallel.
 * in this case the user has to assure, that no race conditions are present in the function (f)
 */
template<class ITYPE>
ITYPE BZIntegrator::bzIntegrate(CPT* cpt, ITYPE (*f)(double,double,void*), void* data, ITYPE init, bool parallel)
{
    const double pi = std::acos(-1);
    ITYPE val = init;
    if(!cpt->H->graph.is_onedim)
    {
        if(cpt->H->graph.poly_bz.size() != 0)
            for(int poly_bz_start = 0; poly_bz_start < cpt->H->graph.poly_bz.size(); poly_bz_start+=4)
            {
                BZData<ITYPE>* bzdata = new BZData<ITYPE>();
                bzdata->data = data;
                bzdata->cpt = cpt;
                bzdata->f = f;
                bzdata->poly_bz_start = poly_bz_start;
                double area = ((cpt->greens->get_H()->graph.poly_bz[1 + poly_bz_start] - 
                    cpt->greens->get_H()->graph.poly_bz[0 + poly_bz_start]).cross(
                    cpt->greens->get_H()->graph.poly_bz[3 + poly_bz_start] - 
                    cpt->greens->get_H()->graph.poly_bz[0 + poly_bz_start])).norm() * 0.5 +
                    ((cpt->greens->get_H()->graph.poly_bz[0 + poly_bz_start] - 
                    cpt->greens->get_H()->graph.poly_bz[2 + poly_bz_start]).cross(
                    cpt->greens->get_H()->graph.poly_bz[3 + poly_bz_start] - 
                    cpt->greens->get_H()->graph.poly_bz[2 + poly_bz_start])).norm() * 0.5;
                
                val += (gauss_legendre_2D_cubeI<ITYPE>(32, bzIntFunction, bzdata, -1.0, 1.0, -1.0, 1.0, init, parallel) * 
                    cpt->greens->get_H()->graph.poly_weight[poly_bz_start/4]) / std::abs(area);
            }
        else
            val = gauss_legendre_2D_cubeI<ITYPE>(32, f, data, -pi, pi, -pi, pi, init, parallel) / (4.0*pi*pi);
    }
    else
        val = gauss_legendreI2<ITYPE>(128, f, data, -pi, pi, init, parallel) / (2.0*pi);
    return val;
}

template double BZIntegrator::bzIntegrate(CPT* cpt, double (*f)(double,double,void*), void* data, double init, bool parallel);
template std::complex<double> BZIntegrator::bzIntegrate(CPT* cpt, std::complex<double> (*f)(double,double,void*), void* data, std::complex<double> init, bool parallel);
template VectorXd BZIntegrator::bzIntegrate(CPT* cpt, VectorXd (*f)(double,double,void*), void* data, VectorXd init, bool parallel);
template VectorXcd BZIntegrator::bzIntegrate(CPT* cpt, VectorXcd (*f)(double,double,void*), void* data, VectorXcd init, bool parallel);
template MatrixXd BZIntegrator::bzIntegrate(CPT* cpt, MatrixXd (*f)(double,double,void*), void* data, MatrixXd init, bool parallel);
template MatrixXcd BZIntegrator::bzIntegrate(CPT* cpt, MatrixXcd (*f)(double,double,void*), void* data, MatrixXcd init, bool parallel);
