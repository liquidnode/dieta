#This code is part of DIETA
#
#Authored by Kirill Alpin

import numpy as np
from dieta_wrapper import DIETA
import numdifftools as nd

USE_3D_OPTIMIZATION = False

#modified from https://codereview.stackexchange.com/questions/38436/python-class-that-implements-the-newton-method
class multivariate_newton(object):

    def __init__(self,func,start_point,step_size=0.5,num_iter=100,tol=1e-5):
        '''
        func: function to be optimized. Takes a vector argument as input and returns
              a scalar output
        step_size: step size in newton method update step
        num_iter: number of iterations for newton method to run
        tol: tolerance to determine convergence
        '''
        self.func=func
        self.start_point=np.array(start_point)
        self.num_iter=num_iter
        self.step_size=step_size
        self.tol=tol
        self.notfound = False


    def optimize(self):
        '''
        perform multivariate newton method for function with vector input
        and scalar output
        '''
        self.notfound = False
        x_t=self.start_point
        #Get an approximation to hessian of function
        H=nd.Hessian(self.func, order=4, step=1e-6)
        #Get an approximation of Gradient of function
        g=nd.Gradient(self.func, order=4, step=1e-6)

        for i in range(self.num_iter):
            gr = g(x_t)
            x_tplus1=x_t-self.step_size*np.dot(np.linalg.inv(H(x_t)),gr)
            #check for convergence
            if max(abs(x_tplus1-x_t))<self.tol:
                break
            x_t=x_tplus1
            print("step\t"+str(np.sqrt(gr.dot(gr))))
            print(x_t)

        self.x=x_tplus1
        self.max_min=self.func(x_t)

        return self

    def critical_point(self):
        '''
        print critical point found in newton_method function. newton_method function
        must be called first.
        '''
        print(self.x)
        
        
mu = -0.5
t =  -0.2
U = 1.0

params = {
        "mu0" : mu,
        "mup0" : 0.0,
        "t" : t,
        "F0" : 0.0,
        "inter" : U/2.0,
        "tap" : 0.0,
        "tacp" : 0.0
        }

dieta = DIETA("graph_2site.txt", params)


def obj_funct(x):
    if not USE_3D_OPTIMIZATION:
        x = np.array([x[0], x[1], 0.0])
    params = {
        "mu0" : mu+x[0],
        "mup0" : -x[0],
        "t" : t,
        "F0" : x[1],
        "inter" : U/2.0,
        "tap" : x[2],
        "tacp" : -x[2]
        }
    
    dieta.update_params(params)
    
    return dieta.get_BSFT_functional()


def get_single_op_densities(x):
    if not USE_3D_OPTIMIZATION:
        x = np.array([x[0], x[1], 0.0])
    params = {
        "mu0" : mu+x[0],
        "mup0" : -x[0],
        "t" : t,
        "F0" : x[1],
        "inter" : U/2.0,
        "tap" : x[2],
        "tacp" : -x[2]
        }
        
    
    dieta.update_params(params)

    return np.array(dieta.get_single_op_densities()).mean()

if USE_3D_OPTIMIZATION:
    x0 = [-0.07023824,  0.01907848, -0.13637898]
else:
    x0 = [-0.10381032,  0.09316732]
fl = open("sweep.txt", "w")
while t <= -0.1:
    res = multivariate_newton(obj_funct, x0)
    res.optimize()
    
    superfluid = get_single_op_densities(res.x)
    print("t\t"+str(t))
    print("Superfluid\t"+str(superfluid))
    fl.write(str(t)+"\t")
    fl.write(str(superfluid)+"\n")
    
    x0 = res.x
    t+=0.005



