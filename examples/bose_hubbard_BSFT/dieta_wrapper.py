#This code is part of DIETA
#
#Authored by Kirill Alpin

import numpy as np
import subprocess
import sys
from scipy.optimize import minimize

class DIETA:
    def __init__(self, graph_file, params):
        self.graph_file = graph_file
        self.params = params
        self.instable = False
        
    def write_params(self):
        lines_out = []  
        
        for k in self.params:
            lines_out.append(k+"\t"+str(self.params[k])+"\n")
    
        f = open("params.txt", "w")
        for c in lines_out:
            f.write(c)
        f.close()
        
    def update_params(self, params):
        self.params = params
        
    def get_BSFT_functional(self):
        self.write_params()
        proc = subprocess.Popen(['../../dieta',self.graph_file,'params.txt','-gp'],stdout=subprocess.PIPE)
        
        (output, err) = proc.communicate()
        p_status = proc.wait()
        output = output.decode('ascii')

        self.instable = output.splitlines()[-2].startswith("Complex eigenvalue")
        
        return float(output.splitlines()[-1].split('\t')[-1])
    
    def get_single_op_densities(self):
        self.write_params()
        proc = subprocess.Popen(['../../dieta',self.graph_file,'params.txt','-sp'],stdout=subprocess.PIPE)
        
        (output, err) = proc.communicate()
        p_status = proc.wait()
        output = output.decode('ascii')

        return [float(a) for a in output.splitlines()[-1].split('\t')]
