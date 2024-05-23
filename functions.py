import numpy as np
import sympy as sp
import matplotlib.pyplot as plt


def Ce(E,v):
    mu = E/(2*(1+v))
    lambda_ = (v*E)/((1+v)*(1-2*v))
    C_e = np.array([[2*mu+lambda_, lambda_, lambda_],[lambda_, 2*mu+lambda_, lambda_],[lambda_, lambda_, 2*mu+lambda_]])

    return C_e

#For calculating True-strain

def true_strain(nominal_strain):
    true_strain_calc = np.zeros_like(nominal_strain)
    for j in range(3):
        for k in range(3):
            if nominal_strain[j,k] >=0:
                true_strain_calc[j,k] = np.log(1+(nominal_strain[j,k]))
            else:
                true_strain_calc[j,k] = -np.log(1-(nominal_strain[j,k]))
    return true_strain_calc

#For calculating von-mises stress

def von_mises_stress(stress):
    von_mises = (0.5*((stress[0][0]-stress[1][1])**2+(stress[1][1]-stress[2][2])**2+(stress[2][2]-stress[0][0])**2))**0.5
    return von_mises

def von_mises_strain(strain):
    von_mises = (0.5*((strain[0][0]-strain[1][1])**2+(strain[1][1]-strain[2][2])**2+(strain[2][2]-strain[0][0])**2))**0.5
    if strain[0][0]+strain[1][1]+strain[2][2]>=0:
        von_mises = von_mises
    else:
        von_mises = -von_mises
    return von_mises

def von_mises_stress_sign(stress):
    von_mises = (0.5*((stress[0][0]-stress[1][1])**2+(stress[1][1]-stress[2][2])**2+(stress[2][2]-stress[0][0])**2))**0.5
    if stress[0][0]+stress[1][1]+stress[2][2]>=0:
        von_mises = von_mises
    else:
        von_mises = -von_mises
    return von_mises
#For calculating n and m

def nandm(tau_c):
    tau_1, tau_2, tau_3, syield, k = sp.symbols('tau_1 tau_2 tau_3 syield k', real = True)
    von = (0.5*((tau_1-tau_2)**2+(tau_2-tau_3)**2+(tau_3-tau_1)**2))**0.5
    #von = (0.5*((tau_1 - tau_2)**2 + tau_2**2 + tau_1**2))**0.5
    f = von - syield - k
    n = sp.Matrix(np.array([[sp.diff(f,tau_1),0,0],[0,sp.diff(f,tau_2),0],[0,0,sp.diff(f,tau_3)]]))
    h = -sp.diff(f,k)
    n_c = n.subs([(tau_1,tau_c[0][0]),(tau_2,tau_c[1][1]),(tau_3,tau_c[2][2])])
    m_c = n.subs([(tau_1,tau_c[0][0]),(tau_2,tau_c[1][1]),(tau_3,tau_c[2][2])])
    return n_c, m_c, h

def calculate_tau_c(tau_e, current_yield):
    x = np.arange(0,1.01,0.01)
    for i in x:
        tau_e_test = i * tau_e
        von_mises = von_mises_stress(tau_e_test)
        if von_mises < current_yield:
            continue
        else:
            tau_c = tau_e_test
            break
    return tau_c

