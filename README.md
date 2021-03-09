# Study-of-the-quantum-Hall-effect-of-the-exxciton-valley-interlaced-in-the-new-two-dimensional-materi
pru pru pru
import numpy as np
import matplotlib.pyplot as plt

#define os valores dos parâmetros

hbar = 1/0.658212 #Planck constant's value
delta1 = 0.0 #mueV
omega = 50.0 #mueV
delta2 = delta1 + 2*omega #mueV
T_e = 10.0 #mueV
Omega = 50.0

#define a função 
def func(P, t):
    p_00, p_01, p_02, p_10, p_11, p_12, p_20, p_21, p_22 = P
    return np.array([-1j * hbar * Omega * (p_10 - p_01),
                     1j * hbar * (delta1 * p_01 + Omega * (p_00 - p_11)+ T_e * p_02),
                     1j * hbar * (p_02/2* (delta1 + delta2) - Omega * p_12 + T_e * p_01),
                     1j * hbar * (- delta1 * p_10 + Omega * (p_11 - p_00) - T_e * p_20), 
                     1j * hbar * (Omega * (p_10 - p_01) + T_e * (p_12 - p_21)),
                     1j * hbar * (p_12/2* (delta2 - delta1) - Omega * p_02 + T_e *(p_11 - p_22)),
                     1j * hbar * (-p_02/2* (delta2 + delta1) + Omega * p_21 - T_e * p_10),
                     1j * hbar * (p_21/2* (delta1 - delta2) + Omega * p_20 + T_e * (p_22 - p_11)), 
                     1j * T_e * hbar * (p_21 - p_12)], dtype = "complex")

#define o runge-kutta de quarta ordem
def rk4(f, P0, t0, tf , n):
    t = np.linspace(t0, tf, n+1)
    # P = np.array((n+1)*[P0])
    P = np.array((n+1)*[P0], dtype = "complex")     #teste
    h = t[1]-t[0]
    for i in range(n):
        k1 = f(P[i], t[i])    
        k2 = f(P[i] + 0.5*h*k1, t[i] + 0.5*h)
        k3 = f(P[i] + 0.5*h*k2, t[i] + h*0.5)
        k4 = f(P[i] + k3*h, t[i] + h)

#----------PLOT-----------------------

#plot the graph
P, t  = rk4(func, np.array([1., 0., 0., 0., 0., 0., 0., 0., 0.], dtype = "complex") , 0. , 1, 10000)
p_00, p_01, p_02, p_10, p_11, p_12, p_20, p_21, p_22 = P.T
plt.plot(t, p_00.real) #'-',color='navy',linewidth=1.5,label = 'omega_{12} = 0')

#Axis range
plt.axis([0.0, 1.0, 0.0 ,1 ])
plt.grid('on')
plt.show()
  

