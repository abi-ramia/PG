# -*- coding: utf-8 -*-

import math
import numpy as np
from scipy import optimize

# =============================================================================
# Propriedades da Água Líquida. Referências [017], [019], [020], [021] e [022].
# =============================================================================

#Massa específica. Entrada em K, saída em kg/m^3.
def rho(T):
    
    t = T - 273.15
    
    a = +999.83952
    b = +16.945176*t
    c = -7.9870401*(10**(-3))*(t**2)
    d = -46.170461*(10**(-6))*(t**3)
    e = +105.56302*(10**(-9))*(t**4)
    f = -280.54253*(10**(-12))*(t**5)
    g = +1
    h = +16.87985*(10**(-3))*t
    
    rho = (a+b+c+d+e+f)/(g+h)
    
    return (rho)

#Viscosidade dinâmica. Entrada em K, saída em Pa.s.
def mu(T):
    
    t = T - 273.15
    
    a = +1002/(10**6)
    b = +20 - t
    c = +96 + t
    d = +1.2378
    e = -1.303*(10**(-3))*((20 - t))
    f = +3.060*(10**(-6))*((20 - t)**2)
    g = +2.550*(10**(-8))*((20 - t)**3)
    
    mu = a*(10**(((b)*(d+e+f+g))/(c)))

    return (mu)

#Viscosidade cinemática. Entrada em K, saída em m^2/s.
def nu(T):
    
    a = mu(T)
    b = rho(T)
    
    nu = a/b
    
    return (nu)

#Condutividade térmica. Entrada em K, saída em W/(m.K).
def lamed(T):
    
    T_aster = T/298.15
    
    a = +0.60650
    b = -1.48445
    c = +4.12292*((T_aster))
    d = -1.63866*((T_aster)**2)
    
    lamed = a*(b+c+d)
    
    return (lamed)

#Calor específico. Entrada em K, saída em J/(kg.K).
def cp(T):
    
    t = T - 273.15
    
    a = 4185.5
    b = 0.996185
    c = 0.0002874*(((t + 100)/(100))**5.26)
    d = 0.011160*(10**(-0.036*t))
    
    cp = a*(b+c+d)
    
    return (cp)

#Número de Reynolds. Entrada em kg/s, m e K, saída em 1.
def Re(m, D, T):
    
    A = math.pi*(D**2)/4
    r = rho(T)
    U = (m/r)/A
    R = U*D/nu(T)
    
    return (R)

#Número de Prandtl. Entrada em K, saída em 1.
def Pr(T):
    
    a = cp(T)*mu(T)
    b = lamed(T)
    
    return (a/b)

#Fator de atrito. Entrada em kg/s, m e K, saída em 1.
def FA(m, D, T, e):
    
    ReD = Re(m, D, T)
    
    def errf(f):
        
        err = (1/math.sqrt(f)) + 2*math.log10(((e/D)/3.7) + 2.51/(ReD*math.sqrt(f)))
        
        return (err)
        
    root = optimize.brentq(errf, 0.001, 1)
    
    f = root
    
    return (f)

#Coeficiente de transferência de calor. Entrada em kg/s, m, K e m, saída em W/(m^2.K).
def hc(m, D, T, e):
    
    R = Re(m, D, T)
    P = Pr(T)
    f = FA(m, D, T, e)
    
    a = (f/8)*(R - 1000)*P
    b = 1 + 12.7*(math.sqrt(f/8))*((P**(2/3)) - 1)
    
    NuD = a/b
    
    k = lamed(T)
    
    h = NuD*k/D
    
    return(h)

# =============================================================================
# Ponto de Orvalho da Água. Referências [025] e [026].
# =============================================================================

#Calcula a temperatura do ponto de orvalho. Entrada em K e 1, saída em K.
def TDP(T, RH):
    
    t = T - 273.15
    
    b = 17.368
    c = 238.88
    d = 234.5
    
    gamma = np.log(RH*np.exp((b - (t/d))*(t/(c + t))))
    
    dp = c*gamma/(b - gamma)
    
    return(dp + 273.15)

def fl(Te, Ta, RH):
    
    if (Te) < (TDP(Ta, RH)):
        return ("Sim")
    return ("Não")