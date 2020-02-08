# -*- coding: utf-8 -*-

import numpy as np

# =============================================================================
# Propriedades do ar a uma atmosfera. Referência [001] - F. J. McQuillan, J.R. 
# Culham and M. M. Yovanovich (1984) - Properties of Dry Air at One Atmosphere.
# =============================================================================

# =============================================================================
# Massa específica. Input em K, output em kg/m^3.
# =============================================================================

def rho(T):
    
    a = 351.99/(T)
    
    b = 344.84/(T**2)
    
    RHO = a + b
    
    return(RHO)

# =============================================================================
# Viscosidade dinâmica. Input em K, output em (N.s)/m^2. Referência [002] -
# Reid, Robert C., Sherwood, Thomas K., (1966) - The Properties of Gases and
# Liquids, McGraw-Hill.
# =============================================================================

def mi(T):
    
    a = 1.4592*(T**1.5)
    
    b = 109.10 + T
    
    MI = (a/b)*(10**(-6))
    
    return(MI)
    
# =============================================================================
# Condutividade térmica. Input em K, outpit em W/(m.K). Referência [002] -
# Reid, RobertC., Sherwood, ThomasK., (1966) - The Properties of Gases and
# Liquids, McGraw-Hill.
# =============================================================================

def lamed(T):
    
    a = 2.3340*(T**1.5)*(10**(-3))
    
    b = 164.54 + T
    
    LAMED = a/b
    
    return(LAMED)

# =============================================================================
# Calor específico. Input em K, output em J/(kg.K).
# =============================================================================
    
def cp(T):
    
    a = 1030.5
    
    b = 0.19975*T
    
    c = 3.9734*(T**2)*(10**(-4))
    
    CP = a - b + c
    
    return(CP)
    
# =============================================================================
# Viscosidade Cinemática. Input em K, output em m^2/s.
# =============================================================================
    
def nu(T):
    
    a = 2.4090*(10**8)/(T**1.5)
    
    b = 2.6737*(10**10)/(T**2.5)
    
    c = a + b
    
    NU = c**(-1)
    
    return(NU)
  
# =============================================================================
# Aceleração gravitacional vezes coeficiente de expansão sobre o produto entre
# a viscosidade cinemática e a difusividade térmica; presente no número de 
# Rayleigh (Ra). Input em K, output em 1/(m^3.K).
# =============================================================================

def psi(T):
    
    a = 6.8568*(10**(-3))
    
    b = 1.5079*(10**(-4))*T
    
    c = 1.5715*(10**(-6))*(T**2)
    
    d = (a - b + c)**2
    
    PSI = (10**6)/d
    
    return(PSI)

# =============================================================================
# Difusividade térmica. Input em K, output em m^2/s.
# =============================================================================

def alpha(T):
    
    a = 4.3274
    
    b = 4.1190*(10**(-2))*T
    
    c = 1.5556*(10**(-4))*(T**2)
    
    d = - a + b + c
    
    ALPHA = d*(10**(-6))
    
    return(ALPHA)

# =============================================================================
# Número de Prandtl. Input em K, output em 1.
# =============================================================================
    
def Pr(T):
    
    a = nu(T)
    
    b = alpha(T)
    
    PR = a/b
    
    return(PR)

# =============================================================================
# Lista de propriedades do ar. Input em K, output em unidades previamente
# definidas.
# Air = [alpha, lambda, mi, nu, rho, psi, cp, Pr, T]
# =============================================================================
    
def Air(T):
    
    Air = np.ones(9)
    
    Air[0] = alpha(T)
    
    Air[1] = lamed(T)
    
    Air[2] = mi(T)
    
    Air[3] = nu(T)
    
    Air[4] = rho(T)
    
    Air[5] = psi(T)
    
    Air[6] = cp(T)
    
    Air[7] = Pr(T)
    
    Air[8] = T
    
    return(Air)