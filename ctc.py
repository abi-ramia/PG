# -*- coding: utf-8 -*-

import numpy as np

import ar

# =============================================================================
# Coeficientes de transferência de calor e fluxo convectivo. Referências [016]
# - N-550 - Projeto de Isolamento Térmico a Alta Temperatura, Anexo E; e [017]
# - Incropera - Fundamentals of Heat & Mass Transfer - 8th Edition.
# =============================================================================



# =============================================================================
# Funções auxiliares.
# =============================================================================

# =============================================================================
# Emissividade de superfícies. Input em string, output 1.
# =============================================================================

def emissividade(nome):

    """
    Preencha com um nome entre aspas.
    Opções:
    'Alumínio';
    'Tinta Preta';
    'Tinta de Alumínio';
    'Chapa de Aço';
    'Tinta Branca';
    'Massa Asfáltica'.
    """

    emissividades = {'Alumínio': 0.2,
                     'Tinta Preta': 0.97,
                     'Tinta de Alumínio': 0.5,
                     'Chapa de Aço' : 0.955,
                     'Tinta Branca' : 0.88,
                     'Massa Asfáltica' : 0.93}

    epsilon = emissividades[nome]

    return(epsilon)

# =============================================================================
# Número de Reynolds de um escoamento. Input em m/s, m e m^2/s, output em 1.
# =============================================================================

def Re(U, Lc, Te, Ta):

    Tf = (Te + Ta)/2

    nu = ar.nu(Tf)

    Re = U*Lc/nu

    return(Re)

# =============================================================================
# Número de Rayleigh. Input em m e K, output em 1.
# =============================================================================

def Ra(Lc, Te, Ta):

    Tf = (Te + Ta)/2

    psi = ar.psi(Tf)

    DT = abs(Te - Ta)

    Ra = psi*(Lc**3)*DT

    return(Ra)

# =============================================================================
# Número de Grashof. Input em m e K, output em 1.
# =============================================================================

def Gr(Lc, Te, Ta):
    
    Tf = (Te + Ta)/2
    
    Ray = Ra(Lc, Te, Ta)
    
    Pra = ar.Pr(Tf)
    
    Gra = Ray/Pra
    
    return (Gra)

# =============================================================================
# Número de Richardson. Input em m/s, m e m^2/s, output em 1.
# =============================================================================

def Ri(U, Lc, Te, Ta):
    
    Gra = Gr(Lc, Te, Ta)
    
    Rey = Re(U, Lc, Te, Ta)
    
    Ric = Gra/(Rey**2)
    
    return (Ric)

# =============================================================================
# Radiação.
# =============================================================================

# =============================================================================
# Coeficiente de transferência de calor por radiação. Input em 1 e K, output em
# W/(m^2.K).
# =============================================================================

def hr(epsilon, Te, Ta):

    sigma = 5.670374e-08

    T = ((Te)**2+(Ta)**2)*(Te + Ta)

    hr = sigma*epsilon*T

    return(hr)

def qr(epsilon, Te, Ta):
    
    h = hr(epsilon, Te, Ta)
    
    return(h*(Te-Ta))



# =============================================================================
# Convecção natural. A variável U é requerida apenas por coerência interna.
# =============================================================================

# =============================================================================
# Coeficiente de transferência de calor por convecção natural em uma placa
# plana vertical. Input em m e K, output em W/(m^2.K).
# =============================================================================

def hc_n_ppv(U, Lc, Te, Ta):
    
    ''' Dimensão característica: Altura '''
    
    Tf = (Te + Ta)/2
    
    PrL = ar.Pr(Tf)
    
    RaL = Ra(Lc, Te, Ta)
    
    if RaL <= (10**9):
        
        a = 0.670*(RaL**0.25)
        
        b = (1 + ((0.492/PrL)**(9/16)))**(4/9)
        
        c = a/b
        
        NuL = 0.68 + c
        
    else:
        
        a = 0.387*(RaL**(1/6))
        
        b = (1 + ((0.492/PrL)**(9/16)))**(8/27)
        
        c = a/b
        
        NuL = (0.825 + c)**2
        
    hc = (NuL*(ar.lamed(Tf)))/Lc
    
    return (hc)

def qc_n_ppv(U, Lc, Te, Ta):
    
    hc = hc_n_ppv(U, Lc, Te, Ta)
    
    return (hc*(Te-Ta))

# =============================================================================
# Coeficiente de transferência de calor por convecção natural em um cilindro
# vertical. Input em m e K, output em W/(m^2.K).
# =============================================================================
    
def hc_n_cv(U, Lc, Te, Ta):
    
    ''' Dimensão característica: Altura '''
       
    return (hc_n_ppv(U, Lc, Te, Ta))

def qc_n_cv(U, Lc, Te, Ta):
    
    hc = hc_n_cv(U, Lc, Te, Ta)
    
    return (hc*(Te - Ta))

# =============================================================================
# Coeficiente de transferência de calor por convecção natural em um cilindro
# horizontal. Input em m e K, output em W/(m^2.K).
# =============================================================================

def hc_n_ch(U, Lc, Te, Ta):
    
    ''' Dimensão característica: Diâmetro Externo '''
    
    Tf = (Te + Ta)/2
    
    PrL = ar.Pr(Tf)
    
    RaL = Ra(Lc, Te, Ta)
    
    a = 0.60
    
    b = 0.387*(RaL**(1.0/6.0))
    
    c = (1 + ((0.559/PrL)**(9.0/16.0)))**(8.0/27.0)
    
    NuL = (a + (b/c))**2
    
    hc = (NuL*(ar.lamed(Tf)))/Lc
    
    return (hc)

def qc_n_ch(U, Lc, Te, Ta):
    
    hc = hc_n_ch(U, Lc, Te, Ta)
    
    return (hc*(Te-Ta))



# =============================================================================
# Convecção forçada.
# =============================================================================

# =============================================================================
# Coeficiente de transferência de calor por convecção forçada em um cilindro.
# Input em m/s, m e K, output em W/(m^2.K).
# =============================================================================

def hc_f_c(U, Lc, Te, Ta):
    
    Tf = (Te + Ta)/2
    
    PrL = ar.Pr(Tf)
    
    ReL = Re(U, Lc, Te, Ta)
    
    a = (ReL/282000)**(5.0/8.0)
    
    b = (1 + a)**(4.0/5.0)
    
    c = 0.62*(ReL**0.5)*(PrL**(1.0/3.0))*b
    
    d = (0.4/PrL)**(2.0/3.0)
    
    e = (1 + d)**(0.25)
    
    f = c/e
    
    NuL = 0.3 + f
    
    hc = (NuL*(ar.lamed(Tf)))/Lc
    
    return (hc)

def qc_f_c(U, Lc, Te, Ta):
    
    hc = hc_f_c(U, Lc, Te, Ta)
    
    return (hc*(Te-Ta))



# =============================================================================
# Convecção mista.
# =============================================================================

# =============================================================================
# Coeficiente de transferência de calor por convecção mista em um cilindro
# horizontal. Input em m/s, m e K, output em W/(m^2.K).
# =============================================================================

def hc_m_ch(U, Lc, Te, Ta):
    
    hf = hc_f_c(U, Lc, Te, Ta)
    
    hn = hc_n_ch(U, Lc, Te, Ta)
    
    hm = ((hf**4)+(hn**4))**(0.25)
    
    return(hm)

def qc_m_ch(U, Lc, Te, Ta):
    
    hm = hc_m_ch(U, Lc, Te, Ta)
    
    return (hm*(Te-Ta))

# =============================================================================
# Coeficiente de transferência de calor por convecção mista em um cilindro
# horizontal. Input em m/s, m e K, output em W/(m^2.K).
# =============================================================================

def hc_m_cv(U, D, H, Te, Ta):
    
    hf = hc_f_c(U, D, Te, Ta)
    
    hn = hc_n_cv(U, H, Te, Ta)
    
    hm = ((hf**4)+(hn**4))**(0.25)
    
    return(hm)

def qc_m_cv(U, D, H, Te, Ta):
    
    hm = hc_m_cv(U, D, H, Te, Ta)
    
    return (hm*(Te-Ta))