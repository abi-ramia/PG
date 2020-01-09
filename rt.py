# -*- coding: utf-8 -*-

import numpy as np

# =============================================================================
# Funções de analogia de circuitos elétricos. Referências [016] - N-550 -
# Projeto de Isolamento Térmico a Alta Temperatura, Anexo E; e [017] -
# Incropera - Fundamentals of Heat & Mass Transfer - 8th Edition.
# =============================================================================



# =============================================================================
# Resistências térmicas em planos. As funções são desenvolvidas levando-se em
# conta que os cálculos serão por unidade de área (para tanques).
# =============================================================================

# =============================================================================
# Resistência térmica de condução por unidade de área em um plano. Input em m e
# W/(m.K), output em m^2.K/W.
# =============================================================================

def rt_cond_plan(L, lmd):
    
    return (L/lmd)

# =============================================================================
# Resistência térmica de convecção por unidade de área em um plano. Input em
# W/(m^2.K), output em m^2.K/W.
# =============================================================================

def rt_conv_plan(h):
    
    return (1.0/h)

# =============================================================================
# Resistência térmica de radiação por unidade de área em um plano. Input em
# W/(m^2.K), output em m^2.K/W.
# =============================================================================

def rt_crad_plan(h):
    
    return (1.0/h)



# =============================================================================
# Resistências térmicas em cilindros. As funções são desenvolvidas levando-se 
# em conta que os cálculos serão por unidade de comprimento (para tubulações).
# =============================================================================

# =============================================================================
# Resistências térmicas de condução por unidade de comprimento em um cilindro. 
# Input em m e W/(m.K), output em m^2.K/W.
# =============================================================================

def rt_cond_cili(De, Di, lmd):
    
    a = np.log(De/Di)
    
    b = 2.0*np.pi*lmd
    
    return (a/b)

# =============================================================================
# Resistência térmica de convecção por unidade de área em um cilindro. Input em
# W/(m^2.K), output em m^2.K/W.
# =============================================================================

def rt_conv_cili(h):
    
    return (1.0/h)

# =============================================================================
# Resistência térmica de radiação por unidade de área em um cilindro. Input em
# W/(m^2.K), output em m^2.K/W.
# =============================================================================

def rt_crad_cili(h):
    
    return (1.0/h)