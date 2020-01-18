# -*- coding: utf-8 -*-

import numpy as np

# =============================================================================
# Funções de analogia de circuitos elétricos. Referências [016] - N-550 -
# Projeto de Isolamento Térmico a Alta Temperatura, Anexo E; e [017] -
# Incropera - Fundamentals of Heat & Mass Transfer - 8th Edition.
# =============================================================================



# =============================================================================
# Resistências térmicas em cilindros. As funções são desenvolvidas levando-se 
# em conta que os cálculos serão por unidade de comprimento (para tubulações).
# =============================================================================

# =============================================================================
# Resistências térmicas de condução por unidade de comprimento em um cilindro. 
# Input em m e W/(m.K), output em m.K/W.
# =============================================================================

def rt_cond_cili(Di, De, lmd):
    
    a = np.log(De/Di)
    
    b = 2.0*np.pi*1*lmd
    
    return (a/b)

# =============================================================================
# Resistência térmica de convecção por unidade de área em um cilindro. Input em
# m e W/(m^2.K), output em m.K/W.
# =============================================================================

def rt_conv_cili(D, h):
    
    return (1.0/(h*np.pi*D))

# =============================================================================
# Resistência térmica de radiação por unidade de área em um cilindro. Input em
# m e W/(m^2.K), output em m.K/W.
# =============================================================================

def rt_crad_cili(D, h):
    
    return (1.0/(h*np.pi*D))