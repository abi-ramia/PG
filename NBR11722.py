# -*- coding: utf-8 -*-

import numpy as np

# =============================================================================
# Propriedades de feltros termoisolantes à base de lã de rocha. Referência
# [013] - ABNT NBR 11722 - Feltros Termoisolantes à Base de Lã de Rocha.
# =============================================================================

def lin(t, t2, t1, l2, l1):
    
    a = l1
    
    b = l2 - l1
    
    x = (t - t1)/(t2 - t1)
    
    return(a + b*x)
    
# =============================================================================
# Feltors de massa específica de 32 kg/m^3.
# =============================================================================

# =============================================================================
# Condutividade térmica. Input em K, output em W/(m.K).
# =============================================================================

def lamed(T):
    
    t = T - 273.15
    
    if t < 50:
        
        return(0.050)
        
    if t <= 100:
        
        lamed = lin(t, 100, 50, 0.067, 0.050)
        
        return(lamed)
        
    if t > 100:
        
        lamed = lin(t, 100, 50, 0.067, 0.050)
        
        return(lamed)

# =============================================================================
# Intervalo de temperatura média para o cálculo da condutividade térmica.
# Output em K.
# =============================================================================

def tmd():
    
    tmd = np.array([50, 100])
    
    tmd = tmd + 273.15
    
    return(tmd)

# =============================================================================
# A temperatura máxima de utilização não é definida pela norma.
# =============================================================================

def tmx():
    
    return(np.nan)