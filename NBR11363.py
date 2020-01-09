# -*- coding: utf-8 -*-

import numpy as np

# =============================================================================
# Propriedades de tubos termoisolantes à base de lã de rocha. Referência
# [011] - ABNT NBR 11363 - Tubos Termoisolantes à Base de Lã de Rocha.
# =============================================================================

def lin(t, t2, t1, l2, l1):
    
    a = l1
    
    b = l2 - l1
    
    x = (t - t1)/(t2 - t1)
    
    return(a + b*x)

# =============================================================================
# Tubos.
# =============================================================================

# =============================================================================
# Condutividade térmica. Input em K, output em W/(m.K).
# =============================================================================

def lamed(T):
    
    t = T - 273.15
    
    if t < 50:
        
        return(0.044)
        
    if t <= 100:
        
        lamed = lin(t, 100, 50, 0.052, 0.044)
        
        return(lamed)
        
    if t <= 150:
        
        lamed = lin(t, 150, 100, 0.061, 0.052)
        
        return(lamed)
        
    if t <= 200:
        
        lamed = lin(t, 200, 150, 0.072, 0.061)
        
        return(lamed)
        
    if t <= 250:
        
        lamed = lin(t, 250, 200, 0.085, 0.072)
        
        return(lamed)
        
    if t <= 300:
        
        lamed = lin(t, 300, 250, 0.100, 0.085)
        
        return(lamed)
        
    if t <= 350:
        
        lamed = lin(t, 350, 300, 0.121, 0.100)
        
        return(lamed)
        
    if t > 350:
        
        lamed = lin(t, 350, 300, 0.121, 0.100)
        
        return(lamed)
        
# =============================================================================
# Intervalo de temperatura média para o cálculo da condutividade térmica.
# Output em K.
# =============================================================================

def tmd():
    
    tmd = np.array([50, 350])
    
    tmd = tmd + 273.15
    
    return(tmd)

# =============================================================================
# A temperatura máxima de utilização não é definida pela norma.
# =============================================================================

def tmx():
    
    return(np.nan)