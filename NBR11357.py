# -*- coding: utf-8 -*-

import numpy as np

# =============================================================================
# Propriedades de tubos termoisolantes à base de lã de vidro. Referência
# [007] - ABNT NBR 11357 - Tubos Termoisolantes à Base de Lã de Vidro.
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
    
    if t < 5:
        
        return(0.032)
    
    if t <= 35:
        
        lamed = lin(t, 35, 5, 0.037, 0.032)
        
        return(lamed)
    
    if t <= 60:
        
        lamed = lin(t, 60, 35, 0.041, 0.037)
        
        return(lamed)
        
    if t <= 85:
        
        lamed = lin(t, 85, 60, 0.045, 0.041)
        
        return(lamed)
        
    if t <= 110:
        
        lamed = lin(t, 110, 85, 0.049, 0.045)
        
        return(lamed)
        
    if t <= 135:
        
        lamed = lin(t, 135, 110, 0.054, 0.049)
        
        return(lamed)
        
    if t <= 160:
        
        lamed = lin(t, 160, 135, 0.061, 0.054)
        
        return(lamed)
        
    if t <= 185:
        
        lamed = lin(t, 185, 160, 0.066, 0.061)
        
        return(lamed)
        
    if t <= 210:
        
        lamed = lin(t, 210, 185, 0.075, 0.066)
        
        return(lamed)
        
    if t <= 235:
        
        lamed = lin(t, 235, 185, 0.086, 0.075)
        
        return(lamed)
        
    if t > 235:
        
        lamed = lin(t, 235, 185, 0.086, 0.075)
        
        return(lamed)

# =============================================================================
# Intervalo de temperatura média para o cálculo da condutividade térmica.
# Output em K.
# =============================================================================

def tmd():
    
    tmd = np.array([5, 235])
    
    tmd = tmd + 273.15
    
    return(tmd)
    
# =============================================================================
# A temperatura máxima de utilização não é definida pela norma.
# =============================================================================

def tmx():
    
    return(np.nan)