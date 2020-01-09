# -*- coding: utf-8 -*-

import numpy as np

# =============================================================================
# Propriedades de mantas termoisolantes à base de lã de vidro. Referência
# [010] - ABNT NBR 11361 - Mantas Termoisolantes à Base de Lã de Vidro.
# =============================================================================

def lin(t, t2, t1, l2, l1):
    
    a = l1
    
    b = l2 - l1
    
    x = (t - t1)/(t2 - t1)
    
    return(a + b*x)

# =============================================================================
# Mantas.
# =============================================================================

# =============================================================================
# Condutividade térmica. Input em K, output em W/(m.K).
# =============================================================================

def lamed(T):
    
    t = T - 273.15
    
    if t < -45:
        
        return(0.031)
        
    if t <= 0:
        
        lamed = lin(t, 0, -45, 0.033, 0.031)
        
        return(lamed)
        
    if t <= 50:
        
        lamed = lin(t, 50, 0, 0.039, 0.033)
        
        return(lamed)
        
    if t <= 100:
        
        lamed = lin(t, 100, 50, 0.051, 0.039)
        
        return(lamed)
        
    if t <= 150:
        
        lamed = lin(t, 150, 100, 0.067, 0.051)
        
        return(lamed)
        
    if t <= 200:
        
        lamed = lin(t, 200, 150, 0.080, 0.067)
        
        return(lamed)
        
    if t > 200:
        
        lamed = lin(t, 200, 150, 0.080, 0.067)
        
        return(lamed)
        
# =============================================================================
# Intervalo de temperatura média para o cálculo da condutividade térmica.
# Output em K.
# =============================================================================

def tmd():
    
    tmd = np.array([-45, 200])
    
    tmd = tmd + 273.15
    
    return(tmd)

# =============================================================================
# A temperatura máxima de utilização não é definida pela norma.
# =============================================================================

def tmx():
    
    return(np.nan)