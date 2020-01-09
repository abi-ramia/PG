# -*- coding: utf-8 -*-

import numpy as np

# =============================================================================
# Propriedades de flocos termoisolantes à base de lã de vidro. Referências
# [009] - ABNT NBR 11360 - Isolantes Térmicos de Lã de Vidro - Flocos e [008] -
# ABNT NBR 11358 - Painéis Termoisolantes à Base de Lã de Vidro.
# =============================================================================

def lin(t, t2, t1, l2, l1):
    
    a = l1
    
    b = l2 - l1
    
    x = (t - t1)/(t2 - t1)
    
    return(a + b*x)

# =============================================================================
# Flocos.
# =============================================================================

# =============================================================================
# Condutividade térmica. Input em K, output em W/(m.K).
# =============================================================================

def lamed(T):
    
    t = T - 273.15
    
    if t < -5:
        
        return(0.035)
        
    if t <= 35:
        
        lamed = lin(t, 35, -5, 0.044, 0.035)
        
        return(lamed)
        
    if t <= 60:
        
        lamed = lin(t, 60, 35, 0.048, 0.044)
        
        return(lamed)
        
    if t <= 85:
        
        lamed = lin(t, 85, 60, 0.056, 0.048)
        
        return(lamed)
        
    if t <= 110:
        
        lamed = lin(t, 110, 85, 0.063, 0.056)
        
        return(lamed)
        
    if t <= 135:
        
        lamed = lin(t, 135, 110, 0.072, 0.063)
        
        return(lamed)
        
    if t > 135:
        
        lamed = lin(t, 135, 110, 0.072, 0.063)
        
        return(lamed)
        
# =============================================================================
# Intervalo de temperatura média para o cálculo da condutividade térmica.
# Output em K.
# =============================================================================

def tmd():
    
    tmd = np.array([-5, 135])
    
    tmd = tmd + 273.15
    
    return(tmd)

# =============================================================================
# A temperatura máxima de utilização não é definida pela norma.
# =============================================================================

def tmx():
    
    return(np.nan)