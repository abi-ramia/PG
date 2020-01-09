# -*- coding: utf-8 -*-

import numpy as np

# =============================================================================
# Propriedades de cimento termoisolantes à base de silicato de Cálcio para
# Rejuntamento. Referência [014] - ABNT NBR 11777 - Cimento Isolante à Base de
# Silicato de Cálcio para Rejuntamento.
# =============================================================================

def lin(t, t2, t1, l2, l1):
    
    a = l1
    
    b = l2 - l1
    
    x = (t - t1)/(t2 - t1)
    
    return(a + b*x)
    
# =============================================================================
# Cimento à base de silicato de cálcio para rejuntamento.
# =============================================================================

# =============================================================================
# Condutividade térmica. Input em K, output em W/(m.K).
# =============================================================================

def lamed(T):
    
    t = T - 273.15
    
    if t < 95:
        
        return(0.0317)
        
    if t <= 150:
        
        lamed = lin(t, 150, 95, 0.331, 0.0317)
        
        return(lamed)
        
    if t <= 200:
        
        lamed = lin(t, 200, 150, 0.346, 0.331)
        
        return(lamed)
        
    if t > 200:
        
        lamed = lin(t, 200, 150, 0.346, 0.331)
        
        return(lamed)

# =============================================================================
# Intervalo de temperatura média para o cálculo da condutividade térmica.
# Output em K.
# =============================================================================

def tmd():
    
    tmd = np.array([95, 200])
    
    tmd = tmd + 273.15
    
    return(tmd)

# =============================================================================
# A temperatura máxima de utilização não é definida pela norma.
# =============================================================================

def tmx():
    
    return(np.nan)