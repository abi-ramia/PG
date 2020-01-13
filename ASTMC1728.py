# -*- coding: utf-8 -*-

import numpy as np

# =============================================================================
# Propriedades de tubos flexíveis de aerogel. 
# Referência [018] - NIA Insulation Materials Spec Chart - ASTM C1728.
# =============================================================================

def lin(t, t2, t1, l2, l1):
    
    a = l1
    
    b = l2 - l1
    
    x = (t - t1)/(t2 - t1)
    
    return(a + b*x)
    
# =============================================================================
# Termoisolantes Tipo III, Grau 2. Temperatura máxima de utilização de 650 °C.
# Intervalo de temperatura média para o cálculo da condutividade térmica de
# 24 °C a 316 °C.
# =============================================================================

# =============================================================================
# Condutividade térmica. Input em K, output em W/(m.K).
# =============================================================================
    
def lamed(T):
    
    t = T - 273.15
    
    t = 32 + 1.8*t
    
    if t < 75:
        
        return(0.16*0.1441314338)
    
    if t <= 200:
        
        lamed = lin(t, 200, 75, 0.18, 0.16)
        
        return(lamed*0.1441314338)
        
    if t <= 400:
        
        lamed = lin(t, 400, 200, 0.24, 0.18)
        
        return(lamed*0.1441314338)
        
    if t <= 600:
        
        lamed = lin(t, 600, 400, 0.33, 0.24)
        
        return(lamed*0.1441314338)
        
    if t > 600:
        
        lamed = lin(t, 600, 400, 0.33, 0.24)
        
        return(lamed*0.1441314338)
    
    return(True)
    
# =============================================================================
# Intervalo de temperatura média para o cálculo da condutividade térmica.
# Output em K.
# =============================================================================

def tmdI():
    
    tmd = np.array([24, 316])
    
    tmd = tmd + 273.15
    
    return(tmd)
    
# =============================================================================
# Temperatura máxima de utilização. Output em K.
# =============================================================================

def tmxI():
    
    tmx = 650 + 273.15
    
    return(tmx)