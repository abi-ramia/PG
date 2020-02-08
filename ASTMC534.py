# -*- coding: utf-8 -*-

import numpy as np

# =============================================================================
# Propriedades de tubos de espuma elastomérica. 
# Referência [018] - NIA Insulation Materials Spec Chart - ASTM C534.
# =============================================================================

def lin(t, t2, t1, l2, l1):
    
    a = l1
    
    b = l2 - l1
    
    x = (t - t1)/(t2 - t1)
    
    return(a + b*x)
    
# =============================================================================
# Termoisolantes Grau III. Temperatura máxima de utilização de 105 °C.
# Intervalo de temperatura média para o cálculo da condutividade térmica de
# -17 °C a 93 °C.
# =============================================================================

# =============================================================================
# Condutividade térmica. Input em K, output em W/(m.K).
# =============================================================================
    
def lamed(T):
    
    t = T - 273.15
    
    t = 32 + 1.8*t
    
    if t < 0:
        
        return(0.26*0.1441314338)
    
    if t <= 75:
        
        lamed = lin(t, 75, 0, 0.28, 0.26)
        
        return(lamed*0.1441314338)
    
    if t <= 200:
        
        lamed = lin(t, 200, 75, 0.31, 0.28)
        
        return(lamed*0.1441314338)
    
    if t > 200:
        
        lamed = lin(t, 200, 75, 0.31, 0.28)
        
        return(lamed*0.1441314338)
    
    return(True)
    
# =============================================================================
# Intervalo de temperatura média para o cálculo da condutividade térmica.
# Output em K.
# =============================================================================

def tmdI():
    
    tmd = np.array([-73, 93])
    
    tmd = tmd + 273.15
    
    return(tmd)
    
# =============================================================================
# Temperatura máxima de utilização. Output em K.
# =============================================================================

def tmxI():
    
    tmx = 105 + 273.15
    
    return(tmx)