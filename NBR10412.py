# -*- coding: utf-8 -*-

import numpy as np

# =============================================================================
# Propriedades de feltros de lamelas termoisolantes à base de lã de vidro. 
# Referência [005] - ABNT NBR 10412 - Isolantes Térmicos de Lã de Vidro - 
# Feltros de Lamelas.
# =============================================================================

def lin(t, t2, t1, l2, l1):
    
    a = l1
    
    b = l2 - l1
    
    x = (t - t1)/(t2 - t1)
    
    return(a + b*x)

# =============================================================================
# Feltros de massa específica de 60 kg/m^3.
# =============================================================================

# =============================================================================
# Condutividade térmica. Input em K, output em W/(m.K).
# =============================================================================

def lamed60(T):
    
    t = T - 273.15
    
    if t < 24:
        
        return(0.040)
    
    if t <= 160:
        
        lamed = lin(t, 160, 24, 0.080, 0.040)
        
        return(lamed)
        
    if t <= 185:
        
        lamed = lin(t, 185, 160, 0.087, 0.080)

        return(lamed)
        
    if t > 185:
        
        lamed = lin(t, 185, 160, 0.087, 0.080)
        
        return(lamed)

# =============================================================================
# Intervalo de temperatura média para o cálculo da condutividade térmica.
# Output em K.
# =============================================================================

def tmd60():
    
    tmd = np.array[24, 185]
    
    tmd = tmd + 273.15
    
    return(tmd)

# =============================================================================
# A temperatura máxima de utilização não é definida pela norma.
# =============================================================================

def tmx60():
    
    return(np.nan)



# =============================================================================
# Feltros de massa específica de 100 kg/m^3.
# =============================================================================

# =============================================================================
# Condutividade térmica. Input em K, output em W/(m.K).
# =============================================================================

def lamed100(T):
    
    t = T - 273.15
    
    if t < 24:
        
        return(0.040)
    
    if t <= 160:
        
        lamed = lin(t, 160, 24, 0.080, 0.040)
        
        return(lamed)
        
    if t <= 185:
        
        lamed = lin(t, 185, 160, 0.087, 0.080)

        return(lamed)
        
    if t > 185:
        
        lamed = lin(t, 185, 160, 0.087, 0.080)
        
        return(lamed)

# =============================================================================
# Intervalo de temperatura média para o cálculo da condutividade térmica.
# Output em K.
# =============================================================================

def tmd100():
    
    tmd = np.array[24, 185]
    
    tmd = tmd + 273.15
    
    return(tmd)

# =============================================================================
# A temperatura máxima de utilização não é definida pela norma.
# =============================================================================

def tmx100():
    
    return(np.nan)