# -*- coding: utf-8 -*-

import numpy as np

# =============================================================================
# Propriedades de painéis termoisolantes à base de lã de vidro. Referência
# [008] - ABNT NBR 11358 - Painéis Termoisolantes à Base de Lã de Vidro.
# =============================================================================

def lin(t, t2, t1, l2, l1):
    
    a = l1
    
    b = l2 - l1
    
    x = (t - t1)/(t2 - t1)
    
    return(a + b*x)

# =============================================================================
# Painéis de massa específica de 20 kg/m^3.
# =============================================================================

# =============================================================================
# Condutividade térmica. Input em K, output em W/(m.K).
# =============================================================================

def lamed20(T):
    
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

def tmd20():
    
    tmd = np.array([-5, 135])
    
    tmd = tmd + 273.15
    
    return(tmd)

# =============================================================================
# A temperatura máxima de utilização não é definida pela norma.
# =============================================================================

def tmx20():
    
    return(np.nan)



# =============================================================================
# Painéis de massa específica de 40 kg/m^3.
# =============================================================================

# =============================================================================
# Condutividade térmica. Input em K, output em W/(m.K).
# =============================================================================

def lamed40(T):
    
    t = T - 273.15
    
    if t < -5:
        
        return(0.031)
        
    if t <= 35:
        
        lamed = lin(t, 35, -5, 0.037, 0.031)
        
        return(lamed)
        
    if t <= 60:
        
        lamed = lin(t, 60, 35, 0.041, 0.037)
        
        return(lamed)
        
    if t <= 85:
        
        lamed = lin(t, 85, 60, 0.047, 0.041)
        
        return(lamed)
        
    if t <= 110:
        
        lamed = lin(t, 110, 85, 0.052, 0.047)
        
        return(lamed)
        
    if t <= 135:
        
        lamed = lin(t, 135, 110, 0.060, 0.052)
        
        return(lamed)
        
    if t <= 160:
        
        lamed = lin(t, 160, 135, 0.069, 0.060)
        
        return(lamed)
        
    if t <= 185:
        
        lamed = lin(t, 185, 160, 0.077, 0.069)
        
        return(lamed)
        
    if t > 185:
        
        lamed = lin(t, 185, 160, 0.077, 0.069)
        
        return(lamed)

# =============================================================================
# Intervalo de temperatura média para o cálculo da condutividade térmica.
# Output em K.
# =============================================================================

def tmd40():
    
    tmd = np.array([-5, 185])
    
    tmd = tmd + 273.15
    
    return(tmd)

# =============================================================================
# A temperatura máxima de utilização não é definida pela norma.
# =============================================================================

def tmx40():
    
    return(np.nan)



# =============================================================================
# Painéis de massa específica de 60 kg/m^3.
# =============================================================================

# =============================================================================
# Condutividade térmica. Input em K, output em W/(m.K).
# =============================================================================

def lamed60(T):
    
    t = T - 273.15
    
    if t < -5:
        
        return(0.032)
        
    if t <= 35:
        
        lamed = lin(t, 35, -5, 0.037, 0.032)
        
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
        
        lamed = lin(t, 235, 210, 0.086, 0.075)
        
        return(lamed)
        
    if t > 235:
        
        lamed = lin(t, 235, 210, 0.086, 0.075)
        
        return(lamed)

# =============================================================================
# Intervalo de temperatura média para o cálculo da condutividade térmica.
# Output em K.
# =============================================================================

def tmd60():
    
    tmd = np.array([-5, 235])
    
    tmd = tmd + 273.15
    
    return(tmd)
    
# =============================================================================
# A temperatura máxima de utilização não é definida pela norma.
# =============================================================================

def tmx60():
    
    return(np.nan)