# -*- coding: utf-8 -*-

import numpy as np

# =============================================================================
# Propriedades de painéis termoisolantes à base de lã de rocha. Referência
# [012] - ABNT NBR 11364 - Painéis Termoisolantes à Base de Lã de Rocha.
# =============================================================================

def lin(t, t2, t1, l2, l1):
    
    a = l1
    
    b = l2 - l1
    
    x = (t - t1)/(t2 - t1)
    
    return(a + b*x)
    
# =============================================================================
# Painéis de massa específica de 32 kg/m^3.
# =============================================================================

# =============================================================================
# Condutividade térmica. Input em K, output em W/(m.K).
# =============================================================================

def lamed32(T):
    
    t = T - 273.15
    
    if t < 0:
        
        return(0.032)
        
    if t <= 24:
        
        lamed = lin(t, 24, 0, 0.042, 0.032)
        
        return(lamed)
        
    if t <= 50:
        
        lamed = lin(t, 50, 24, 0.049, 0.042)
        
        return(lamed)
        
    if t <= 75:
        
        lamed = lin(t, 75, 50, 0.057, 0.049)
        
        return(lamed)
        
    if t <= 100:
        
        lamed = lin(t, 100, 75, 0.067, 0.057)
        
        return(lamed)
        
    if t > 100:
        
        lamed = lin(t, 100, 75, 0.067, 0.057)
        
        return(lamed)
        
# =============================================================================
# Intervalo de temperatura média para o cálculo da condutividade térmica.
# Output em K.
# =============================================================================

def tmd32():
    
    tmd = np.array([0, 100])
    
    tmd = tmd + 273.15
    
    return(tmd)

# =============================================================================
# A temperatura máxima de utilização não é definida pela norma.
# =============================================================================

def tmx32():
    
    return(np.nan)



# =============================================================================
# Painéis de massa específica de 48 kg/m^3.
# =============================================================================

# =============================================================================
# Condutividade térmica. Input em K, output em W/(m.K).
# =============================================================================

def lamed48(T):
    
    t = T - 273.15
    
    if t < 0:
        
        return(0.032)
        
    if t <= 24:
        
        lamed = lin(t, 24, 0, 0.039, 0.032)
        
        return(lamed)
        
    if t <= 50:
        
        lamed = lin(t, 50, 24, 0.045, 0.039)
        
        return(lamed)
        
    if t <= 75:
        
        lamed = lin(t, 75, 50, 0.055, 0.045)
        
        return(lamed)
        
    if t <= 100:
        
        lamed = lin(t, 100, 75, 0.065, 0.055)
        
        return(lamed)
        
    if t > 100:
        
        lamed = lin(t, 100, 75, 0.065, 0.055)
        
        return(lamed)
    
# =============================================================================
# Intervalo de temperatura média para o cálculo da condutividade térmica.
# Output em K.
# =============================================================================

def tmd48():
    
    tmd = np.array([0, 100])
    
    tmd = tmd + 273.15
    
    return(tmd)
 
# =============================================================================
# A temperatura máxima de utilização é não definida pela norma.
# =============================================================================

def tmx48():
    
    return(np.nan)
    
    
    
# =============================================================================
# Painéis de massa específica de 64 kg/m^3.
# =============================================================================

# =============================================================================
# Condutividade térmica. Input em K, output em W/(m.K).
# =============================================================================

def lamed64(T):
    
    t = T - 273.15
    
    if t < 0:
        
        return(0.032)
        
    if t <= 24:
        
        lamed = lin(t, 24, 0, 0.035, 0.032)
        
        return(lamed)
        
    if t <= 50:
        
        lamed = lin(t, 50, 24, 0.045, 0.035)
        
        return(lamed)
        
    if t <= 75:
        
        lamed = lin(t, 75, 50, 0.051, 0.045)
        
        return(lamed)
        
    if t <= 100:
        
        lamed = lin(t, 100, 75, 0.057, 0.051)
        
        return(lamed)
        
    if t <= 150:
        
        lamed = lin(t, 150, 100, 0.067, 0.057)
        
        return(lamed)
        
    if t <= 200:
        
        lamed = lin(t, 200, 150, 0.085, 0.067)
        
        return(lamed)
        
    if t > 200:
        
        lamed = lin(t, 200, 150, 0.085, 0.067)
        
        return(lamed) 

# =============================================================================
# Intervalo de temperatura média para o cálculo da condutividade térmica.
# Output em K.
# =============================================================================

def tmd64():
    
    tmd = np.array([0, 200])
    
    tmd = tmd + 273.15
    
    return(tmd)

# =============================================================================
# A temperatura máxima de utilização não definida pela norma.
# =============================================================================
    
def tmx64():
    
    return(np.nan)
    
    
    
# =============================================================================
# Painéis de massa específica de 80 kg/m^3.
# =============================================================================

# =============================================================================
# Condutividade térmica. Input em K, output em W/(m.K).
# =============================================================================

def lamed80(T):
    
    t = T - 273.15
    
    if t < 50:
        
        return(0.043)
        
    if t <= 75:
        
        lamed = lin(t, 75, 50, 0.051, 0.043)
        
        return(lamed)
        
    if t <= 100:
        
        lamed = lin(t, 100, 75, 0.047, 0.051)
        
        return(lamed)
        
    if t <= 150:
        
        lamed = lin(t, 150, 100, 0.065, 0.047)
        
        return(lamed)
        
    if t <= 200:
        
        lamed = lin(t, 200, 150, 0.082, 0.065)
        
        return(lamed)
        
    if t <= 250:
        
        lamed = lin(t, 250, 200, 0.096, 0.082)
        
        return(lamed)
        
    if t <= 300:
        
        lamed = lin(t, 300, 250, 0.112, 0.096)
        
        return(lamed)
        
    if t > 300:
        
        lamed = lin(t, 300, 250, 0.112, 0.096)
        
        return(lamed)

# =============================================================================
# Intervalo de temperatura média para o cálculo da condutividade térmica.
# Output em K.
# =============================================================================

def tmd80():
    
    tmd = np.array([50, 300])
    
    tmd = tmd + 273.15
    
    return(tmd)
  
# =============================================================================
# A temperatura máxima de utilização não é definida pela norma.
# =============================================================================

def tmx80():
    
    return(np.nan)
    
    
    
# =============================================================================
# Painéis de massa específica de 96 kg/m^3.
# =============================================================================

# =============================================================================
# Condutividade térmica. Input em K, output em W/(m.K).
# =============================================================================

def lamed96(T):
    
    t = T - 273.15
    
    if t < 50:
        
        return(0.043)
        
    if t <= 75:
        
        lamed = lin(t, 75, 50, 0.050, 0.043)
        
        return(lamed)
        
    if t <= 100:
        
        lamed = lin(t, 100, 75, 0.045, 0.050)
        
        return(lamed)
        
    if t <= 150:
        
        lamed = lin(t, 150, 100, 0.060, 0.045)
        
        return(lamed)
        
    if t <= 200:
        
        lamed = lin(t, 200, 150, 0.075, 0.060)
        
        return(lamed)
        
    if t <= 250:
        
        lamed = lin(t, 250, 200, 0.088, 0.075)
        
        return(lamed)
        
    if t <= 300:
        
        lamed = lin(t, 300, 250, 0.095, 0.088)
        
        return(lamed)
        
    if t > 300:
        
        lamed = lin(t, 300, 250, 0.095, 0.088)
        
        return(lamed)

# =============================================================================
# Intervalo de temperatura média para o cálculo da condutividade térmica.
# Output em K.
# =============================================================================

def tmd96():
    
    tmd = np.array([50, 300])
    
    tmd = tmd + 273.15
    
    return(tmd)
    
# =============================================================================
# A temperatura máxima de utilização não é definida pela norma.
# =============================================================================

def tmx96():
    
    return(np.nan) 
    
    
    
# =============================================================================
# Painéis de massa específica de 128 kg/m^3.
# =============================================================================

# =============================================================================
# Condutividade térmica. Input em K, output em W/(m.K).
# =============================================================================

def lamed128(T):
    
    t = T - 273.15
    
    if t < 100:
        
        return(0.045)
        
    if t <= 150:
        
        lamed = lin(t, 150, 100, 0.058, 0.045)
        
        return(lamed)
        
    if t <= 200:
        
        lamed = lin(t, 200, 150, 0.070, 0.058)
        
        return(lamed)
        
    if t <= 250:
        
        lamed = lin(t, 250, 200, 0.085, 0.070)
        
        return(lamed)
        
    if t <= 300:
        
        lamed = lin(t, 300, 250, 0.103, 0.085)
        
        return(lamed)
        
    if t <= 350:
        
        lamed = lin(t, 350, 300, 0.126, 0.103)
        
        return(lamed)
        
    if t > 350:
        
        lamed = lin(t, 350, 300, 0.126, 0.103)
        
        return(lamed)

# =============================================================================
# Intervalo de temperatura média para o cálculo da condutividade térmica.
# Output em K.
# =============================================================================

def tmd128():
    
    tmd = np.array([100, 350])
    
    tmd = tmd + 273.15
    
    return(tmd)
 
# =============================================================================
# A temperatura máxima de utilização não é definida pela norma.
# =============================================================================

def tmx128():
    
    return(np.nan)  