# -*- coding: utf-8 -*-

import numpy as np

# =============================================================================
# Propriedades de termoisolantes pré-moldados de silicato de cálcio. 
# Referência [006] - ABNT NBR 10662 - Isolantes Térmicos Pré-Moldados de 
# Silicato de Cálcio.
# =============================================================================

def lin(t, t2, t1, l2, l1):
    
    a = l1
    
    b = l2 - l1
    
    x = (t - t1)/(t2 - t1)
    
    return(a + b*x)
    
# =============================================================================
# Termoisolantes Tipo I. Massa específica de 240 kg/m^3. Temperatura máxima de 
# utilização de 650 °C. Intervalo de temperatura média para o cálculo da 
# condutividade térmica de 38 °C a 371 °C.
# =============================================================================

# =============================================================================
# Condutividade térmica. Input em K, output em W/(m.K).
# =============================================================================
    
def lamedI(T):
    
    t = T - 273.15
    
    if t < 38:
        
        return(0.059)
    
    if t <= 93:
        
        lamed = lin(t, 93, 38, 0.065, 0.059)
        
        return(lamed)
        
    if t <= 149:
        
        lamed = lin(t, 149, 93, 0.072, 0.065)
        
        return(lamed)
        
    if t <= 260:
        
        lamed = lin(t, 260, 149, 0.087, 0.072)
        
        return(lamed)
        
    if t <= 371: 
        
        lamed = lin(t, 371, 260, 0.102, 0.087)
        
        return(lamed)
        
    if t > 371:
        
        lamed = lin(t, 371, 260, 0.102, 0.087)
        
        return(lamed)
    
# =============================================================================
# Intervalo de temperatura média para o cálculo da condutividade térmica.
# Output em K.
# =============================================================================

def tmdI():
    
    tmd = np.array([38, 371])
    
    tmd = tmd + 273.15
    
    return(tmd)
    
# =============================================================================
# Temperatura máxima de utilização. Output em K.
# =============================================================================

def tmxI():
    
    tmx = 650 + 273.15
    
    return(tmx)



# =============================================================================
# Termoisolantes Tipo II. Massa específica de 240 kg/m^3. Temperatura máxima de 
# utilização de 815 °C. Intervalo de temperatura média para o cálculo da 
# condutividade térmica de 38 °C a 371 °C.
# =============================================================================

# =============================================================================
# Condutividade térmica. Input em K, output em W/(m.K).
# =============================================================================
    
def lamedII(T):
    
    t = T - 273.15
    
    if t < 38:
        
        return(0.059)
    
    if t <= 93:
        
        lamed = lin(t, 93, 38, 0.065, 0.059)
        
        return(lamed)
        
    if t <= 149:
        
        lamed = lin(t, 149, 93, 0.072, 0.065)
        
        return(lamed)
        
    if t <= 260:
        
        lamed = lin(t, 260, 149, 0.087, 0.072)
        
        return(lamed)
        
    if t <= 371: 
        
        lamed = lin(t, 371, 260, 0.102, 0.087)
        
        return(lamed)
        
    if t > 371:
        
        lamed = lin(t, 371, 260, 0.102, 0.087)
        
        return(lamed)
    
# =============================================================================
# Intervalo de temperatura média para o cálculo da condutividade térmica.
# Output em K.
# =============================================================================

def tmdII():
    
    tmd = np.array([38, 371])
    
    tmd = tmd + 273.15
    
    return(tmd)
    
# =============================================================================
# Temperatura máxima de utilização. Output em K.
# =============================================================================

def tmxII():
    
    tmx = 815 + 273.15
    
    return(tmx)



# =============================================================================
# Termoisolantes Tipo III. Massa específica de 352 kg/m^3. Temperatura máxima
# de utilização de 927 °C. Intervalo de temperatura média para o cálculo da 
# condutividade térmica de 38 °C a 538 °C.
# =============================================================================

# =============================================================================
# Condutividade térmica. Input em K, output em W/(m.K).
# =============================================================================
    
def lamedIII(T):
    
    t = T - 273.15
    
    if t < 38:
        
        return(0.072)
        
    if t <= 93:
        
        lamed = lin(t, 93, 38, 0.078, 0.072)
        
        return(lamed)
        
    if t <= 149:
        
        lamed = lin(t, 149, 93, 0.084, 0.078)
        
        return(lamed)
        
    if t <= 260:
        
        lamed = lin(t, 260, 149, 0.092, 0.084)
        
        return(lamed)
        
    if t <= 371:
        
        lamed = lin(t, 371, 260, 0.101, 0.092)
        
        return(lamed)
        
    if t <= 427:
        
        lamed = lin(t, 427, 371, 0.105, 0.101)
        
        return(lamed)
        
    if t <= 482:
        
        lamed = lin(t, 482, 427, 0.108, 0.105)
        
        return(lamed)
        
    if t <= 538:
        
        lamed = lin(t, 538, 482, 0.111, 0.108)
        
        return(lamed)
        
    if t > 538: 
        
        lamed = lin(t, 538, 482, 0.111, 0.108)
        
        return(lamed)
    
# =============================================================================
# Intervalo de temperatura média para o cálculo da condutividade térmica.
# Output em K.
# =============================================================================

def tmdIII():
    
    tmd = np.array([38, 538])
    
    tmd = tmd + 273.15
    
    return(tmd)
    
# =============================================================================
# Temperatura máxima de utilização. Output em K.
# =============================================================================

def tmxIII():
    
    tmx = 927 + 273.15
    
    return(tmx)