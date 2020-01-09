# -*- coding: utf-8 -*-

import numpy as np

# =============================================================================
# Propriedades de painéis termoisolantes à base de lã de cerâmica. 
# Referência [004] - ABNT NBR 9909 - Isolantes Térmicos de Lã de Cerâmica - 
# Painéis.
# =============================================================================

# =============================================================================
# A condutividade térmica é igual para as três classes, portanto só será
# definida uma vez.
# =============================================================================

def lin(t, t2, t1, l2, l1):
    
    a = l1
    
    b = l2 - l1
    
    x = (t - t1)/(t2 - t1)
    
    return(a + b*x)

# =============================================================================
# Condutividade térmica para painéis de massa específica de 64 kg/m^3. Input em 
# K, output em W/(m.K).
# =============================================================================

def lamed64(T):
    
    t = T - 273.15
    
    if t < 24:
        
        return(0.072)
    
    if t <= 93:
        
        lamed = lin(t, 93, 24, 0.078, 0.072)
        
        return(lamed)
        
    if t <= 204:
        
        lamed = lin(t, 204, 93, 0.098, 0.078)

        return(lamed)
        
    if t <= 427:
        
        lamed = lin(t, 427, 204, 0.163, 0.098)
        
        return(lamed)
        
    if t <= 649:
        
        lamed = lin(t, 649, 427, 0.275, 0.163)
        
        return(lamed)
        
    if t <= 871:
        
        lamed = lin(t, 871, 649, 0.446, 0.275)
        
        return(lamed)
        
    if t <= 1093:
        
        lamed = lin(t, 1093, 871, 0.676, 0.446)
        
        return(lamed)
        
    if t <= 1371:
        
        lamed = lin(t, 1371, 1093, 1.061, 0.676)
        
        return(lamed)
        
    if t > 1371:
        
        lamed = lin(t, 1371, 1093, 1.061, 0.676)
        
        return(lamed)

# =============================================================================
# Condutividade térmica para painéis de massa específica de 96 kg/m^3. Input em 
# K, output em W/(m.K).
# =============================================================================

def lamed96(T):
    
    t = T - 273.15
    
    if t < 24:
        
        return(0.060)
    
    if t <= 93:
        
        lamed = lin(t, 93, 24, 0.068, 0.060)
        
        return(lamed)
        
    if t <= 204:
        
        lamed = lin(t, 204, 93, 0.086, 0.068)

        return(lamed)
        
    if t <= 427:
        
        lamed = lin(t, 427, 204, 0.149, 0.086)
        
        return(lamed)
        
    if t <= 649:
        
        lamed = lin(t, 649, 427, 0.211, 0.149)
        
        return(lamed)
        
    if t <= 871:
        
        lamed = lin(t, 871, 649, 0.332, 0.211)
        
        return(lamed)
        
    if t <= 1093:
        
        lamed = lin(t, 1093, 871, 0.493, 0.332)
        
        return(lamed)
        
    if t <= 1371:
        
        lamed = lin(t, 1371, 1093, 0.759, 0.493)
        
        return(lamed)
        
    if t > 1371:
        
        lamed = lin(t, 1371, 1093, 0.759, 0.493)
        
        return(lamed)

# =============================================================================
# Condutividade térmica para painéis de massa específica de 128 kg/m^3. Input 
# em K, output em W/(m.K).
# =============================================================================

def lamed128(T):
    
    t = T - 273.15
    
    if t < 24:
        
        return(0.059)
    
    if t <= 93:
        
        lamed = lin(t, 93, 24, 0.066, 0.059)
        
        return(lamed)
        
    if t <= 204:
        
        lamed = lin(t, 204, 93, 0.081, 0.066)

        return(lamed)
        
    if t <= 427:
        
        lamed = lin(t, 427, 204, 0.146, 0.081)
        
        return(lamed)
        
    if t <= 649:
        
        lamed = lin(t, 649, 427, 0.179, 0.146)
        
        return(lamed)
        
    if t <= 871:
        
        lamed = lin(t, 871, 649, 0.273, 0.179)
        
        return(lamed)
        
    if t <= 1093:
        
        lamed = lin(t, 1093, 871, 0.395, 0.273)
        
        return(lamed)
        
    if t <= 1371:
        
        lamed = lin(t, 1371, 1093, 0.589, 0.395)
        
        return(lamed)
        
    if t > 1371:
        
        lamed = lin(t, 1371, 1093, 0.589, 0.395)
        
        return(lamed)

# =============================================================================
# Condutividade térmica para painéis de massa específica de 160 kg/m^3. Input
# em K, output em W/(m.K).
# =============================================================================

def lamed160(T):
    
    t = T - 273.15
    
    if t < 24:
        
        return(0.058)
    
    if t <= 93:
        
        lamed = lin(t, 93, 24, 0.065, 0.058)
        
        return(lamed)
        
    if t <= 204:
        
        lamed = lin(t, 204, 93, 0.079, 0.065)

        return(lamed)
        
    if t <= 427:
        
        lamed = lin(t, 427, 204, 0.141, 0.079)
        
        return(lamed)
        
    if t <= 649:
        
        lamed = lin(t, 649, 427, 0.161, 0.141)
        
        return(lamed)
        
    if t <= 871:
        
        lamed = lin(t, 871, 649, 0.243, 0.161)
        
        return(lamed)
        
    if t <= 1093:
        
        lamed = lin(t, 1093, 871, 0.347, 0.243)
        
        return(lamed)
        
    if t <= 1371:
        
        lamed = lin(t, 1371, 1093, 0.508, 0.347)
        
        return(lamed)
        
    if t > 1371:
        
        lamed = lin(t, 1371, 1093, 0.508, 0.347)
        
        return(lamed)

# =============================================================================
# Condutividade térmica para painéis de massa específica de 192 kg/m^3. Input
# em K, output em W/(m.K).
# =============================================================================

def lamed192(T):
    
    t = T - 273.15
    
    if t < 24:
        
        return(0.053)
    
    if t <= 93:
        
        lamed = lin(t, 93, 24, 0.062, 0.053)
        
        return(lamed)
        
    if t <= 204:
        
        lamed = lin(t, 204, 93, 0.076, 0.062)

        return(lamed)
        
    if t <= 427:
        
        lamed = lin(t, 427, 204, 0.138, 0.076)
        
        return(lamed)
        
    if t <= 649:
        
        lamed = lin(t, 649, 427, 0.147, 0.138)
        
        return(lamed)
        
    if t <= 871:
        
        lamed = lin(t, 871, 649, 0.215, 0.147)
        
        return(lamed)
        
    if t <= 1093:
        
        lamed = lin(t, 1093, 871, 0.300, 0.215)
        
        return(lamed)
        
    if t <= 1371:
        
        lamed = lin(t, 1371, 1093, 0.429, 0.300)
        
        return(lamed)
        
    if t > 1371:
        
        lamed = lin(t, 1371, 1093, 0.429, 0.300)
        
        return(lamed)

# =============================================================================
# O intervalo de temperatura média para o cálculo da condutividade térmica é
# igual para as três classes, portanto só será definido uma vez. Output em K.
# =============================================================================

def tmd():
    
    tmd = np.array[24, 1371]
    
    tmd = tmd + 273.15
    
    return(tmd)

# =============================================================================
# Temperatura máxima de utilização para painéis classe A. Output em K.
# =============================================================================

def tmxA():
    
    tmx = 1260
    
    tmx = tmx + 273.15
    
    return(tmx)
    
# =============================================================================
# Temperatura máxima de utilização para painéis classe B. Output em K.
# =============================================================================

def tmxB():
    
    tmx = 1426
    
    tmx = tmx + 273.15
    
    return(tmx)

# =============================================================================
# Temperatura máxima de utilização para painéis classe C. Output em K.
# =============================================================================

def tmxC():
    
    tmx = 815
    
    tmx = tmx + 273.15
    
    return(tmx)