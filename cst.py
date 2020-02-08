# -*- coding: utf-8 -*-

import numpy as np

# =============================================================================
# Estimativa de Custos de Isolantes. Referência [023].
# =============================================================================



# =============================================================================
# Poliisocianurato. Referência [024].
# =============================================================================

#Dados amostrais.
def samples_C591():
    LDi = [0.5, 5, 11]
    LDe = [2.5, 10, 19]
    LD = list(zip(LDi, LDe))
    LA = [np.pi*(x[1]**2 - x[0]**2)/4 for x in LD]
    LAM = [x*25.4*25.4/1000000 for x in LA]
    LC = [10*1.890, 10*12.340, 10*34.450]
    return (LAM, LC)

#Interpolação dos ln.
def inter_C591(LL):
    LAM = LL[0]
    LC = LL[1]
    LnLAM = [np.log(x) for x in LAM]
    LnLC = [np.log(x) for x in LC]
    z = np.polyfit(LnLAM, LnLC, 1)
    p = np.poly1d(z)
    return (p)

#Estimação do custo.
def cst_C591(Di, E):
    A = np.pi*(((Di + 2*E)**2) - (Di**2))/4
    LnA = np.log(A)
    poly = inter_C591(samples_C591())
    LnC = poly(LnA)
    C = np.exp(LnC)
    return (C)



# =============================================================================
# Vidro Celular. Referência [024].
# =============================================================================

#Dados amostrais.
def samples_C552():
    LDi = [0.5, 5, 21]
    LDe = [1.5, 10, 29]
    LD = list(zip(LDi, LDe))
    LA = [np.pi*(x[1]**2 - x[0]**2)/4 for x in LD]
    LAM = [x*25.4*25.4/1000000 for x in LA]
    LC = [10*1.050, 10*9.720, 10*50.930]
    return (LAM, LC)

#Interpolação dos ln.
def inter_C552(LL):
    LAM = LL[0]
    LC = LL[1]
    LnLAM = [np.log(x) for x in LAM]
    LnLC = [np.log(x) for x in LC]
    z = np.polyfit(LnLAM, LnLC, 1)
    p = np.poly1d(z)
    return (p)

#Estimação do custo.
def cst_C552(Di, E):
    A = np.pi*(((Di + 2*E)**2) - (Di**2))/4
    LnA = np.log(A)
    poly = inter_C552(samples_C552())
    LnC = poly(LnA)
    C = np.exp(LnC)
    return (C)



# =============================================================================
# Silicato de Cálcio. Referência [024].
# =============================================================================

#Dados amostrais.
def samples_NBR10662():
    LDi = [3, 6, 24]
    LDe = [5, 10, 30]
    LD = list(zip(LDi, LDe))
    LA = [np.pi*(x[1]**2 - x[0]**2)/4 for x in LD]
    LAM = [x*25.4*25.4/1000000 for x in LA]
    LC = [10*3.000, 10*8.830, 10*41.840]
    return (LAM, LC)

#Interpolação dos ln.
def inter_NBR10662(LL):
    LAM = LL[0]
    LC = LL[1]
    LnLAM = [np.log(x) for x in LAM]
    LnLC = [np.log(x) for x in LC]
    z = np.polyfit(LnLAM, LnLC, 1)
    p = np.poly1d(z)
    return (p)

#Estimação do custo.
def cst_NBR10662(Di, E):
    A = np.pi*(((Di + 2*E)**2) - (Di**2))/4
    LnA = np.log(A)
    poly = inter_NBR10662(samples_NBR10662())
    LnC = poly(LnA)
    C = np.exp(LnC)
    return (C)



# =============================================================================
# Espuma Elastomérica.
# =============================================================================

#Dados amostrais.
def samples_C534():
    LDi = [0.875, 1.125, 2.125]
    LDe = [2.375, 2.625, 3.625]
    LD = list(zip(LDi, LDe))
    LA = [np.pi*(x[1]**2 - x[0]**2)/4 for x in LD]
    LAM = [x*25.4*25.4/1000000 for x in LA]
    LC = [10*1.749, 10*1.940, 10*2.650]
    return (LAM, LC)

#Interpolação dos ln.
def inter_C534(LL):
    LAM = LL[0]
    LC = LL[1]
    LnLAM = [np.log(x) for x in LAM]
    LnLC = [np.log(x) for x in LC]
    z = np.polyfit(LnLAM, LnLC, 1)
    p = np.poly1d(z)
    return (p)

#Estimação do custo.
def cst_C534(Di, E):
    A = np.pi*(((Di + 2*E)**2) - (Di**2))/4
    LnA = np.log(A)
    poly = inter_C534(samples_C534())
    LnC = poly(LnA)
    C = np.exp(LnC)
    return (C)



# =============================================================================
# Materias sem informações disponíveis.
# =============================================================================

#Custos não disponíveis.
def cst_null(Di, E):
    return(np.nan)