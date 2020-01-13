# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from scipy import optimize
import ctc
import rt
import ASTMC1728
import NBR9688
import NBR9909
import NBR10412
import NBR10662
import NBR11357
import NBR11358
import NBR11360
import NBR11361
import NBR11363
import NBR11364
import NBR11722
import NBR11777
import NBR13047

# =============================================================================
# Análise de Tubulações.
# =============================================================================

# =============================================================================
# Tubulações em convecção forçada.
# =============================================================================

#Reduz o número de variáveis a um, a temperatura da face externa do isolante.
def generate_err_tubes_for(Ti, Ta, Di, U, E, eps, flmd):
    
    def err_tubes_for(Te):
        
        return (Te - Ti + (ctc.qc_f_c(U, Di + 2*E, Te, Ta) + ctc.qr(eps, Te, Ta))*rt.rt_cond_cili(Di + 2*E, Di, flmd((Te+Ti)/2)))
    
    return (err_tubes_for)
        
        
#Função principal.
def iso_tubes_for(Ti, Ta, Di, U, eps):
    
    #Lista de espessuras consideradas.
    LE = [(12.7*x)/1000 for x in range(1, 31)]
    
    #Lista de espessuras consideradas arredondadas.
    LE_Disp =  [(int(12.7*x)) for x in range(1,31)]
    
    #Lista de espessuras consideradas em polegadas.
    LE_Disp_Imp = [0.5*x for x in range (1, 31)]
    
    #Lista de funções de condutividade térmica.
    LLMD = [ASTMC1728.lamed,
            NBR10662.lamedI, NBR10662.lamedII, NBR10662.lamedIII,
            NBR10412.lamed60, NBR10412.lamed100,
            NBR11357.lamed,
            NBR11363.lamed,
            NBR11722.lamed]
    
    #Lista de nomes de isolantes.
    Lnm = ['Aerogel',
           'Silicato de Cálcio Tipo I',
           'Silicato de Cálcio Tipo II',
           'Silicato de Cálcio Tipo III',
           'Feltro de Lamelas de Lã de Vidro D60',
           'Feltro de Lamelas de Lã de Vidro D100',
           'Tubo de Lã de Vidro',
           'Tubo de Lã de Rocha',
           'Feltro de Lamelas de Lã de Rocha']
    
    #Número de espessuras consideradas.
    ne = len(LE)
    
    #Número de isolantes considerados.
    ni = len(LLMD)
    
    #Lista de nomes para o DataFrame.
    LNM = []
    
    for d in Lnm:
        
        LNM += [d]*ne
    
    #Lista de soluções para a temperatura na face externa.
    LTe = []
    
    #Lista de soluções para o fluxo de calor.
    Lq = []
    
    for flmd in LLMD:
        for E in LE:
            err = generate_err_tubes_for(Ti, Ta, Di, U, E, eps, flmd)
            root = optimize.brentq(err, Ta, Ti)
            LTe = LTe + [root]
            Lq = Lq + [ctc.qc_f_c(U, Di + 2*E, root, Ta) + ctc.qr(eps, root, Ta)]
    
    #Adição do caso sem isolante e organização para o DataFrame.
    LNM = ['Sem Isolante'] + LNM
    LE_Disp_DF = [0] + LE_Disp*ni
    LE_Disp_Imp_DF = [0] + LE_Disp_Imp*ni
    LTe = [Ti] + LTe
    Lq = [ctc.qc_f_c(U, Di, Ti, Ta) + ctc.qr(eps, Ti, Ta)] + Lq
    
    #Lista de soluções para a temperatura na face externa em °C.
    Lte = [(Te - 273.15) for Te in LTe]
    
    #Organização dos dados em um DataFrame.
    Disp = pd.DataFrame({'Material' : LNM,
                         'Espessura [mm]' : LE_Disp_DF,
                         'Espessura [pol]' : LE_Disp_Imp_DF,
                         'Temperatura [K]' : LTe,
                         'Temperatura [°C]' : Lte,
                         'Fluxo de Calor [W/m^2]': Lq})
    
    return Disp



# =============================================================================
# Tubulações horizontais em convecção natural.
# =============================================================================

#Reduz o número de variáveis a um, a temperatura da face externa do isolante.
def generate_err_tubes_nat_h(Ti, Ta, Di, E, eps, flmd):
    
    def err_tubes_nat_h(Te):
        
        return (Te - Ti + (ctc.qc_n_ch(Di + 2*E, Te, Ta) + ctc.qr(eps, Te, Ta))*rt.rt_cond_cili(Di + 2*E, Di, flmd((Te+Ti)/2)))
    
    return (err_tubes_nat_h)
        
        
#Função principal.
def iso_tubes_nat_h(Ti, Ta, Di, eps):
    
    #Lista de espessuras consideradas.
    LE = [(12.7*x)/1000 for x in range(1, 31)]
    
    #Lista de espessuras consideradas arredondadas.
    LE_Disp = [(int(12.7*x)) for x in range(1,31)]
    
    #Lista de espessuras consideradas em polegadas.
    LE_Disp_Imp = [0.5*x for x in range (1, 31)]
    
    #Lista de funções de condutividade térmica.
    LLMD = [ASTMC1728.lamed,
            NBR10662.lamedI, NBR10662.lamedII, NBR10662.lamedIII,
            NBR10412.lamed60, NBR10412.lamed100,
            NBR11357.lamed,
            NBR11363.lamed,
            NBR11722.lamed]
    
    #Lista de nomes de isolantes.
    Lnm = ['Aerogel',
           'Silicato de Cálcio Tipo I',
           'Silicato de Cálcio Tipo II',
           'Silicato de Cálcio Tipo III',
           'Feltro de Lamelas de Lã de Vidro D60',
           'Feltro de Lamelas de Lã de Vidro D100',
           'Tubo de Lã de Vidro',
           'Tubo de Lã de Rocha',
           'Feltro de Lamelas de Lã de Rocha']
    
    #Número de espessuras consideradas.
    ne = len(LE)
    
    #Número de isolantes considerados.
    ni = len(LLMD)
    
    #Lista de nomes para o DataFrame.
    LNM = []
    
    for d in Lnm:
        
        LNM += [d]*ne
    
    #Lista de soluções para a temperatura na face externa.
    LTe = []
    
    #Lista de soluções para o fluxo de calor.
    Lq = []
    
    for flmd in LLMD:
        for E in LE:
            err = generate_err_tubes_nat_h(Ti, Ta, Di, E, eps, flmd)
            root = optimize.brentq(err, Ta, Ti)
            LTe = LTe + [root]
            Lq = Lq + [ctc.qc_n_ch(Di + 2*E, root, Ta) + ctc.qr(eps, root, Ta)]
    
    #Adição do caso sem isolante e organização para o DataFrame.
    LNM = ['Sem Isolante'] + LNM
    LE_Disp_DF = [0] + LE_Disp*ni
    LE_Disp_Imp_DF = [0] + LE_Disp_Imp*ni
    LTe = [Ti] + LTe
    Lq = [ctc.qc_n_ch(Di, Ti, Ta) + ctc.qr(eps, Ti, Ta)] + Lq
    
    #Lista de soluções para a temperatura na face externa em °C.
    Lte = [(Te - 273.15) for Te in LTe]
    
    #Organização dos dados em um DataFrame.
    Disp = pd.DataFrame({'Material' : LNM,
                         'Espessura [mm]' : LE_Disp_DF,
                         'Espessura [pol]' : LE_Disp_Imp_DF,
                         'Temperatura [K]' : LTe,
                         'Temperatura [°C]' : Lte,
                         'Fluxo de Calor [W/m^2]': Lq})
    
    return Disp



# =============================================================================
# Tubulações verticais em convecção natural.
# =============================================================================

#Reduz o número de variáveis a um, a temperatura da face externa do isolante.
def generate_err_tubes_nat_v(Ti, Ta, H, Di, E, eps, flmd):
    
    def err_tubes_nat_v(Te):
        
        return (Te - Ti + (ctc.qc_n_cv(H, Te, Ta, Di + 2*E) + ctc.qr(eps, Te, Ta))*rt.rt_cond_cili(Di + 2*E, Di, flmd((Te+Ti)/2)))
    
    return (err_tubes_nat_v)



#Função principal.
def iso_tubes_nat_v(Ti, Ta, H, Di, eps):
    
    #Lista de espessuras consideradas.
    LE = [(12.7*x)/1000 for x in range(1,31)]
    
    #Lista de espessuras consideradas arredondadas.
    LE_Disp = [(int(12.7*x)) for x in range(1,31)]
    
    #Lista de espessuras consideradas em polegadas.
    LE_Disp_Imp = [0.5*x for x in range(1,31)]
    
    #Lista de funções de condutividade térmica.
    LLMD = [ASTMC1728.lamed,
            NBR10662.lamedI, NBR10662.lamedII, NBR10662.lamedIII,
            NBR10412.lamed60, NBR10412.lamed100,
            NBR11357.lamed,
            NBR11363.lamed,
            NBR11722.lamed]
    
    #Lista de nomes de isolantes.
    Lnm = ['Aerogel',
           'Silicato de Cálcio Tipo I',
           'Silicato de Cálcio Tipo II',
           'Silicato de Cálcio Tipo III',
           'Feltro de Lamelas de Lã de Vidro D60',
           'Feltro de Lamelas de Lã de Vidro D100',
           'Tubo de Lã de Vidro',
           'Tubo de Lã de Rocha',
           'Feltro de Lamelas de Lã de Rocha']
    
    #Número de espessuras consideradas.
    ne = len(LE)
    
    #Número de isolantes considerados.
    ni = len(LLMD)
    
    #Lista de nomes para o DataFrame.
    LNM = []
    
    for d in Lnm:
        
        LNM += [d]*ne
    
    #Lista de soluções para a temperatura na face externa.
    LTe = []
    
    #Lista de soluções para o fluxo de calor.
    Lq = []
    
    for flmd in LLMD:
        for E in LE:
            err = generate_err_tubes_nat_v(Ti, Ta, H, Di, E, eps, flmd)
            root = optimize.brentq(err, Ta, Ti)
            LTe = LTe + [root]
            Lq = Lq + [ctc.qc_n_cv(H, root, Ta, Di + 2*E) + ctc.qr(eps, root, Ta)]
    
    #Adição do caso sem isolante e organização para o DataFrame.
    LNM = ['Sem Isolante'] + LNM
    LE_Disp_DF = [0] + LE_Disp*ni
    LE_Disp_Imp_DF = [0] + LE_Disp_Imp*ni
    LTe = [Ti] + LTe
    Lq = [ctc.qc_n_cv(H, Ti, Ta, Di) + ctc.qr(eps, Ti, Ta)] + Lq
    
    #Lista de soluções para a temperatura na face externa em °C.
    Lte = [(Te - 273.15) for Te in LTe]
    
    #Organização dos dados em um DataFrame.
    Disp = pd.DataFrame({'Material' : LNM,
                         'Espessura [mm]' : LE_Disp_DF,
                         'Espessura [pol]' : LE_Disp_Imp_DF,
                         'Temperatura [K]' : LTe,
                         'Temperatura [°C]' : Lte,
                         'Fluxo de Calor [W/m^2]': Lq})
    
    return Disp



# =============================================================================
# Custo de energia perdida.
# =============================================================================

#Custo de energia perdida trazido para valor atual em $/(ano.m^2)
def CE_VA(q, N, F, eta, n, i, delta):
    
    CE = (3600*q*N*F)/eta
    
    j = ((1 + i)/(1 + delta)) - 1
    
    f = (((1 + j)**n) - 1)/(j*((1 + j)**n))
    
    CEVA = f*CE
    
    return (CEVA)

# =============================================================================
# Tubulações quaisquer.
# =============================================================================

def iso_tubes(Ti, Ta, Di, H, U, eps, N, F, eta, n, i, delta):
    
    if (U != 0):
        
        Disp = iso_tubes_for(Ti, Ta, Di, U, eps)
        
    if ((U == 0) and (H != 0)):
        
        Disp = iso_tubes_nat_v(Ti, Ta, H, Di, eps)
        
    if ((U == 0) and (H == 0)):
        
        Disp = iso_tubes_nat_h(Ti, Ta, Di, eps)
        
    Lq = Disp['Fluxo de Calor [W/m^2]']
    
    Lq = np.array(Lq)
    
    Lq = CE_VA(Lq, N, F, eta, n, i, delta)
    
    Disp['Custo de Energia Perdida [$/(ano.m^2)]'] = Lq
        
    return Disp