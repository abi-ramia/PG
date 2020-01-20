# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from scipy import optimize
import ctc
import rt
import ASTMC1728
import NBR10412
import NBR10662
import NBR11357
import NBR11363
import NBR11722



# =============================================================================
# Tubulações horizontais em convecção combinada.
# =============================================================================

#Reduz o número de variáveis a um, a temperatura na interface com o ar.
def generate_err_tubes_for_h(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd):
    
    Di = de
    
    def err_tubes_for(Te):
        
        qdpp = ctc.qc_m_ch(U, Di + 2*E, Te, Ta) + ctc.qr(eps, Te, Ta)
        
        qdp = qdpp*(np.pi*(Di + 2*E))
        
        Tdi = Ti - qdp*rt.rt_conv_cili(di, h_fld)
        
        Tde = Tdi - qdp*rt.rt_cond_cili(di, de, lmd_tube)
        
        TDi = Tde
        
        #TODO Caso se adicione manta de proteção, mude esta linha.
        TDe = Te
        
        lmd_iso = flmd((TDe + TDi)/2)
        
        Te_calc = Ti - qdp*(rt.rt_conv_cili(di, h_fld) + rt.rt_cond_cili(di, de, lmd_tube) + rt.rt_cond_cili(Di, Di + 2*E, lmd_iso))
        
        return (Te - Te_calc)
    
    return (err_tubes_for)

#Caso sem isolante.
def generate_err_tubes_for_h_si(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps):
    
    Di = de
    
    def err_tubes_for(Te):
        
        qdpp = ctc.qc_m_ch(U, Di, Te, Ta) + ctc.qr(eps, Te, Ta)
        
        qdp = qdpp*(np.pi*(Di))
        
        Tdi = Ti - qdp*rt.rt_conv_cili(di, h_fld)
        
        Tde = Tdi - qdp*rt.rt_cond_cili(di, de, lmd_tube)
        
        Te_calc = Tde
        
        return (Te - Te_calc)
    
    return (err_tubes_for)

# =============================================================================
# Tubulações verticais em convecção combinada.
# =============================================================================

#Reduz o número de variáveis a um, a temperatura na interface com o ar.
def generate_err_tubes_for_v(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd):
    
    Di = de
    
    def err_tubes_for(Te):
        
        qdpp = ctc.qc_m_cv(U, Di + 2*E, H, Te, Ta) + ctc.qr(eps, Te, Ta)
        
        qdp = qdpp*(np.pi*(Di + 2*E))
        
        Tdi = Ti - qdp*rt.rt_conv_cili(di, h_fld)
        
        Tde = Tdi - qdp*rt.rt_cond_cili(di, de, lmd_tube)
        
        TDi = Tde
        
        #TODO Caso se adicione manta de proteção, mude esta linha.
        TDe = Te
        
        lmd_iso = flmd((TDe + TDi)/2)
        
        Te_calc = Ti - qdp*(rt.rt_conv_cili(di, h_fld) + rt.rt_cond_cili(di, de, lmd_tube) + rt.rt_cond_cili(Di, Di + 2*E, lmd_iso))
        
        return (Te - Te_calc)
    
    return (err_tubes_for)

#Caso sem isolante.
def generate_err_tubes_for_v_si(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps):
    
    Di = de
    
    def err_tubes_for(Te):
        
        qdpp = ctc.qc_m_cv(U, Di, H, Te, Ta) + ctc.qr(eps, Te, Ta)
        
        qdp = qdpp*(np.pi*(Di))
        
        Tdi = Ti - qdp*rt.rt_conv_cili(di, h_fld)
        
        Tde = Tdi - qdp*rt.rt_cond_cili(di, de, lmd_tube)
        
        Te_calc = Tde
        
        return (Te - Te_calc)
    
    return (err_tubes_for)

# =============================================================================
# Tubulações horizontais em convecção natural.
# =============================================================================

#Reduz o número de variáveis a um, a temperatura na interface com o ar.
def generate_err_tubes_nat_h(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd):
    
    Di = de
    
    def err_tubes_nat_h(Te):
        
        qdpp = ctc.qc_n_ch(U, Di + 2*E, Te, Ta) + ctc.qr(eps, Te, Ta)
        
        qdp = qdpp*(np.pi*(Di + 2*E))
        
        Tdi = Ti - qdp*rt.rt_conv_cili(di, h_fld)
        
        Tde = Tdi - qdp*rt.rt_cond_cili(di, de, lmd_tube)
        
        TDi = Tde
        
        #TODO Caso se adicione manta de proteção, mude esta linha.
        TDe = Te
        
        lmd_iso = flmd((TDe + TDi)/2)
        
        Te_calc = Ti - qdp*(rt.rt_conv_cili(di, h_fld) + rt.rt_cond_cili(di, de, lmd_tube) + rt.rt_cond_cili(Di, Di + 2*E, lmd_iso))
        
        return (Te - Te_calc)
    
    return (err_tubes_nat_h)

#Caso em isolante.
def generate_err_tubes_nat_h_si(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps):
    
    Di = de
    
    def err_tubes_nat_h(Te):
        
        qdpp = ctc.qc_n_ch(U, Di, Te, Ta) + ctc.qr(eps, Te, Ta)
        
        qdp = qdpp*(np.pi*(Di))
        
        Tdi = Ti - qdp*rt.rt_conv_cili(di, h_fld)
        
        Tde = Tdi - qdp*rt.rt_cond_cili(di, de, lmd_tube)
        
        Te_calc = Tde
        
        return (Te - Te_calc)
    
    return (err_tubes_nat_h)

# =============================================================================
# Tubulações verticais em convecção natural.
# =============================================================================

#Reduz o número de variáveis a um, a temperatura na interface com o ar.
def generate_err_tubes_nat_v(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd):
    
    Di = de
    
    def err_tubes_nat_v(Te):
        
        qdpp = ctc.qc_n_cv(U, H, Te, Ta) + ctc.qr(eps, Te, Ta)
        
        qdp = qdpp*(np.pi*(Di + 2*E))
        
        Tdi = Ti - qdp*rt.rt_conv_cili(di, h_fld)
        
        Tde = Tdi - qdp*rt.rt_cond_cili(di, de, lmd_tube)
        
        TDi = Tde
        
        #TODO Caso se adicione manta de proteção, mude esta linha.
        TDe = Te
        
        lmd_iso = flmd((TDe + TDi)/2)
        
        Te_calc = Ti - qdp*(rt.rt_conv_cili(di, h_fld) + rt.rt_cond_cili(di, de, lmd_tube) + rt.rt_cond_cili(Di, Di + 2*E, lmd_iso))
        
        return (Te - Te_calc)
    
    return (err_tubes_nat_v)

#Caso sem isoalnte
def generate_err_tubes_nat_v_si(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps):
    
    Di = de
    
    def err_tubes_nat_v(Te):
        
        qdpp = ctc.qc_n_cv(U, H, Te, Ta) + ctc.qr(eps, Te, Ta)
        
        qdp = qdpp*(np.pi*(Di))
        
        Tdi = Ti - qdp*rt.rt_conv_cili(di, h_fld)
        
        Tde = Tdi - qdp*rt.rt_cond_cili(di, de, lmd_tube)
        
        Te_calc = Tde
        
        return (Te - Te_calc)
    
    return (err_tubes_nat_v)

# =============================================================================
# Função principal.
# =============================================================================

def iso_tubes(di, de, Ti, Ta, h_fld, lmd_tube, U, H, z, eps, fase_change, m, c_or_h):
    
    Di = de
    
    #Lista de espessuras consideradas.
    LE = [(12.7*x)/1000 for x in range(1, 31)]
    
    #Lista de espessuras consideradas arredondadas.
    LE_Disp =  [(int(12.7*x)) for x in range(1,31)]
    
    #Lista de espessuras consideradas em polegadas.
    LE_Disp_Imp = [0.5*x for x in range (1, 31)]
    
    #Lista de funções de condutividade térmica.
    LLMD = [ASTMC1728.lamed,
            NBR10662.lamedII, NBR10662.lamedIII,
            NBR10412.lamed60, NBR10412.lamed100,
            NBR11357.lamed,
            NBR11363.lamed,
            NBR11722.lamed]
    
    #Lista de nomes de isolantes.
    Lnm = ['Aerogel',
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
    
    #Lista de espessuras para o DataFrame.
    LE_Disp_DF = [0] + LE_Disp*ni
    LE_Disp_Imp_DF = [0] + LE_Disp_Imp*ni
    
    #Lista de nomes para o DataFrame.
    LNM = ['Sem Isolante']
    
    for d in Lnm:
        
        LNM += [d]*ne
    
    #Lista de soluções para a temperatura na face externa.
    LTe = []
    
    #Lista de soluções para a condutividade térmica.
    Lslmd = [np.NaN]
    
    #Lista de soluções para o fluxo de calor na face externa.
    Lq = []
    
    #Lista de diâmetros externos.
    LDe = [de]
    
    #Lista de taxas de formação de condensado.
    LFC = []
    
    #Lista de resistências térmicas.
    LR = []
    
    #Lista de variações de temperatura do fluido.
    LVT = []
    
    #Lista de temperaturas na saída.
    LTS = []
    
    if ((U != 0) and (H == 0)):
        
        errf0 = generate_err_tubes_for_h_si
        
        err0 = errf0(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps)
        
        root0 = optimize.brentq(err0, Ta, Ti)
        
        LTe = [root0]
        
        Lq = [(ctc.qc_m_ch(U, Di, root0, Ta) + ctc.qr(eps, root0, Ta))*(np.pi*(Di))]
        
        LR = [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+((((rt.rt_conv_cili(de,ctc.hc_m_ch(U,de,root0,Ta)))**(-1))+((rt.rt_crad_cili(de,ctc.hr(eps,root0,Ta)))**(-1)))**(-1)))/z]
        
        errf = generate_err_tubes_for_h
        
        for flmd in LLMD:
            for E in LE:
                err = errf(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd)
                root = optimize.brentq(err, Ta, Ti)
                LTe = LTe + [root]
                Lslmd = Lslmd + [flmd((root + Ti)/2)]
                Lq = Lq + [(ctc.qc_m_ch(U, Di + 2*E, root, Ta) + ctc.qr(eps, root, Ta))*(np.pi*(Di + 2*E))]
                LDe = LDe + [Di + 2*E]
                LR = LR + [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+rt.rt_cond_cili(Di,Di + 2*E,flmd((root+Ti)/2))+((((rt.rt_conv_cili(Di + 2*E,ctc.hc_m_ch(U,Di + 2*E,root,Ta)))**(-1))+((rt.rt_crad_cili(Di + 2*E,ctc.hr(eps,root,Ta)))**(-1)))**(-1)))/z]
        
    if ((U != 0) and (H != 0)):
        
        errf0 = generate_err_tubes_for_v_si
        
        err0 = errf0(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps)
        
        root0 = optimize.brentq(err0, Ta, Ti)
        
        LTe = [root0]
        
        Lq = [(ctc.qc_m_cv(U, Di, H, root0, Ta) + ctc.qr(eps, root0, Ta))*(np.pi*(Di))]
        
        LR = [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+((((rt.rt_conv_cili(Di,ctc.hc_m_cv(U,Di,H,root0,Ta)))**(-1))+((rt.rt_crad_cili(Di,ctc.hr(eps,root0,Ta)))**(-1)))**(-1)))/z]
        
        errf = generate_err_tubes_for_v
        
        for flmd in LLMD:
            for E in LE:
                err = errf(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd)
                root = optimize.brentq(err, Ta, Ti)
                LTe = LTe + [root]
                Lslmd = Lslmd + [flmd((root + Ti)/2)]
                Lq = Lq + [(ctc.qc_m_cv(U, Di + 2*E, H, root, Ta) + ctc.qr(eps, root, Ta))*(np.pi*(Di + 2*E))]
                LDe = LDe + [Di + 2*E]
                LR = LR + [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+rt.rt_cond_cili(Di,Di + 2*E,flmd((root+Ti)/2))+((((rt.rt_conv_cili(Di + 2*E,ctc.hc_m_cv(U,Di + 2*E,H,root,Ta)))**(-1))+((rt.rt_crad_cili(Di + 2*E,ctc.hr(eps,root,Ta)))**(-1)))**(-1)))/z]
        
    if ((U == 0) and (H == 0)):
        
        errf0 = generate_err_tubes_nat_h_si
        
        err0 = errf0(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps)
        
        root0 = optimize.brentq(err0, Ta, Ti)
        
        LTe = [root0]
        
        Lq = [(ctc.qc_n_ch(U, Di, root0, Ta) + ctc.qr(eps, root0, Ta))*(np.pi*(Di))]
        
        LR = [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+((((rt.rt_conv_cili(de,ctc.hc_n_ch(U,de,root0,Ta)))**(-1))+((rt.rt_crad_cili(de,ctc.hr(eps,root0,Ta)))**(-1)))**(-1)))/z]
        
        errf = generate_err_tubes_nat_h
        
        for flmd in LLMD:
            for E in LE:
                err = errf(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd)
                root = optimize.brentq(err, Ta, Ti)
                LTe = LTe + [root]
                Lslmd = Lslmd + [flmd((root + Ti)/2)]
                Lq = Lq + [(ctc.qc_n_ch(U, Di + 2*E, root, Ta) + ctc.qr(eps, root, Ta))*(np.pi*(Di + 2*E))]
                LDe = LDe + [Di + 2*E]
                LR = LR + [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+rt.rt_cond_cili(Di,Di + 2*E,flmd((root+Ti)/2))+((((rt.rt_conv_cili(Di + 2*E,ctc.hc_n_ch(U,Di + 2*E,root,Ta)))**(-1))+((rt.rt_crad_cili(Di + 2*E,ctc.hr(eps,root,Ta)))**(-1)))**(-1)))/z]
        
    if ((U == 0) and (H != 0)):
        
        errf0 = generate_err_tubes_nat_v_si
        
        err0 = errf0(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps)
        
        root0 = optimize.brentq(err0, Ta, Ti)
        
        LTe = [root0]
        
        Lq = [(ctc.qc_n_cv(U, H, root0, Ta) + ctc.qr(eps, root0, Ta))*(np.pi*(Di))]
        
        LR = [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+((((rt.rt_conv_cili(de,ctc.hc_n_cv(U,H,root0,Ta)))**(-1))+((rt.rt_crad_cili(de,ctc.hr(eps,root0,Ta)))**(-1)))**(-1)))/z]
        
        errf = generate_err_tubes_nat_v
        
        for flmd in LLMD:
            for E in LE:
                err = errf(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd)
                root = optimize.brentq(err, Ta, Ti)
                LTe = LTe + [root]
                Lslmd = Lslmd + [flmd((root + Ti)/2)]
                Lq = Lq + [(ctc.qc_n_cv(U, H, root, Ta) + ctc.qr(eps, root, Ta))*(np.pi*(Di + 2*E))]
                LDe = LDe + [Di + 2*E]
                LR = LR + [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+rt.rt_cond_cili(Di,Di + 2*E,flmd((root+Ti)/2))+((((rt.rt_conv_cili(Di + 2*E,ctc.hc_n_cv(U,H,root,Ta)))**(-1))+((rt.rt_crad_cili(Di + 2*E,ctc.hr(eps,root,Ta)))**(-1)))**(-1)))/z]            
            
    #Diâmetros externos em milímetros.
    LDe_Disp = [1000*D for D in LDe]
    
    #Lista de soluções para a temperatura na face externa em °C.
    Lte = [(Te - 273.15) for Te in LTe]
    
    #Nome da coluna de queda de temperatura ou de formação de condensado.
    if (fase_change):
        noc = 'Formação de Condensado [kg/s]'
        LFC = list(map(lambda x: (x*z)/c_or_h, Lq))
        Lnoc = LFC
    else:
        noc = 'Variação de Temperatura do Fluido [°C]'
        LVT = list(map(lambda x: Ti - (Ta - (Ta - Ti)*np.exp(-1/(m*c_or_h*x))),LR))
        Lnoc = LVT
    
    #Organização dos dados em um DataFrame.
    Disp = pd.DataFrame({'Material' : LNM,
                         'Espessura [mm]' : LE_Disp_DF,
                         'Espessura [pol]' : LE_Disp_Imp_DF,
                         'Diâmetro Externo [mm]' : LDe_Disp,
                         'Temperatura [K]' : LTe,
                         'Temperatura [°C]' : Lte,
                         'Condutividade Térmica [W/(m.k)]': Lslmd,
                         'Fluxo de Calor [W/m]': Lq,
                          noc : Lnoc})
    
    if (not(fase_change)):
        LTS = list(map(lambda x: Ti - 273.15 - x, LVT))
        Disp['Temperatura do Fluido na Saída [°C]'] = LTS
    
    return Disp



# =============================================================================
# Custos.
# =============================================================================

#Custo de energia perdida trazido para valor atual em $/(ano.m^2).
def CE_VA(q, N, F, eta, n, i, delta):
    
    CE = (3600*q*N*F)/eta
    
    j = ((1 + i)/(1 + delta)) - 1
    
    f = (((1 + j)**n) - 1)/(j*((1 + j)**n))
    
    CEVA = f*CE
    
    return (CEVA)

#Custo de investimanto no isolamento em $/(ano.m^2).
def CI_VA(nome, Di, De):
    
    return (True)

#Custo de manutenção do isolamento em $/(ano.m^2).
def CM_VA(CIVA, tm, n, i):
    
    CM = tm*CIVA
    
    f = (((1+i)**n)-1)/(i*((1+i)**n))
    
    CMVA = f*CM
    
    return (CMVA)



# =============================================================================
# Tubulações quaisquer.
# =============================================================================

#def iso_tubes(Ti, Ta, Di, H, U, eps, N, F, eta, n, i, delta):
#    
#    if (U != 0):
#        
#        Disp = iso_tubes_for(Ti, Ta, Di, U, eps)
#        
#    if ((U == 0) and (H != 0)):
#        
#        Disp = iso_tubes_nat_v(Ti, Ta, H, Di, eps)
#        
#    if ((U == 0) and (H == 0)):
#        
#        Disp = iso_tubes_nat_h(Ti, Ta, Di, eps)
#    
#    Disp['Custo de Energia Perdida [$/(ano.m^2)]'] = Disp['Fluxo de Calor (Face Externa) [W/m^2]'].apply(lambda x : CE_VA(x, N, F, eta, n, i, delta))
#    
#    return Disp