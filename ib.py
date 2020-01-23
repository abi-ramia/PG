# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from scipy import optimize
import ctc
import rt
import ASTMC552
import ASTMC591
import ASTMC610
import ASTMC1728
import NBR9688
import NBR10412
import NBR10662
import NBR11357
import NBR11363
import NBR11722

#TODO Adicionar resistência do revestimento de proteção.
#TODO Adicionar análise econômica.
#TODO Pivotear por menor custo.

# =============================================================================
# Tubulações horizontais em convecção combinada.
# =============================================================================

#Reduz o número de variáveis a um, a temperatura na interface com o ar.
def generate_err_tubes_for_h(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd, R_rev):
    
    Di = de
    
    def err_tubes_for(Te):
        
        qdpp = ctc.qc_m_ch(U, Di + 2*E, Te, Ta) + ctc.qr(eps, Te, Ta)
        
        qdp = qdpp*(np.pi*(Di + 2*E))
        
        Tdi = Ti - qdp*rt.rt_conv_cili(di, h_fld)
        
        Tde = Tdi - qdp*rt.rt_cond_cili(di, de, lmd_tube)
        
        TDi = Tde
        
        TDe = Te + R_rev*qdp
        
        lmd_iso = flmd((TDe + TDi)/2)
        
        Te_calc = Ti - qdp*(rt.rt_conv_cili(di, h_fld) + rt.rt_cond_cili(di, de, lmd_tube) + rt.rt_cond_cili(Di, Di + 2*E, lmd_iso) + R_rev)
        
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
def generate_err_tubes_for_v(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd, R_rev):
    
    Di = de
    
    def err_tubes_for(Te):
        
        qdpp = ctc.qc_m_cv(U, Di + 2*E, H, Te, Ta) + ctc.qr(eps, Te, Ta)
        
        qdp = qdpp*(np.pi*(Di + 2*E))
        
        Tdi = Ti - qdp*rt.rt_conv_cili(di, h_fld)
        
        Tde = Tdi - qdp*rt.rt_cond_cili(di, de, lmd_tube)
        
        TDi = Tde
        
        TDe = Te + R_rev*qdp
        
        lmd_iso = flmd((TDe + TDi)/2)
        
        Te_calc = Ti - qdp*(rt.rt_conv_cili(di, h_fld) + rt.rt_cond_cili(di, de, lmd_tube) + rt.rt_cond_cili(Di, Di + 2*E, lmd_iso) + R_rev)
        
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
def generate_err_tubes_nat_h(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd, R_rev):
    
    Di = de
    
    def err_tubes_nat_h(Te):
        
        qdpp = ctc.qc_n_ch(U, Di + 2*E, Te, Ta) + ctc.qr(eps, Te, Ta)
        
        qdp = qdpp*(np.pi*(Di + 2*E))
        
        Tdi = Ti - qdp*rt.rt_conv_cili(di, h_fld)
        
        Tde = Tdi - qdp*rt.rt_cond_cili(di, de, lmd_tube)
        
        TDi = Tde
        
        TDe = Te + R_rev*qdp
        
        lmd_iso = flmd((TDe + TDi)/2)
        
        Te_calc = Ti - qdp*(rt.rt_conv_cili(di, h_fld) + rt.rt_cond_cili(di, de, lmd_tube) + rt.rt_cond_cili(Di, Di + 2*E, lmd_iso) + R_rev)
        
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
def generate_err_tubes_nat_v(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd, R_rev):
    
    Di = de
    
    def err_tubes_nat_v(Te):
        
        qdpp = ctc.qc_n_cv(U, H, Te, Ta) + ctc.qr(eps, Te, Ta)
        
        qdp = qdpp*(np.pi*(Di + 2*E))
        
        Tdi = Ti - qdp*rt.rt_conv_cili(di, h_fld)
        
        Tde = Tdi - qdp*rt.rt_cond_cili(di, de, lmd_tube)
        
        TDi = Tde
        
        TDe = Te + R_rev*qdp
        
        lmd_iso = flmd((TDe + TDi)/2)
        
        Te_calc = Ti - qdp*(rt.rt_conv_cili(di, h_fld) + rt.rt_cond_cili(di, de, lmd_tube) + rt.rt_cond_cili(Di, Di + 2*E, lmd_iso) + R_rev)
        
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

def iso_tubes(di, de, Ti, Ta, h_fld, lmd_tube, U, H, z, eps, fase_change, m, c_or_h, ts_max, Dt_max, R_rev):
    
    Di = de
    
    #Lista de espessuras consideradas.
    LE1 = [(6.35*x)/1000 for x in range(1, 9)]
    LE2 = [(6.35*x)/1000 for x in range(1, 17)]
    LE3 = [(12.7*x)/1000 for x in range(1, 31)]
    
    #Lista de espessuras consideradas arredondadas.
    LE_Disp1 =  [(int(x*1000)) for x in LE1]
    LE_Disp2 =  [(int(x*1000)) for x in LE2]
    LE_Disp3 =  [(int(x*1000)) for x in LE3]
    
    #Lista de espessuras consideradas em polegadas.
    LE_Disp_Imp1 = [0.25*x for x in range (1, 9)]
    LE_Disp_Imp2 = [0.25*x for x in range (1, 17)]
    LE_Disp_Imp3 = [0.50*x for x in range (1, 31)]
    
    #Lista de funções de condutividade térmica.
    #Fibra cerâmica.
    LLMD1 = [NBR9688.lamed64,
             NBR9688.lamed96,
             NBR9688.lamed128,
             NBR9688.lamed160,
             NBR9688.lamed192]
    #Isolantes flexíveis.
    LLMD2 = [ASTMC1728.lamed,
             NBR10412.lamed60,
             NBR10412.lamed100,
             NBR11357.lamed,
             NBR11363.lamed,
             NBR11722.lamed]
    #Isolantes rígidos.
    LLMD3 = [NBR10662.lamedII,
             NBR10662.lamedIII,
             ASTMC610.lamed,
             ASTMC552.lamed,
             ASTMC591.lamed]
    
    #Lista de nomes de isolantes.
    #Fibra cerâmica.
    Lnm1 = ['Manta de Fibra \n Cerâmica D64',
            'Manta de Fibra \n Cerâmcia D96',
            'Manta de Fibra \n Cerâmica D128',
            'Manta de Fibra \n Cerâmica D160',
            'Manta de Fibra \n Cerâmica D192']
    #Isolantes flexívies.
    Lnm2 = ['Aerogel',
            'Feltro de Lamelas de \n Lã de Vidro D60',
            'Feltro de Lamelas de \n Lã de Vidro D100',
            'Tubo de Lã de Vidro',
            'Tubo de Lã de Rocha',
            'Feltro de Lamelas de \n Lã de Rocha']
    #Isolantes rígidos.
    Lnm3 = ['Silicato de Cálcio \n Tipo II',
            'Silicato de Cálcio \n Tipo III',
            'Perlita Expandida',
            'Vidro Celular',
            'Poliisocianurato']
    
    #Número de espessuras consideradas.
    ne1 = len(LE1)
    ne2 = len(LE2)
    ne3 = len(LE3)
    
    #Número de isolantes considerados.
    ni1 = len(LLMD1)
    ni2 = len(LLMD2)
    ni3 = len(LLMD3)    
    
    #Lista de espessuras para o DataFrame.
    LE_Disp_DF = [0] + LE_Disp1*ni1 + LE_Disp2*ni2 + LE_Disp3*ni3
    LE_Disp_Imp_DF = [0] + LE_Disp_Imp1*ni1 + LE_Disp_Imp2*ni2 + LE_Disp_Imp3*ni3
    
    #Lista de nomes para o DataFrame.
    LNM = ['Sem \n Isolante']
    
    for d in Lnm1: 
        LNM += [d]*ne1
    for d in Lnm2:
        LNM += [d]*ne2
    for d in Lnm3:
        LNM += [d]*ne3
    
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
        
        for flmd in LLMD1:
            for E in LE1:
                err = errf(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd, R_rev)
                root = optimize.brentq(err, Ta, Ti)
                LTe = LTe + [root]
                Lslmd = Lslmd + [flmd((root + Ti)/2)]
                Lq = Lq + [(ctc.qc_m_ch(U, Di + 2*E, root, Ta) + ctc.qr(eps, root, Ta))*(np.pi*(Di + 2*E))]
                LDe = LDe + [Di + 2*E]
                LR = LR + [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+rt.rt_cond_cili(Di,Di + 2*E,flmd((root+Ti)/2))+R_rev+((((rt.rt_conv_cili(Di + 2*E,ctc.hc_m_ch(U,Di + 2*E,root,Ta)))**(-1))+((rt.rt_crad_cili(Di + 2*E,ctc.hr(eps,root,Ta)))**(-1)))**(-1)))/z]
        for flmd in LLMD2:
            for E in LE2:
                err = errf(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd, R_rev)
                root = optimize.brentq(err, Ta, Ti)
                LTe = LTe + [root]
                Lslmd = Lslmd + [flmd((root + Ti)/2)]
                Lq = Lq + [(ctc.qc_m_ch(U, Di + 2*E, root, Ta) + ctc.qr(eps, root, Ta))*(np.pi*(Di + 2*E))]
                LDe = LDe + [Di + 2*E]
                LR = LR + [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+rt.rt_cond_cili(Di,Di + 2*E,flmd((root+Ti)/2))+R_rev+((((rt.rt_conv_cili(Di + 2*E,ctc.hc_m_ch(U,Di + 2*E,root,Ta)))**(-1))+((rt.rt_crad_cili(Di + 2*E,ctc.hr(eps,root,Ta)))**(-1)))**(-1)))/z]
        for flmd in LLMD3:
            for E in LE3:
                err = errf(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd, R_rev)
                root = optimize.brentq(err, Ta, Ti)
                LTe = LTe + [root]
                Lslmd = Lslmd + [flmd((root + Ti)/2)]
                Lq = Lq + [(ctc.qc_m_ch(U, Di + 2*E, root, Ta) + ctc.qr(eps, root, Ta))*(np.pi*(Di + 2*E))]
                LDe = LDe + [Di + 2*E]
                LR = LR + [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+rt.rt_cond_cili(Di,Di + 2*E,flmd((root+Ti)/2))+R_rev+((((rt.rt_conv_cili(Di + 2*E,ctc.hc_m_ch(U,Di + 2*E,root,Ta)))**(-1))+((rt.rt_crad_cili(Di + 2*E,ctc.hr(eps,root,Ta)))**(-1)))**(-1)))/z]
        
    if ((U != 0) and (H != 0)):
        
        errf0 = generate_err_tubes_for_v_si
        
        err0 = errf0(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps)
        
        root0 = optimize.brentq(err0, Ta, Ti)
        
        LTe = [root0]
        
        Lq = [(ctc.qc_m_cv(U, Di, H, root0, Ta) + ctc.qr(eps, root0, Ta))*(np.pi*(Di))]
        
        LR = [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+((((rt.rt_conv_cili(Di,ctc.hc_m_cv(U,Di,H,root0,Ta)))**(-1))+((rt.rt_crad_cili(Di,ctc.hr(eps,root0,Ta)))**(-1)))**(-1)))/z]
        
        errf = generate_err_tubes_for_v
        
        for flmd in LLMD1:
            for E in LE1:
                err = errf(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd, R_rev)
                root = optimize.brentq(err, Ta, Ti)
                LTe = LTe + [root]
                Lslmd = Lslmd + [flmd((root + Ti)/2)]
                Lq = Lq + [(ctc.qc_m_cv(U, Di + 2*E, H, root, Ta) + ctc.qr(eps, root, Ta))*(np.pi*(Di + 2*E))]
                LDe = LDe + [Di + 2*E]
                LR = LR + [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+rt.rt_cond_cili(Di,Di + 2*E,flmd((root+Ti)/2))+R_rev+((((rt.rt_conv_cili(Di + 2*E,ctc.hc_m_cv(U,Di + 2*E,H,root,Ta)))**(-1))+((rt.rt_crad_cili(Di + 2*E,ctc.hr(eps,root,Ta)))**(-1)))**(-1)))/z]
        for flmd in LLMD2:
            for E in LE2:
                err = errf(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd, R_rev)
                root = optimize.brentq(err, Ta, Ti)
                LTe = LTe + [root]
                Lslmd = Lslmd + [flmd((root + Ti)/2)]
                Lq = Lq + [(ctc.qc_m_cv(U, Di + 2*E, H, root, Ta) + ctc.qr(eps, root, Ta))*(np.pi*(Di + 2*E))]
                LDe = LDe + [Di + 2*E]
                LR = LR + [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+rt.rt_cond_cili(Di,Di + 2*E,flmd((root+Ti)/2))+R_rev+((((rt.rt_conv_cili(Di + 2*E,ctc.hc_m_cv(U,Di + 2*E,H,root,Ta)))**(-1))+((rt.rt_crad_cili(Di + 2*E,ctc.hr(eps,root,Ta)))**(-1)))**(-1)))/z]
        for flmd in LLMD3:
            for E in LE3:
                err = errf(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd, R_rev)
                root = optimize.brentq(err, Ta, Ti)
                LTe = LTe + [root]
                Lslmd = Lslmd + [flmd((root + Ti)/2)]
                Lq = Lq + [(ctc.qc_m_cv(U, Di + 2*E, H, root, Ta) + ctc.qr(eps, root, Ta))*(np.pi*(Di + 2*E))]
                LDe = LDe + [Di + 2*E]
                LR = LR + [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+rt.rt_cond_cili(Di,Di + 2*E,flmd((root+Ti)/2))+R_rev+((((rt.rt_conv_cili(Di + 2*E,ctc.hc_m_cv(U,Di + 2*E,H,root,Ta)))**(-1))+((rt.rt_crad_cili(Di + 2*E,ctc.hr(eps,root,Ta)))**(-1)))**(-1)))/z]
        
    if ((U == 0) and (H == 0)):
        
        errf0 = generate_err_tubes_nat_h_si
        
        err0 = errf0(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps)
        
        root0 = optimize.brentq(err0, Ta, Ti)
        
        LTe = [root0]
        
        Lq = [(ctc.qc_n_ch(U, Di, root0, Ta) + ctc.qr(eps, root0, Ta))*(np.pi*(Di))]
        
        LR = [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+((((rt.rt_conv_cili(de,ctc.hc_n_ch(U,de,root0,Ta)))**(-1))+((rt.rt_crad_cili(de,ctc.hr(eps,root0,Ta)))**(-1)))**(-1)))/z]
        
        errf = generate_err_tubes_nat_h
        
        for flmd in LLMD1:
            for E in LE1:
                err = errf(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd, R_rev)
                root = optimize.brentq(err, Ta, Ti)
                LTe = LTe + [root]
                Lslmd = Lslmd + [flmd((root + Ti)/2)]
                Lq = Lq + [(ctc.qc_n_ch(U, Di + 2*E, root, Ta) + ctc.qr(eps, root, Ta))*(np.pi*(Di + 2*E))]
                LDe = LDe + [Di + 2*E]
                LR = LR + [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+rt.rt_cond_cili(Di,Di + 2*E,flmd((root+Ti)/2))+R_rev+((((rt.rt_conv_cili(Di + 2*E,ctc.hc_n_ch(U,Di + 2*E,root,Ta)))**(-1))+((rt.rt_crad_cili(Di + 2*E,ctc.hr(eps,root,Ta)))**(-1)))**(-1)))/z]
        for flmd in LLMD2:
            for E in LE2:
                err = errf(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd, R_rev)
                root = optimize.brentq(err, Ta, Ti)
                LTe = LTe + [root]
                Lslmd = Lslmd + [flmd((root + Ti)/2)]
                Lq = Lq + [(ctc.qc_n_ch(U, Di + 2*E, root, Ta) + ctc.qr(eps, root, Ta))*(np.pi*(Di + 2*E))]
                LDe = LDe + [Di + 2*E]
                LR = LR + [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+rt.rt_cond_cili(Di,Di + 2*E,flmd((root+Ti)/2))+R_rev+((((rt.rt_conv_cili(Di + 2*E,ctc.hc_n_ch(U,Di + 2*E,root,Ta)))**(-1))+((rt.rt_crad_cili(Di + 2*E,ctc.hr(eps,root,Ta)))**(-1)))**(-1)))/z]
        for flmd in LLMD3:
            for E in LE3:
                err = errf(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd, R_rev)
                root = optimize.brentq(err, Ta, Ti)
                LTe = LTe + [root]
                Lslmd = Lslmd + [flmd((root + Ti)/2)]
                Lq = Lq + [(ctc.qc_n_ch(U, Di + 2*E, root, Ta) + ctc.qr(eps, root, Ta))*(np.pi*(Di + 2*E))]
                LDe = LDe + [Di + 2*E]
                LR = LR + [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+rt.rt_cond_cili(Di,Di + 2*E,flmd((root+Ti)/2))+R_rev+((((rt.rt_conv_cili(Di + 2*E,ctc.hc_n_ch(U,Di + 2*E,root,Ta)))**(-1))+((rt.rt_crad_cili(Di + 2*E,ctc.hr(eps,root,Ta)))**(-1)))**(-1)))/z]
        
    if ((U == 0) and (H != 0)):
        
        errf0 = generate_err_tubes_nat_v_si
        
        err0 = errf0(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps)
        
        root0 = optimize.brentq(err0, Ta, Ti)
        
        LTe = [root0]
        
        Lq = [(ctc.qc_n_cv(U, H, root0, Ta) + ctc.qr(eps, root0, Ta))*(np.pi*(Di))]
        
        LR = [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+((((rt.rt_conv_cili(de,ctc.hc_n_cv(U,H,root0,Ta)))**(-1))+((rt.rt_crad_cili(de,ctc.hr(eps,root0,Ta)))**(-1)))**(-1)))/z]
        
        errf = generate_err_tubes_nat_v
        
        for flmd in LLMD1:
            for E in LE1:
                err = errf(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd, R_rev)
                root = optimize.brentq(err, Ta, Ti)
                LTe = LTe + [root]
                Lslmd = Lslmd + [flmd((root + Ti)/2)]
                Lq = Lq + [(ctc.qc_n_cv(U, H, root, Ta) + ctc.qr(eps, root, Ta))*(np.pi*(Di + 2*E))]
                LDe = LDe + [Di + 2*E]
                LR = LR + [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+rt.rt_cond_cili(Di,Di + 2*E,flmd((root+Ti)/2))+R_rev+((((rt.rt_conv_cili(Di + 2*E,ctc.hc_n_cv(U,H,root,Ta)))**(-1))+((rt.rt_crad_cili(Di + 2*E,ctc.hr(eps,root,Ta)))**(-1)))**(-1)))/z]
        for flmd in LLMD2:
            for E in LE2:
                err = errf(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd, R_rev)
                root = optimize.brentq(err, Ta, Ti)
                LTe = LTe + [root]
                Lslmd = Lslmd + [flmd((root + Ti)/2)]
                Lq = Lq + [(ctc.qc_n_cv(U, H, root, Ta) + ctc.qr(eps, root, Ta))*(np.pi*(Di + 2*E))]
                LDe = LDe + [Di + 2*E]
                LR = LR + [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+rt.rt_cond_cili(Di,Di + 2*E,flmd((root+Ti)/2))+R_rev+((((rt.rt_conv_cili(Di + 2*E,ctc.hc_n_cv(U,H,root,Ta)))**(-1))+((rt.rt_crad_cili(Di + 2*E,ctc.hr(eps,root,Ta)))**(-1)))**(-1)))/z]
        for flmd in LLMD3:
            for E in LE3:
                err = errf(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd, R_rev)
                root = optimize.brentq(err, Ta, Ti)
                LTe = LTe + [root]
                Lslmd = Lslmd + [flmd((root + Ti)/2)]
                Lq = Lq + [(ctc.qc_n_cv(U, H, root, Ta) + ctc.qr(eps, root, Ta))*(np.pi*(Di + 2*E))]
                LDe = LDe + [Di + 2*E]
                LR = LR + [(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube)+rt.rt_cond_cili(Di,Di + 2*E,flmd((root+Ti)/2))+R_rev+((((rt.rt_conv_cili(Di + 2*E,ctc.hc_n_cv(U,H,root,Ta)))**(-1))+((rt.rt_crad_cili(Di + 2*E,ctc.hr(eps,root,Ta)))**(-1)))**(-1)))/z]
            
    #Diâmetros externos em milímetros.
    LDe_Disp = [1000*D for D in LDe]
    
    #Lista de soluções para a temperatura na face externa em °C.
    Lte = [(Te - 273.15) for Te in LTe]
    
    if ts_max > 0:
        for i in range(len(Lq) - 1, 0, -1):
            if Lte[i] > ts_max:
                del Lq[i]
                del LNM[i]
                del LE_Disp_DF[i]
                del LE_Disp_Imp_DF[i]
                del LDe_Disp[i]
                del Lte[i]
                del Lslmd[i]
                del LR[i]
    
    if Dt_max > 0 and not(fase_change):
        LVT = list(map(lambda x: Ti - (Ta - (Ta - Ti)*np.exp(-1/(m*c_or_h*x))),LR))
        for i in range(len(LVT) - 1, 0, -1):
            if LVT[i] > Dt_max:
                del LVT[i]
                del Lq[i]
                del LNM[i]
                del LE_Disp_DF[i]
                del LE_Disp_Imp_DF[i]
                del LDe_Disp[i]
                del Lte[i]
                del Lslmd[i]
                del LR[i]
    
    #Organização dos dados em um DataFrame.
    Disp = pd.DataFrame({'Material' : LNM,
                         'Espessura \n [mm]' : LE_Disp_DF,
                         'Espessura \n [pol]' : LE_Disp_Imp_DF,
                         'Diâmetro \n Externo [mm]' : LDe_Disp,
                         'Temperatura na \n Face Externa [°C]' : Lte,
                         'Condutividade Térmica \n do Isolante [W/(m.k)]': Lslmd})
    
    if (not(fase_change)):
        LVT = list(map(lambda x: Ti - (Ta - (Ta - Ti)*np.exp(-1/(m*c_or_h*x))),LR))
        Disp['Variação de Temperatura \n do Fluido [°C]'] = LVT
        LTS = list(map(lambda x: Ti - 273.15 - x, LVT))
        Disp['Temperatura do Fluido \n na Saída [°C]'] = LTS
    
    if (fase_change):
        LFC = list(map(lambda x: (x*z)/c_or_h, Lq))
        Disp['Formação de \n Condensado [kg/s]'] = LFC
    
    Disp['Fluxo de Calor \n [W/m]'] = Lq
    
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