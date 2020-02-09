# -*- coding: utf-8 -*-
import ib
import ctc
import rt
import ar

#Diâmetro interno da tubulação em m.
di = 0.254

#Diâmetro externo da tubulação em m.
de = 0.300

#Temperatura do fluido em K.
Ti = 250 + 273.15

#Temperatura do ar ambiente em K.
Ta = 24 + 273.15

#Coeficiente de convecção do escoamento interno à tubulação em SI.
h_fld = 75

#Condutividade térmica do material da tubulaçãoem SI.
lmd_tube = 60

#Velocidade do vento em m/s. Caso convecção natural, U = 0.
U = 6

#Altura da tubulação em m, caso vertical.
H = 0

#comprimento da tubulação em m.
z = 1000

#Emissividade da superfície. Pode ser importada de ctc.py.
eps = ctc.emissividade('Alumínio')

#Mudança [True] ou não [False] de fase.
fase_change = False

#Vazão mássica de fluido em kg/s.
m = 35

#Calor específico [J/(kg.K)] se não houver mudança de fase; entalpia de vaporizaçao se houver [J/kg].
c_or_h = 4180

#Temperatura máxima admissível da superfície externa [°C]. Se == 0, ignorada.
ts_max = 0

#Variação máxima admissível de temperatura do fluido [°C]. Se == 0 ou houver mudança de fase, ignorada.
Dt_max = 0

#Resitência térmica da seção de revestimento protetivo em m.K/W.
R_rev = 1e-5

#Taxa máxima admissível de transferência de calor em W. Se == 0, ignorada.
Q_max = 0

(N, F, eta, n, i, delta, tm) = (400, 3.5e-9, 0.80, 10, 0.15, 0.08, 0.02)

if True:
    result = ib.iso_tubes(di, de, Ti, Ta, h_fld, lmd_tube, U, H, z, eps, fase_change, m, c_or_h, ts_max, Dt_max, R_rev, Q_max, N, F, eta, n, i, delta, tm)
    result.to_excel("resultado.xlsx", sheet_name='IsoBrás', index = False, startcol = 0, freeze_panes = (2,2))