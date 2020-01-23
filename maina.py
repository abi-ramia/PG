# -*- coding: utf-8 -*-
import iba
import ctc
import wtr

#Diâmetro interno da tubulação em m.
di = 0.254

#Diâmetro externo da tubulação em m.
de = 0.300

#Temperatura do fluido em K.
Ti = 4 + 273.15

#Temperatura do ar ambiente em K.
Ta = 32 + 273.15

#Condutividade térmica do material da tubulaçãoem SI.
lmd_tube = 60

#Velocidade do vento em m/s. Caso convecção natural, U = 0.
U = 3

#Altura da tubulação em m, caso vertical.
H = 0

#comprimento da tubulação em m.
z = 500

#Emissividade da superfície. Pode ser importada de ctc.py.
eps = ctc.emissividade('Alumínio')

#Vazão mássica de água em kg/s.
m = 35

#Coeficiente de convecção do escoamento interno à tubulação em SI.
h_fld = wtr.hc(m, di, Ti, 0)

#Calor específico em J/(kg.K) da água. Pode ser importado de wtr.py.
c = wtr.cp(Ti)

#Variação máxima admissível de temperatura da água em °C. Se == 0.
Dt_max = 0

if True:
    result = iba.iso_tubes(di, de, Ti, Ta, h_fld, lmd_tube, U, H, z, eps, m, c, Dt_max)
    #result.to_excel("resultado.xlsx", sheet_name='Todos_os_Isolantes')