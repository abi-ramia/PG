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
Ti = 240 + 273.15

#Temperatura do ar ambiente em K.
Ta = 24 + 273.15

#Coeficiente de convecção do escoamento interno à tubulação em SI.
h_fld = 10

#Condutividade térmica do material da tubulaçãoem SI.
lmd_tube = 60

#Velocidade do vento em m/s. Caso convecção natural, U = 0.
U = 2.4

#Altura da tubulação em m, caso vertical.
H = 0

#comprimento da tubulação em m.
z = 1500

#Emissividade da superfície. Pode ser importada de ctc.py.
eps = ctc.emissividade('Alumínio')

#Mudança ou não de fase.
fase_change = False

#Vazão mássica de fluido em kg/s.
m = 35

#Calor específico [J/(kg.K)] se não houver mudança de fase; entalpia de vaporizaçao se houver [J/kg].
c_or_h = 2090

if __name__ == '__main__':
    result = ib.iso_tubes(di, de, Ti, Ta, h_fld, lmd_tube, U, H, z, eps, fase_change, m, c_or_h)
    #result.to_excel("resultado.xlsx", sheet_name='Todos_os_Isolantes')