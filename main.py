# -*- coding: utf-8 -*-
import ib
import ctc

#Temperatura do fluido em K.
Ti = 300 + 273.15

#Temperatura do ar ambiente em K.
Ta = 24 + 273.15

#Diâmetro externo da tubulação em m.
Di = 0.300

#Altura da tubulação em m, caso vertical.
H = 7

#Velocidade do vento em m/s. Caso convecção natural, U = 0.
U = 3

#Emissividade da superfície. Pode ser importada de ctc.py.
eps = ctc.emissividade('Alumínio')

#Número de horas de funcionamento por ano.
N = 2500

#Custo do combustível em $/J.
F = 1e-9

#Eficiência do sistema de conversão de combustível.
eta = 0.90

#Vida útil do isolamento em anos.
n = 5

#Taxa de atratividade anual.
i = 0.15

#Taxa de crescimento diferenciado do custo de energia.
delta = 0.08

if __name__ == '__main__':
    result = ib.iso_tubes(Ti, Ta, Di, H, U, eps, N, F, eta, n, i, delta)
    result.to_excel("resultado.xlsx", sheet_name='Todos_os_Isolantes')