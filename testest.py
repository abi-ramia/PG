# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 15:33:58 2020

@author: Tiago
"""

import ctc
import numpy as np
import matplotlib.pyplot as plt
import tikzplotlib

# =============================================================================
# DESENHA NU VS U
# =============================================================================

Te = 273.15 + 120
Ta = 273.15 + 24
Lc = 0.254
LLc = np.linspace(0.254, 0.254, 101)
LU = np.linspace(0, 10, 101)

LN = ctc.hc_n_ch(LU, LLc, Te, Ta)
LF = ctc.hc_f_c(LU, Lc, Te, Ta)
LC = ctc.hc_m_ch(LU, Lc, Te, Ta)

LN = LN/LC[-1]
LF = LF/LC[-1]
LC = LC/LC[-1]

plt.plot(LU, LN, 'k--', label = '$Nu_{n}$')
plt.plot(LU, LF, 'k:', label = '$Nu_{f}$')
plt.plot(LU, LC, 'k-.', label = '$Nu_{c}$')

plt.legend()

plt.xlabel('Velocidade do Vento [$\\frac{m}{s}$]')
plt.ylabel('Número de Nusselt Normalizado [1]')

tikzplotlib.save("test.tex")