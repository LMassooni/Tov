import matplotlib.pyplot as plt
import numpy as np
from rungekutta4 import rk4
from scipy import interpolate

pi = round(np.pi, 10)
CONV = 0.00026115  # 1 fm-4 = 2.6115e-4 km-2
# 939.0 * (1 MeV/fm3) = 939.0*1.3234e-6 km-2 = 1.2426726e-3 km-2
CONVBARYON = 0.0012426726
Massa_sol = 1.4766  # MSOL = 1.4766 km = 1.116e60 MeV

aux = []
aux2 = []


def Eos(eos):
    with open(eos) as f:
        for line in f:
            x = line.strip()
            aux.append(x)
        for i in range(len(aux)):
            a = str(aux[i])
            b = a.split()
            aux2.append(b)
        col1 = []
        col2 = []
        col3 = []
        for j in range(len(aux2)):
            col1.append(float(aux2[j][0]))
            col2.append(float(aux2[j][1]))
            col3.append(float(aux2[j][2]))
        return col1, col2, col3


x = input(
    'Digite apenas o nome do arquivo da equação de estado (exemplo: eos147.dat): ')
eos = np.array(Eos(x))

# Atribuindo as colunas da equação de estado e
# Passando de fm^-4 pra km^-2
densidade_barionica = eos[0]
densidade_barionica = np.multiply(densidade_barionica, CONVBARYON)
densidade_energia = eos[1]
densidade_energia = np.multiply(densidade_energia, CONV)
pressao = eos[2]
pressao = np.multiply(pressao, CONV)


Spline_barionica = interpolate.interp1d(
    pressao, densidade_barionica, kind='cubic', fill_value="extrapolate")
Spline_energia = interpolate.interp1d(
    pressao, densidade_energia, kind='cubic', fill_value="extrapolate")


pressao_max = pressao[-1]
pressao_min = pressao[0]
ec = densidade_energia[0]
bc = densidade_barionica[0]
raio = 1
dr = 0.01
dp = 0.01*CONV


def dmdr(d_en, raio):

    massa = 4*pi*Spline_energia(d_en)*raio**2

    return massa


def dprdr(d_en, pressao, raio):

    pressao = - (Spline_energia(pressao) + pressao) * \
        (d_en+4*pi*raio**3*pressao)/(raio**2 - 2*raio*d_en)

    return pressao


def dbdr(d_en, pressao, raio):

    densidade = Spline_barionica(pressao)*4*pi*raio**2/(np.sqrt(1-2*d_en/raio))

    return densidade


i = 0
massa_ = []
d_en_ = []
d_bar_ = []
raio_ = []

while pressao_min < pressao_max:
    ci = [(4.0/3.0)*pi*dr*dr*dr*ec, pressao_min, (4.0/3.0)*pi*dr*dr*dr*bc]
    tov_resolved = rk4(dmdr, dprdr, dbdr, ci, dr, raio, pressao[0])
    massa_.append(tov_resolved[0])
    d_bar_.append(tov_resolved[1])
    raio_.append(tov_resolved[2])

    pressao_min = pressao_min + dp
    ec = Spline_energia(pressao_min)
    bc = Spline_barionica(pressao_min)
    print('Estrela:', i)
    i += 1

massa_max = (np.amax(massa_))
pos = massa_.index(massa_max)
massa_max = massa_max*Massa_sol
d_barionica = (d_bar_[pos])*CONVBARYON
raio_max = raio_[pos]

print(raio_)
