import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import RK45
from scipy.interpolate import interp1d

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
            col1.append(aux2[j][0])
            col2.append(aux2[j][1])
            col3.append(aux2[j][2])
        return col1, col2, col3


x = input(
    'Digite apenas o nome do arquivo das equações de estado (exemplo: eos147.dat): ')
eos = Eos(x)





def Tov(raio, eos):
    dydr = []
    dydr[0] = 4.0*pi*raio*raio*energyDensitySpline(eos[1])  # dm/dr
    dydr[1] = - (energyDensitySpline(eos[1]) + eos[1]) * (eos[0] +
                                                          4*pi*raio*raio*raio*y[1])/(raio*raio - 2*raio*eos[0])  # dpr/dr
    dydr[2] = baryonDensitySpline(
        eos[1])*4.0*pi*raio*raio/(np.sqrt(1-2*eos[0]/raio))  # d(baryonic density)/dr


Tov_resolved = RK45(Tov(0.1, eos), 0, eos[0, 0, 0])

