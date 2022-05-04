import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate, interpolate


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

#Atribuindo as colunas da equação de estado e
#Passando de fm^-4 pra km^-2
densidade_barionica = eos[0]
densidade_barionica = np.multiply(densidade_barionica, CONVBARYON)
densidade_energia = eos[1]
densidade_energia = np.multiply(densidade_energia, CONV)
pressao = eos[2]
pressao = np.multiply(pressao, CONV)


Spline_barionica = interpolate.interp1d(pressao, densidade_barionica, kind='cubic', fill_value="extrapolate")
Spline_energia = interpolate.interp1d(pressao, densidade_energia, kind='cubic', fill_value="extrapolate")


pressao_max = pressao[-1]
pressao_min = pressao[0]
ec = densidade_energia[0]
bc = densidade_barionica[0]
raio = 0.01
dr = 0.01
dp = 0.01*CONV



def Tov(raio, y):
    dydr = [0, 0, 0]
    dydr[0] = 4.0*pi*raio**2*Spline_energia(ci[0])  # dm/dr
    dydr[1] = -(Spline_energia(ci[1]) + ci[1]) * (ci[0] + 4*pi*raio*raio*raio*ci[1])/(raio*raio - 2*raio*ci[0])  # dpr/dr
    dydr[2] = Spline_barionica(ci[1])*4.0*pi*raio*raio/(np.sqrt(1-2*ci[0]/raio))  # d(baryonic density)/dr

    return dydr
i = 0 
while pressao_min < pressao_max:
    ci = [4.0/3.0*pi*dr*dr*dr*ec, pressao_min, 4.0/3.0*pi*dr*dr*dr*bc]
    Tov_resolved = integrate.solve_ivp(Tov,(raio, pressao_min>=pressao_max),ci, method='RK45', raio_eval=dr)
    dmdr, dprdr, dbdr = Tov_resolved.y
    pressao_min = pressao_min + dp
    ec = Spline_energia(pressao_min)
    bc = Spline_barionica(pressao_min)
    i+=1

massa_max = np.amax(dmdr)
pos = np.where(massa_max)
dprdr_max = dprdr[pos]
dbdr_max = dbdr[pos]
print(massa_max, dbdr_max, dbdr_max)
