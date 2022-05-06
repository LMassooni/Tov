import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate, interpolate

                        ###################### COMO USAR ##########################
###             A única exigência para o código é instalar as bibliotecas numpy e scipy
###             com o comando pip install numpy, scipy, ou python3 pip install numpy, scipy
###             e ter o arquivo da equação de estado (.dat) na mesma pasta do código.


# Constantes e Conversões
pi = round(np.pi, 10)
CONV = 0.00026115  # 1 fm-4 = 2.6115e-4 km-2
# 939.0 * (1 MeV/fm3) = 939.0*1.3234e-6 km-2 = 1.2426726e-3 km-2
CONVBARYON = 0.0012426726
Massa_sol = 1.4766  # MSOL = 1.4766 km = 1.116e60 MeV

### Vai receber a lista de números da equação de estado e vai tratar os dados, jogando pra uma array com
### componente 0 = densidade barionica, componente 1 = densidade de energia e componente 2 = pressão
def Eos(eos):
    aux = []
    aux2 = []
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

## Montando as splines

Spline_barionica = interpolate.interp1d(pressao, densidade_barionica, kind='cubic', fill_value="extrapolate")
Spline_energia = interpolate.interp1d(pressao, densidade_energia, kind='cubic', fill_value="extrapolate")


pressao_max = pressao[-1]
pressao_min = pressao[0]
ec = densidade_energia[0]
bc = densidade_barionica[0]
raio = 0.01
dr = 0.01
dp = 0.01*CONV

### função com as equações de tov
def Tov(raio, y):
    dydr = [0, 0, 0, 0]
    dydr[0] = 4.0*pi*raio**2*Spline_energia(ci[0])  # dm/dr
    dydr[1] = -(Spline_energia(ci[1]) + ci[1]) * (ci[0] + 4*pi*raio*raio*raio*ci[1])/(raio*raio - 2*raio*ci[0])  # dpr/dr
    dydr[2] = Spline_barionica(ci[1])*4.0*pi*raio*raio/(np.sqrt(1-2*ci[0]/raio))  # d(baryonic density)/dr
    dydr[3] = raio

    return dydr
i = 0

### O while vai integrar a tov das menores pressões até as maiores da equação de estado
### ci são as condições iniciais, onde a variavel pressao_min ira aumentar de dp a cada estrela
### onde ci[0] = densidade de energia, ci[1] é a pressao inicial, ci[2] densidade barionica e ci[3] raio 
### integrate.solve_ivp irá receber como argumentos a função a se resolver
### os limites, raio (que foi definido anteriormente), até qnd a pressao minima chegar na pressao maxima
### as condições iniciais, o método, e raio_eval = tamanho do passo (dr)
### Tov_resolved.y recebe os resultados da integração e salva em 4 variaveis
### a cada integral resolvida, redefine as condições iniciais (pressao_min, bc e ec)
while pressao_min < pressao_max:
    ci = [4.0/3.0*pi*dr*dr*dr*ec, pressao_min, 4.0/3.0*pi*dr*dr*dr*bc, 0]
    Tov_resolved = integrate.solve_ivp(Tov,(raio, pressao_min>=pressao_max),ci, method='RK45', raio_eval=dr)
    dmdr, dprdr, dbdr, r = Tov_resolved.y
    pressao_min = pressao_min + dp
    ec = Spline_energia(pressao_min)
    bc = Spline_barionica(pressao_min)
    print("Estrela:", i)
    i+=1


### np.amax vai pegar o valor máximo da array dmdr. np.where vai identificar a posição do máximo
### para conseguir as características da estrela de massa maior.
massa_max = np.amax(dmdr)
m = massa_max/Massa_sol     ### Retornando as unidades para massa solares
pos = np.where(massa_max)
dprdr_max = dprdr[pos]
dbdr_max = dbdr[pos]
raio_max = r[pos]
db = dbdr_max/CONVBARYON    ### Retornando para fm^-4
d = dprdr_max/CONV 

print(m, d, db, raio_max)
