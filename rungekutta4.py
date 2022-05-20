import numpy as np


def rk4(dmdr, dprdr, dbdr, ci, dr, raio, pressao_min):
    massa = ci[0]
    pressao = ci[1]
    d_bar = ci[2]
    while pressao >= pressao_min:
        massa = dmdr(pressao, raio)
        pressao = dprdr(massa, pressao, raio)
        d_bar = dbdr(massa, pressao, raio)
        raio = raio + dr
    return np.array([massa, d_bar, raio])
