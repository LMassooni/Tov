import numpy as np


def rk4(dmdr, dprdr, dbdr, ci, dr, raio, pressao_min):
    d_en = ci[0]
    pressao = ci[1]
    d_bar = ci[2]
    while pressao >= pressao_min:
        massa = dmdr(d_en, raio)
        d_en = dprdr(massa, d_en, raio)
        print(d_en)
        d_bar = dbdr(massa, d_en, raio)
        raio = raio + dr
        pressao = abs(d_en)
    return np.array([massa, d_bar, raio])
