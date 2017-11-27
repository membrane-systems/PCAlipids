import math
import numpy as np
from scipy.stats import ks_2samp, gaussian_kde
import matplotlib.pyplot as plt
def load_data(filename):
    file = open(filename, 'r')
    data = []
    for line in file:
        if line.find('@') == -1 and line.find('&') == -1:
            #print(line.split())
            data.append(float(line.split()[0]))
        if line.find('&') != -1:
            break
    file.close()
    data = np.array(data, dtype = np.float64)
    return data


def split_data_by_lip(data, N_lip):
    data_spl = []
    for i in range(N_lip):
        data_spl.append(data[i * len(data) // N_lip:(i + 1) * len(data) // N_lip])
    return data_spl


def split_data_by_chunks(data, tau):
    L_tau = math.floor(len(data) / tau)
    N_chunks = math.floor(len(data) / L_tau)
    data_spl = []
    for i in range(N_chunks):
        data_spl.append(data[i * L_tau:(i + 1) * L_tau])
    return data_spl, N_chunks, L_tau


def linspace(min_, max_, N):
    return np.linspace(min_, max_, N)


def KDE(data, grid):
    return gaussian_kde(data)(grid)


def KSS(data_1, data_2):
    return ks_2samp(data_1, data_2)[0]


def grid_len(N_samples, cutoff):
    return [1.5 ** i  for i in [2, 4, 6, 10, 14, 16, 18, 20]][cutoff:]


def main(filename, N_lip, N_bins, N_samples, cutoff = 0):
    data = load_data(filename)
    y = linspace(min(data), max(data), N_bins)
    dist_ideal = KDE(data, y)
    data = split_data_by_lip(data, N_lip) # Разделение по липидам
    KSS_time = []
    T = []
    buff_KSS = 0.
    for data_an in data:
        dist_an = KDE(data_an, y)
        buff_KSS += KSS(dist_ideal, dist_an)
    KSS_time.append(buff_KSS / N_lip)
    T.append(len(data[0]) * 0.01)
    print(KSS_time, T)
    print(grid_len(N_samples, cutoff))
    for tau in grid_len(N_samples, cutoff):
        print(tau)
        for data_ in data:
            data_an, N_chunks, L_tau = split_data_by_chunks(data_, tau)
            if tau > 1:
                for data_an_ in data_an:
                    dist_an = KDE(data_an_, y)
                    buff_KSS += KSS(dist_ideal, dist_an)
        KSS_time.append(buff_KSS / N_lip / N_chunks)
        T.append(L_tau * 0.01)
        print(KSS_time, T)
    return KSS_time, T


file = open('KSS_vs_T.xvg', 'w')
for i in range(1, 51):
    KSS_time, T = main('/home/kmustafin/sim/1030GC36_DOPC_T_310/pca/proj_%s.xvg' % i, 128, 101, 24)
    file.write(str(KSS_time) + ' ' + str(T))
    file.write('\n')
    file.write('\n')
    plt.loglog(T, KSS_time, color = [0, 0, 1 - i / 100 * 2])
file.close()
plt.xlabel(Time (ns))
plt.ylabel('K-S statistics')
plt.show()