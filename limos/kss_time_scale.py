import math
import numpy as np
from scipy.stats import ks_2samp, gaussian_kde
import matplotlib.pyplot as plt
from multiprocessing import Pool
import time


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
    delta = (grid[1] - grid[0])
    dist = [0.] * len(grid)
    for value in data:
        dist[int(round((value - grid[0]) / delta))] += 1.
    return dist
    # return gaussian_kde(data)(grid)


def KSS(data_1, data_2):
    sum_1 = data_1[0]
    cum_1 = [data_1[0]]
    for i in range(1, len(data_1)):
        sum_1 += data_1[i]
        cum_1.append(cum_1[i - 1] + data_1[i])
    for i in range(len(cum_1)):
        cum_1[i] /= sum_1


    sum_2 = data_2[0]
    cum_2 = [data_2[0]]
    for i in range(1, len(data_2)):
        sum_2 += data_2[i]
        cum_2.append(cum_2[i - 1] + data_2[i])
    for i in range(len(cum_2)):
        cum_2[i] /= sum_2


    max_ = -1
    for i in range(len(cum_1)):
        max_ = max(max_, abs(cum_1[i] - cum_2[i]))
    # print(max_)
    # print(data_1)
    # print(data_2)
    # print(sum_1)
    # print(sum_2)
    # plt.plot(cum_1)
    # plt.plot(cum_2)
    # plt.show()
    return max_

    # return ks_2samp(data_1, data_2)[0]


def grid_len(N_samples, cutoff):
    return [1.5 ** i  for i in [2, 4, 6, 10, 14, 16, 18, 20, 21, 22, 23]][cutoff:] # [2, 4, 6, 10, 14, 16, 18, 20]


def main(filename):     #, N_lip, N_bins, N_samples, cutoff = 0):
    N_lip = 128
    N_bins = 51
    N_samples = 24
    cutoff = 0
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
    buff_KSS = 0.0
    T.append(len(data[0]) * 0.01)
    # print(KSS_time, T)
    # print(grid_len(N_samples, cutoff))
    for tau in grid_len(N_samples, cutoff):
        # print(tau)
        for data_ in data:
            data_an, N_chunks, L_tau = split_data_by_chunks(data_, tau)
            if tau > 1:
                for data_an_ in data_an:
                    dist_an = KDE(data_an_, y)
                    buff_KSS += KSS(dist_ideal, dist_an)
        KSS_time.append(buff_KSS / N_lip / N_chunks)
        buff_KSS = 0.0
        T.append(L_tau * 0.01)
    print(filename)
    return KSS_time, T



filenames = []
N = 100
for i in range(1, N + 1):
    filenames.append('proj_%s.xvg' % i)


TIME = time.time()

with Pool(8) as p:
    data = p.map(main, filenames)
    # print(data)
    for idx, obj in  enumerate(data):
        T = obj[1]
        KSS_time = obj[0]
        plt.loglog(T, KSS_time, color = [0, 0, 1 - idx / len(data)])
# file.close()
print(time.time() - TIME)
plt.ylim([0.005, 1])
plt.xlabel('Time (ns)')
plt.ylabel('K-S statistics')
plt.show()
