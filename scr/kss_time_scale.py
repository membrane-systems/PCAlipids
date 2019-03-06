import math
import numpy as np
from scipy.stats import ks_2samp, gaussian_kde
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Pool
import sys, os


def get_nearest_value(iterable, value):
    for idx, x in enumerate(iterable):
        if x > value:
            break
    B = idx
    if idx - 1 == -1:
        A = idx + 1
    else:
        A = idx - 1
    # print(A,B)
    return A, B


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
    if type(data) is np.ndarray:
        for value in data:
            dist[int(round((value - grid[0]) / delta))] += 1.
    else:
        dist[int(round((data - grid[0]) / delta))] += 1.
    return dist
    # return gaussian_kde(data)(grid)


def cum_dist(dist):
    sum_ = dist[0]
    cum_ = [dist[0]]
    for i in range(1, len(dist)):
        sum_ += dist[i]
        cum_.append(cum_[i - 1] + dist[i])
    for i in range(len(cum_)):
        cum_[i] /= sum_
    return cum_


def KSS(cum_1, cum_2):
    max_ = 0
    for i in range(len(cum_1)):
        max_ = max(max_, abs(cum_1[i] - cum_2[i]))
    return max_

    # return ks_2samp(data_1, data_2)[0]


def KDE_for_step(data, grid):
    delta = (grid[1] - grid[0])
    dist = [0.] * len(grid)
    bin_idx = []
    if type(data) is np.ndarray:
        for value in data:
            dist[int(round((value - grid[0]) / delta))] += 1.
            if int(round((value - grid[0]) / delta)) not in bin_idx:
                bin_idx.append(int(round((value - grid[0]) / delta)))
    else:
        dist[int(round((data - grid[0]) / delta))] += 1.
        bin_idx.append(int(round((data - grid[0]) / delta)))
    return dist, sorted(bin_idx)


def KSS_for_step(dist_ideal, dist_an, bin_idx):
    last_value = 0
    max_ = 0
    for idx in bin_idx:
        max_ = max(max_, abs(dist_an[idx] - dist_ideal[idx]), abs(last_value - dist_ideal[idx]))
        last_value = dist_an[idx]
    return max_


def grid_len(N_samples, max_power):
    return [1.5 ** i  for i in range(2, int(max_power) - 1)]


def calc(filename, N_lip, timestep):     #, N_lip, N_bins, N_samples, cutoff = 0):
    N_lip = 128
    N_bins = 51
    N_samples = 24
    cutoff = 0
    data = load_data(filename)
    gN = len(data)
    y = linspace(min(data), max(data), N_bins)
    dist_ideal = KDE(data, y)
    cum_ideal = cum_dist(dist_ideal)
    data = split_data_by_lip(data, N_lip)
    max_power = math.log(len(data[0]), 1.5)
    # print(max_power)
    KSS_time = []
    T = []
    buff_KSS = 0.
    for data_an in data:
        dist_an = KDE(data_an, y)
        cum_an = cum_dist(dist_an)
        buff_KSS += KSS(cum_ideal, cum_an)
    KSS_time.append(buff_KSS / N_lip)
    buff_KSS = 0.0
    T.append(len(data[0]) * timestep)
    for tau in grid_len(N_samples, max_power):
        for data_ in data:
            data_an, N_chunks, L_tau = split_data_by_chunks(data_, tau)
            for data_an_ in data_an:
                if math.floor(len(data_) / tau) * 2 > N_bins:
                    dist_an = KDE(data_an_, y)
                    cum_an = cum_dist(dist_an)
                    buff_KSS += KSS(cum_ideal, cum_an)
                else:
                    # dist_an, bin_idx = KDE_for_step(data_an_, y)
                    # cum_an = cum_dist(dist_an)
                    # buff_KSS += KSS_for_step(cum_ideal, cum_an, bin_idx)
                    dist_an = KDE(data_an_, y)
                    cum_an = cum_dist(dist_an)
                    buff_KSS += KSS(cum_ideal, cum_an)
        KSS_time.append(buff_KSS / N_lip / N_chunks)
        buff_KSS = 0.0
        T.append(L_tau * timestep)

    
    # data = load_data(filename)
    # gN = len(data)
    # print(gN)
    # data = split_data_by_lip(data, N_lip)

    buff_KSS = 0.
    for data_lip in data:
        for i in range(len(data_lip)):
            dist_an, bin_idx = KDE_for_step(data_lip[i], y)
            cum_an = cum_dist(dist_an)
            buff_KSS += KSS_for_step(cum_ideal, cum_an, bin_idx)
    if buff_KSS / gN > 0.75:
        KSS_time.append(0.75)
    else:
        KSS_time.append(buff_KSS / gN)
    T.append(timestep)


    print(filename.split('/')[-1] + ' - processed')   
    # print(KSS_time, T) 
    return KSS_time, T

def main(file_name, range_fnames, N_lips, timestep, file_out):
    PATH = os.getcwd() + '/'
    # if both file or range of files are not defined 
    if not file_name and not range_fnames:
        print("The projections files have to be provided.\n\
        Use pcalipids.py projdist -h for help")
    # if range is defined
    elif range_fnames:
        # define file name for the first projection file in list
        firstFile = range_fnames[:range_fnames.find('-')]
        # define file name for the last projection file in list
        lastFile  = range_fnames[range_fnames.find('-')+1:]
        # find first PC
        first = int(firstFile[firstFile.rfind('_')+1:firstFile.rfind('.')])
        # find last PC
        last = int(lastFile[lastFile.rfind('_')+1:lastFile.rfind('.')])
        # define filemask
        mask = firstFile[:firstFile.rfind('_')]
        # define filerez
        rez = firstFile[firstFile.rfind('.'):]
        input_data = [(PATH + mask + "_" + str(i) + rez, \
            N_lips, timestep) \
            for i in range(first,last+1)]
    else:
        input_data = [(file_name, N_lips, timestep)]

    name = file_out[:file_out.rfind('.')]
    rez  = file_out[file_out.rfind('.'):]

    with Pool(8) as p:
        data = p.starmap(calc, input_data)
        for idx, obj in  enumerate(data[:len(input_data)]):
            T = obj[1]
            KSS_time = obj[0]
            plt.loglog(T, KSS_time, color = [0, 0, 1 - idx / len(input_data)])
    plt.ylim([0.005, 0.75])
    plt.xlabel('Time (ns)')
    plt.ylabel('K-S statistics')
    plt.savefig(name+'_values_vs_t'+'.png')
    plt.clf()

    file = open(name+'_values_vs_t'+rez, 'w')
    for j in range(len(data[0][0])):
        for i in range(len(data)):
            if i==0:
                file.write(str(data[i][1][j]) + ' ' + str(str(data[i][0][j])))
            else:
                file.write(' ' + str(str(data[i][0][j])))
        file.write("\n")
    file.close()

    file = open(name + '_relax_times_vs_pc' + rez, 'w')
    file.write('### KSS ###\n')
    file.write('E**2\n')

    handles = []

    PC = []
    T_relax = []

    for i, value in enumerate(data[:]):
        T = np.log(value[1])
        KSS = np.log(value[0])
        A,B = get_nearest_value(KSS, np.log(0.75 * math.e ** (-2)))
        a = (T[B] - T[A]) / (KSS[B] - KSS[A])
        b = T[A] - a * KSS[A]
        t_relax = a * np.log(0.75 * math.e ** (-2)) + b
        PC.append(i + 1)
        T_relax.append(math.e ** t_relax)

    p = plt.loglog(PC, T_relax, label = r'$\tau_2$',\
     color = 'blue', linestyle = '--')
    handle, = p 
    handles.append(handle)

    for i in range(len(PC)):
        file.write(str(PC[i]) + ' ' + str(T_relax[i]))
        file.write('\n')
    file.write('\n')

    PC = []
    T_relax = []

    for i, value in enumerate(data[:]):
        T = np.log(value[1])
        KSS = np.log(value[0])
        A, B = get_nearest_value(KSS, np.log(0.75 * math.e ** (-1)))
        a = (T[B] - T[A]) / (KSS[B] - KSS[A])
        b = T[A] - a * KSS[A]
        t_relax = a * np.log(0.75 * math.e ** (-1)) + b
        PC.append(i + 1)
        T_relax.append(math.e ** t_relax)

    file.write('E**1\n')
    for i in range(len(PC)):
        file.write(str(PC[i]) + ' ' + str(T_relax[i]))
        file.write('\n')
    file.close()


    p = plt.loglog(PC, T_relax, label = r'$\tau_1$', color = 'blue',\
     linestyle = '-')
    handle, = p

    handles.append(handle)
    plt.ylim([0.1, 10**4])
    plt.ylabel('Relaxation time (ns)')
    plt.xlabel('Component')
    plt.legend(handles = handles)
    plt.savefig(name + '_relax_times_vs_pc.png')
