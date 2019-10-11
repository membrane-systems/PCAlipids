import math
import numpy as np
from scipy.stats import ks_2samp, gaussian_kde
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Pool
import sys, os
import time
from numba import jit
import mdtraj as md

def sample(data, N_bins): #data to integers from 0 to N_bins-1
    m, M = np.min(data), np.max(data)
    x = (N_bins-0.0000001)/(M - m)
    return np.array((data - m)*x, dtype = int)

def load_data(filename):
    return np.loadtxt(filename, dtype=np.float64, comments=('@', '#', '&'))

def grid(tt, N_lip, timestep): #len(grids) = max_power, gr[0] - L_tau, gr[1] - number of parts, gr[2] - cut_fr
    max_power = int(math.log(tt, 1.5))
    grids = np.zeros((max_power, 3), dtype='int')
    for i in range(max_power):
        gr = grids[i]
        gr[0] = tt / 1.5**i
        N_chunks = tt // gr[0]
        gr[2] = gr[0] * N_chunks
        gr[1] = N_chunks * N_lip
    timescale = grids[:, 0] * timestep
    return grids, timescale

def prepare_data(sdat, N_bins, N_lip, fr):  #after sampling converts value in data to [cum_ideal[value], cum_ideal[value-1] or 0]
    sdat = sample(sdat, N_bins)
    N = len(sdat)
    t=time.clock()
    M = np.zeros((N, 2))
    cum = 0
    for i in range(N_bins):
        cuma = np.mean(sdat<=i) #cum_ideal
        bol = sdat==i
        M[bol, 0] = cum
        M[bol, 1] = cuma
        cum = cuma
    print('array prepared', time.clock() - t) 
    return M.reshape((N_lip, fr, 2))


def calc(M, grids):     #calculates KSS for prerared array M
    t=time.clock()
    max_power = len(grids)
    KSS_time = np.zeros(max_power)
    for g in range(max_power):
        gr = grids[g]
        M_an = M[:, :gr[2], :].reshape((gr[1], gr[0], 2))*gr[0]
        M_an = np.sort(M_an, axis=1)
        M_an[:, :, 0] = np.subtract(M_an[:, :, 0], np.arange(gr[0]))
        M_an[:, :, 1] = 1 - np.subtract(M_an[:, :, 1], np.arange(gr[0]))
        k = np.amax(M_an, axis=(1, 2))
        KSS_time[g] = np.mean(k)/gr[0]
    print('KSS done', time.clock() - t) 
    return KSS_time

def produce(filename, N_lip, fr, grids, N_bins): #loads, prepares and calculates KSS
    return calc(prepare_data(load_data(filename), N_bins, N_lip, fr), grids)

def get_nearest_value(iterable, value):
    for idx, x in enumerate(iterable):
        if x > value:
            break
    B = idx
    if idx - 1 == -1:
        A = idx + 1
    else:
        A = idx - 1
    return A, B


def plotter(data, timescale, power, n_PC):
    T_relax = []

    for i, value in enumerate(data):
        T = np.log(timescale)
        KSS = np.log(value)
        c = np.log(0.75 * math.e ** (-power))
        A,B = get_nearest_value(KSS, c)
        a = (T[B] - T[A]) / (KSS[B] - KSS[A])
        b = T[A] - a * KSS[A]
        t_relax = a * c + b
        T_relax.append(math.e ** t_relax)

    p = plt.loglog(np.arange(1, n_PC + 1), T_relax, label = r'$\tau_'+str(power)+'$',\
     color = 'blue', linestyle = '-'*power)
    handle, = p 
    return(T_relax, handle)

def main(file_name, range_fnames, fr, N_lip, N_bins, traj, timestep, file_out):
    if N_lip is None:
        N_lip = md.load(traj).n_residues #calculate N_lip from pdb file
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
        n_PC = last - first + 1
        filenames = [PATH + mask + "_" + str(i) + rez for i in range(first,last+1)]
    else:
        n_PC = 1
        filenames = [file_name]
        firstFile = file_name
    if fr is None:
        fr = len(load_data(firstFile))//N_lip #calculate fr from projection file, very slow!
    grids, timescale = grid(fr, N_lip, timestep) #one for all projections

    name = file_out[:file_out.rfind('.')]
    rez  = file_out[file_out.rfind('.'):]
    max_power = len(grids)

    input_data = [(filenames[i], N_lip, fr, grids, N_bins) for i in range(n_PC)] #prepare to Pool
    with Pool(8) as p:
        data = p.starmap(produce, input_data) #parallel calculations
    for idx, KSS_time in  enumerate(data):
        plt.loglog(timescale, KSS_time, color = [0, 0, 1 - idx / n_PC])
    ttt=time.clock()
    plt.ylim([0.005, 0.75])
    plt.xlabel('Time (ns)')
    plt.ylabel('K-S statistics')
    plt.savefig(name+'_values_vs_t'+'.png')
    plt.clf()
    np.savetxt(name+'_values_vs_t'+rez, np.vstack((timescale, data)).T, fmt='%-8.3f', header='# time in first column, KSS in other columns', footer='&', comments='')

    handles = [0, 0]
    T_relax = np.zeros((3, n_PC))
    T_relax[0] = np.arange(1, n_PC+1)
    T_relax[1], handles[0] = plotter(data, timescale, 2, n_PC)
    T_relax[2], handles[1] = plotter(data, timescale, 1, n_PC)
    plt.ylim([0.1, 10**4])
    plt.ylabel('Relaxation time (ns)')
    plt.xlabel('Component')
    plt.legend(handles = handles)
    plt.savefig(name + '_relax_times_vs_pc.png')

    np.savetxt(name + '_relax_times_vs_pc' + rez, T_relax.T, fmt='%-8.3f', header='# PC in first column, E**2 in second column, E**1 in third column', footer='&', comments='')
    print('finished', time.clock() - ttt)
