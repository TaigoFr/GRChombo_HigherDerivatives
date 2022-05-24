# on a cluster with paralellism, run as
# 'mpiexec -n 8 -hosts=i01r01c01s01 python C_slices.py'
# host is the node running the script, 8 is the number of cores
# I realized making this '8' bigger does not help that much

import numpy as np
import yt
import matplotlib.pyplot as plt
import os
from scipy.interpolate import make_interp_spline

#################################################################
# USER DATA

location = '../'
is_corner = False
jump = 1 # plot every 'jump' files

fontsize = 18
fontsizeBig = 21

#################################################################

def get_AH_radius():
    print("Reading AH data")
    stats_AH1 = np.loadtxt(location + "data/stats_AH1.dat")
    stats_AH3 = np.array([])
    times = []
    radii = []
    i1 = 0
    i3 = 0
    while True:
        if i1 < len(stats_AH1):
            row1 = stats_AH1[i1]
            time1 = row1[0]
        else:
            # time1 = 1e10 # big number
            break # finished
        if i3 < len(stats_AH3):
            row3 = stats_AH3[i3]
            time3 = row3[0]
        else:
            time3 = 1e10 # big number
            # break # finished

        time = min(time1,time3)
        print("On time", time, end="\r")

        rows = []
        if time1 <= time3:
            rows.append([1, row1])
            i1 += 1
        if time3 <= time1:
            rows.append([3, row3])
            i3 += 1

        compare_radius = []
        for index, row in rows:
            area = row[2]
            if np.isnan(area):
                continue
            file = int(row[1])
            coords = np.loadtxt(location + f"data/coords/coords_AH{index}_{file:06d}.dat")
            compare_radius.append( np.min(coords[:,2]) )

        ah_radius_min = np.max(compare_radius)
        times.append(time)
        radii.append(ah_radius_min)

    return times, radii

# use AH Finder
times, radii = get_AH_radius()
AH_radius_vs_time = make_interp_spline(times, radii, k=3)

# or use a priori function
# AH_radius_vs_time = lambda t: 1

yt.enable_parallelism()

# Loading dataset
def read_hdf5(location, is_corner=False):
    ds = yt.load(location+'*p_*.hdf5') # plot files

    center = np.array([0,0,0])
    if not is_corner:
        center = np.array(ds[0].domain_right_edge)
        center = center - np.max(center)/2

    print("CENTER = {}".format(center), flush=True)

    # ds = ds[0:4]

    return ds, center

ds, center = read_hdf5(location + "hdf5/", is_corner)

allC = []
allCphys = []
allCminusCphys = []
for i in range(0, len(ds), jump):
    file = ds[i]
    time = file.current_time
    
    ray = file.r[[center[0] + 0, center[1] + 0, center[2] + AH_radius_vs_time(time)]]
    C_point = abs(ray["C"][0])
    Cphys_point = abs(ray["Cphys"][0])

    allC.append([time, C_point])
    allCphys.append([time, Cphys_point])
    allCminusCphys.append([time, abs(Cphys_point / C_point - 1) * 100])

[times, C] = np.transpose(allC)
[times2, Cphys] = np.transpose(allCphys)
[times3, CminusCphys] = np.transpose(allCminusCphys)

print("Times:")
print(list(times))
print("C", C)
print("Cphys", Cphys)
print("CminusCphys", CminusCphys)

def do_plot(times, C, Cphys, use_log, name_end):
    fig, ax = plt.subplots(figsize=(12, 9))
    if len(Cphys)>0:
        ax.plot(times, Cphys, label=r"$C_{phys}$")
    ax.plot(times, C, label=r"$C$")
    
    ax.legend(fontsize=fontsize)

    minC = min(np.min(C), np.min(Cphys)) if len(Cphys)>0 else np.min(C)
    maxC = max(np.max(C), np.max(Cphys)) if len(Cphys)>0 else np.max(C)
    if use_log:
        ax.set_ylim([minC * 0.9 , maxC * 1.1])
    else:
        ax.set_ylim([minC * (1.1 if minC < 0 else 0.9), maxC * (0.9 if maxC < 0 else 1.1)])

    plt.xlabel("Time (M)", fontsize=fontsize)

    if use_log:
        plt.yscale("log")

    ax.tick_params(axis='both', which='major', labelsize=fontsize)

    time = file.current_time
    fig.suptitle(r'$C$ vs $C_{phys}$', fontsize=fontsizeBig, position=(0.5,0.93))

    plt.savefig("CvsCphys_AH" + name_end + ".png", bbox_inches = 'tight')
    plt.close()

do_plot(times, C, Cphys, False, "")
do_plot(times, C, Cphys, True, "_log")
do_plot(times, CminusCphys, [], False, "_diff")

