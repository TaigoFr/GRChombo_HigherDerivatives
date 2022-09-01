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
z_symmetry = True

fontsize = 18
fontsizeBig = 21

#################################################################

def get_AH_radius():
    print("Reading AH data")
    stats_AH1 = np.loadtxt(location + "data/stats_AH1.dat")
    hasAH3 = False
    try:
        stats_AH3 = np.loadtxt(location + "data/stats_AH3.dat")
        hasAH3 = True
    except:
        stats_AH3 = []
    times = []
    radii = []
    i1 = 0
    i3 = 0
    while True:
        if i1 < len(stats_AH1):
            row1 = stats_AH1[i1]
            time1 = row1[0]
        else:
            if hasAH3:
                time1 = 1e10 # big number
            else:
                break # finished
        if i3 < len(stats_AH3):
            row3 = stats_AH3[i3]
            time3 = row3[0]
        else:
            if hasAH3:
                break
            else:
                time3 = 1e10 # big number

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

# or use a priori function given merger time at X (1200 below)
# AH_radius_vs_time = lambda t: 0.5 if t < 1200 else 1

punctures = np.loadtxt(location + "data/punctures.dat")

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

allNCC_p = []
allNCC_m = []
allNCC_p_z4 = []
allNCC_m_z4 = []
for i in range(0, len(ds), jump):
    file = ds[i]
    time = file.current_time
    
    puncture = punctures[punctures[:,0] == time][0][1:4]
    if z_symmetry:
        puncture[2] = 0 # force numeric 0
    radius = np.linalg.norm(puncture-center)

    # radial +
    point = center + (puncture-center) * (radius + AH_radius_vs_time(time)) / radius
    # radial -
    #point = center + (puncture-center) * (radius - AH_radius_vs_time(time)) / radius

    #tangent = np.array([-(puncture-center)[1],(puncture-center)[0],(puncture-center)[2]]) / radius
    # tangent +
    #point = puncture + tangent * AH_radius_vs_time(time)
    # tangent -
    # point = puncture - tangent * AH_radius_vs_time(time)

    ray = file.r[point]
    NCC_p = ray["NCC_plus"][0]
    NCC_m = ray["NCC_minus"][0]
    NCC_p_z4 = ray["NCC_Z4_plus"][0]
    NCC_m_z4 = ray["NCC_Z4_minus"][0]

    allNCC_p.append([time, NCC_p])
    allNCC_m.append([time, NCC_m])
    allNCC_p_z4.append([time, NCC_p_z4])
    allNCC_m_z4.append([time, NCC_m_z4])

[times, NCC_p] = np.transpose(allNCC_p)
[times2, NCC_m] = np.transpose(allNCC_m)
[times3, NCC_p_z4] = np.transpose(allNCC_p_z4)
[times4, NCC_m_z4] = np.transpose(allNCC_m_z4)

print("Times:")
print(list(times))
print("NCC_p", NCC_p)
print("NCC_m", NCC_m)
print("NCC_p_z4", NCC_p_z4)
print("NCC_m_z4", NCC_m_z4)

def do_plot(times, plots, labels, name_end):
    fig, ax = plt.subplots(figsize=(12, 9))
    for time, plot, label in zip(times, plots, labels):
        ax.plot(time, plot, label=label)
    
    ax.legend(fontsize=fontsize)

    minVal = np.min([np.min(plot) for plot in plots])
    maxVal = np.max([np.max(plot) for plot in plots])
    
    ax.set_ylim([minVal * (1.1 if minVal < 0 else 0.9), maxVal * (0.9 if maxVal < 0 else 1.1)])

    plt.xlabel("Time (M)", fontsize=fontsize)

    ax.tick_params(axis='both', which='major', labelsize=fontsize)

    time = file.current_time
    fig.suptitle(r'$NCC$', fontsize=fontsizeBig, position=(0.5,0.93))

    plt.savefig("NCC_AH" + name_end + ".png", bbox_inches = 'tight')
    plt.close()

do_plot([times, times2, times3, times4], [NCC_p, NCC_m, NCC_p_z4, NCC_m_z4], ["NCC_plus", "NCC_minus", "NCC_Z4_plus", "NCC_Z4_minus"], "")

