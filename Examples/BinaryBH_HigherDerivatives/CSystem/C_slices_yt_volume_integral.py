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

width = 20.
location = '../' # base folder with /data and /hdf5 subfolders
jump = 10 # plot every 'jump' files
z_symmetry = True

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

def _C_volume_without_AH(field, data):
    unit = file.length_unit.to('code_length')
    points = np.transpose([data["index", "x"] - unit * puncture[0], data["index", "y"] - unit * puncture[1], data["index", "z"] - unit * puncture[2]])
    points_norm = [np.linalg.norm(point) for point in points]
    points_invalid = np.array([norm < AH_radius_vs_time(time) and norm > width for norm in points_norm])
    data = data["C"]**2
    data[points_invalid] = 0;
    return data

def _Cphys_volume_without_AH(field, data):
    unit = file.length_unit.to('code_length')
    points = np.transpose([data["index", "x"] - unit * puncture[0], data["index", "y"] - unit * puncture[1], data["index", "z"] - unit * puncture[2]])
    points_norm = [np.linalg.norm(point) for point in points]
    points_invalid = np.array([norm < AH_radius_vs_time(time) and norm > width for norm in points_norm])
    data = data["Cphys"]**2
    data[points_invalid] = 0;
    return data

def _CvsCphys_volume_without_AH(field, data):
    unit = file.length_unit.to('code_length')
    points = np.transpose([data["index", "x"] - unit * puncture[0], data["index", "y"] - unit * puncture[1], data["index", "z"] - unit * puncture[2]])
    points_norm = [np.linalg.norm(point) for point in points]
    points_invalid = np.array([norm < AH_radius_vs_time(time) and norm > width for norm in points_norm])
    data = np.abs(data["C"] - data["Cphys"])**2
    data[points_invalid] = 0;
    return data

# to calculate volume
def _vol_eff(field, data):
    unit = file.length_unit.to('code_length')
    points = np.transpose([data["index", "x"] - unit * puncture[0], data["index", "y"] - unit * puncture[1], data["index", "z"] - unit * puncture[2]])
    points_norm = [np.linalg.norm(point) for point in points]
    points_invalid = np.array([norm < AH_radius_vs_time(time) and norm > width for norm in points_norm])
    chi = data["chi"]
    chi[points_invalid] = 0.
    chi[np.logical_not(points_invalid)] = 1.
    return chi

ds, center = read_hdf5(location + "hdf5/")

fontsize = 18
fontsizeBig = 21

allAverageCWithoutAH = []
allAverageCphysWithoutAH = []
allAverageCminusCphysWithoutAH = []
for i in range(0, len(ds), jump):
    file = ds[i]
    time = file.current_time
    
    puncture = punctures[punctures[:,0] == time][0][1:4]

    file.add_field("C_volume_without_AH", _C_volume_without_AH, units = "", sampling_type = "cell")
    file.add_field("Cphys_volume_without_AH", _Cphys_volume_without_AH, units = "", sampling_type = "cell")
    file.add_field("CvsCphys_volume_without_AH", _CvsCphys_volume_without_AH, units = "", sampling_type = "cell")
    file.add_field("vol_eff",_vol_eff, units = "")

    if z_symmetry:
        puncture[2] = 0 # force numeric 0
    radius = np.linalg.norm(puncture)

    high = file.domain_right_edge
    # overestimating the region required
    point_max = puncture / radius * (radius + width * 2)
    point_min = puncture / radius * (radius - width * 2)
    print(point_min, point_max)
    domain = file.r[point_min[0]:point_max[0], point_min[1]:point_max[1], point_min[2]:max(point_max[2],width * 2)]

    # percentage of domain without counting the inside of a 'chi' contour
    vol = domain.mean("vol_eff", weight="cell_volume")
    L2C = np.sqrt(domain.mean("C_volume_without_AH", weight="cell_volume"))
    L2Cphys = np.sqrt(domain.mean("Cphys_volume_without_AH", weight="cell_volume"))
    L2CvsCphys = np.sqrt(domain.mean("CvsCphys_volume_without_AH", weight="cell_volume"))

    L2C = L2C / vol
    L2Cphys = L2Cphys / vol
    L2CvsCphys = L2CvsCphys / vol

    allAverageCWithoutAH.append([float(time), L2C])
    allAverageCphysWithoutAH.append([float(time), L2Cphys])
    allAverageCminusCphysWithoutAH.append([float(time), L2CvsCphys])

[times_C, C] = np.transpose(allAverageCWithoutAH)
[times_Cphys, Cphys] = np.transpose(allAverageCphysWithoutAH)
[times_CvsCphys, CvsCphys] = np.transpose(allAverageCminusCphysWithoutAH)

print("Times:")
print(times_C)
print("C")
print(C)
print("Cphys")
print(Cphys)
print("CvsCphys")
print(CvsCphys)

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
    fig.suptitle(r'$C$ vs $C_{phys} (volume integral)$', fontsize=fontsizeBig, position=(0.5,0.93))

    plt.savefig("CvsCphys_volume" + name_end + ".png", bbox_inches = 'tight')
    plt.close()

do_plot(times, C, Cphys, False, "")
do_plot(times, C, Cphys, True, "_log")
do_plot(times, CminusCphys, [], False, "_diff")
