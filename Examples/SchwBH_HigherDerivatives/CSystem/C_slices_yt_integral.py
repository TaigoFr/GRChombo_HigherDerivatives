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

rmin = 0.
rmax = 5.
location = '../'
is_corner = False
jump = 10 # plot every 'jump' files

fontsize = 18
fontsizeBig = 21

# Plot several widths:
# width_min_deviation = -3 # 20 * 2^-3
# width_max_deviation =  1 # 20 * 2^1 

# Only 1 width:
width_min_deviation = 0
width_max_deviation = 0

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

def _C_minus_Cphys(field, data):
    return np.abs(data["C"] - data["Cphys"])

allAverageCminusCphys = []
allAverageCminusCphysWithoutAH = []
for i in range(0, len(ds), jump):
    file = ds[i]
    time = file.current_time
    file.add_field("C_minus_Cphys", _C_minus_Cphys, units = "")
    
    currentAverageCminusCphys = []
    currentAverageCminusCphysWithoutAH = []
    for p in np.arange(width_min_deviation, width_max_deviation + 1, 1):
        # compute integral differences between C and Cphys
        power = np.power(2., p)
        ray = file.r[[center[0] + 0, center[1] + 0, center[2] + 0]:[center[0] + 0, center[1] + 0, center[2] + rmax*power]]
        averageCminusCphys = ray.mean("C_minus_Cphys", weight="cell_volume")
        averageC = abs(ray.mean("C", weight="cell_volume"))
        currentAverageCminusCphys.append([rmax*power, averageCminusCphys / averageC * 100])

        # now without AH
        ray = file.r[[center[0] + 0, center[1] + 0, center[2] + AH_radius_vs_time(time)]:[center[0] + 0, center[1] + 0, center[2] + rmax*power]]
        averageCminusCphys = ray.mean("C_minus_Cphys", weight="cell_volume")
        averageC = abs(ray.mean("C", weight="cell_volume"))
        currentAverageCminusCphysWithoutAH.append([rmax*power, averageCminusCphys / averageC * 100])

    allAverageCminusCphys.append([float(time), np.array(currentAverageCminusCphys)])
    allAverageCminusCphysWithoutAH.append([float(time), np.array(currentAverageCminusCphysWithoutAH)])

[times_noAH, rsAndValues_noAH] = np.transpose(allAverageCminusCphys)
[times_withAH, rsAndValues_withAH] = np.transpose(allAverageCminusCphysWithoutAH)

print("Times:")
print(list(times_noAH))
print("Radius + CvsCphys (from r=0)")
print(rsAndValues_noAH if len(rsAndValues_noAH[0])>1 else [row[0,1] for row in rsAndValues_noAH])
print("Radius + CvsCphys (from r=r_AH)")
print(rsAndValues_withAH if len(rsAndValues_withAH[0])>1 else [row[0,1] for row in rsAndValues_withAH])

fig, ax = plt.subplots(figsize=(12, 9))

for i in range(len(rsAndValues_noAH[0])):
    values = [rsAndValues_noAHInTime[i][1] for rsAndValues_noAHInTime in rsAndValues_noAH]
    ax.plot(times_noAH, values, label=r"$r_{max}$ = " + f"{allAverageCminusCphys[0][1][i][0]}", linewidth=2, color=plt.cm.RdYlBu(i/len(rsAndValues_noAH[0])))

for i in range(len(rsAndValues_withAH[0])):
    values = [rsAndValues_withAHInTime[i][1] for rsAndValues_withAHInTime in rsAndValues_withAH]
    ax.plot(times_withAH, values, '--', linewidth=2, color=plt.cm.RdYlBu(i/len(rsAndValues_withAH[0])))

ax.legend(fontsize=fontsize)

plt.xlabel("Time (M)", fontsize=fontsize)
ax.tick_params(axis='both', which='major', labelsize=fontsize)

# fig.suptitle(r'$C$ vs $C_{phys}$ ( t = %.2fM )' % time, fontsize=fontsizeBig, position=(0.5,0.93))

plt.savefig("CvsCphys_integral", bbox_inches = 'tight')
plt.close()
