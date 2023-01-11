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
location1 = '../../test_v1_better_params_with_eps_ramp/' # always compare with this
location2 = '../'
jump = 1 # plot every 'jump' files
z_symmetry = True

# Plot several widths:
# width_min_deviation = -3 # 20 * 2^-3
# width_max_deviation =  1 # 20 * 2^1 

# Only 1 width:
width_min_deviation = 0
width_max_deviation = 0

#################################################################

def get_AH_radius(location):
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
times1, radii1 = get_AH_radius(location1)
times2, radii2 = get_AH_radius(location2)
AH_radius_vs_time1 = make_interp_spline(times1, radii1, k=3)
AH_radius_vs_time2 = make_interp_spline(times2, radii2, k=3)

# or use a priori function
# AH_radius_vs_time = lambda t: 1

AH_centers1 = np.loadtxt(location1 + "data/stats_AH1.dat")
AH_centers2 = np.loadtxt(location2 + "data/stats_AH1.dat")

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

ds1, center1 = read_hdf5(location1 + "hdf5/")
ds2, center2 = read_hdf5(location2 + "hdf5/")

fontsize = 18
fontsizeBig = 21

# some global stuff to use in '_C_diff_two_files'
AH_center1 = AH_centers1[0][17:20]
AH_center2 = AH_centers2[0][17:20]

variable_comparing = "Cphys"
# variable_comparing = "C"

def _C_diff_two_files(field, data):
    Cphys1 = data[variable_comparing]

    # find the same points on the 2nd file
    unit = file2.length_unit.to('code_length')
    edge = file2.domain_right_edge
    minmax = lambda v, i : np.maximum(np.minimum(v, edge[i]), 0)

    points = np.transpose([
        minmax( data["index", "x"] - unit * AH_center1[0] + unit * AH_center2[0] , 0) ,
        minmax( data["index", "y"] - unit * AH_center1[1] + unit * AH_center2[1] , 1) ,
        minmax( data["index", "z"] - unit * AH_center1[2] + unit * AH_center2[2] , 2) ,
    ])
    print(points.shape)
    # YT can pass the whole domain when defining the field, or subsets
    if hasattr(points[0][0], "__len__"):
        Cphys2 = [[[file2.r[point][variable_comparing][0] for point in pointy] for pointy in pointz] for pointz in points]
    else:
        Cphys2 = [file2.r[point][variable_comparing][0] for point in points]
        # print([data["index", "x"][0], data["index", "y"][0], data["index", "z"][0]])
        # print(Cphys1[0])
        # print(points[0])
        # print(Cphys2[0])

    return np.abs(Cphys1 - Cphys2)

allAverageCminusCphys = []
allAverageCminusCphysWithoutAH = []
for i in range(0, min(len(ds1), len(ds2)), jump):
    file1 = ds1[i]
    file2 = ds2[i]
    time = file1.current_time # assume same!
    
    file1.add_field("C_diff_two_files", _C_diff_two_files, units = "", sampling_type = "cell")
    
    AH_center1 = AH_centers1[AH_centers1[:,0] == time][0][17:20]
    AH_center2 = AH_centers2[AH_centers2[:,0] == time][0][17:20]
    print("AH_centers", AH_center1, AH_center2)

    if z_symmetry:
        AH_center1[2] = 0 # force numeric 0
        AH_center2[2] = 0 # force numeric 0

    currentAverageCminusCphys = []
    currentAverageCminusCphysWithoutAH = []
    for p in np.arange(width_min_deviation, width_max_deviation + 1, 1):
        # compute integral differences between C and Cphys
        power = np.power(2., p)
        ray = file1.r[[AH_center1[0] + 0, AH_center1[1] + 0, AH_center1[2] + 0]:[AH_center1[0] + 0, AH_center1[1] + 0, AH_center1[2] + width*power]]
        averageCminusCphys = ray.mean("C_diff_two_files", weight="cell_volume")
        averageC = abs(ray.mean("C", weight="cell_volume"))
        currentAverageCminusCphys.append([width*power, averageCminusCphys / averageC * 100])

        # now without AH
        ray = file1.r[[AH_center1[0] + 0, AH_center1[1] + 0, AH_center1[2] + AH_radius_vs_time1(time)]:[AH_center1[0] + 0, AH_center1[1] + 0, AH_center1[2] + width*power]]
        averageCminusCphys = ray.mean("C_diff_two_files", weight="cell_volume")
        averageC = abs(ray.mean("C", weight="cell_volume"))
        currentAverageCminusCphysWithoutAH.append([width*power, averageCminusCphys / averageC * 100])

        print(averageCminusCphys, averageC, currentAverageCminusCphysWithoutAH[-1])

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

np.save(f'data_C_slice_yt_integral_compare_two_runs_{variable_comparing}.npy', [times_noAH, rsAndValues_noAH, rsAndValues_withAH])
[times_noAH, rsAndValues_noAH, rsAndValues_withAH] = np.load(f'data_C_slice_yt_integral_compare_two_runs_{variable_comparing}.npy', allow_pickle=True)

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

plt.savefig(f"CvsCphys_integral_compare_two_runs_{variable_comparing}", bbox_inches = 'tight')
plt.close()
