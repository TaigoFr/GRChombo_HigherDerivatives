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

width = 2.
location = '../' # base folder with /data and /hdf5 subfolders
jump = 1 # plot every 'jump' files
z_symmetry = True

#################################################################

yt.enable_parallelism()



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
            break #time1 = 1e10 # big number
        if i3 < len(stats_AH3):
            row3 = stats_AH3[i3]
            time3 = row3[0]
        else:
            time3 =1e10 #break # finished

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

print(times, radii)

AH_radius_vs_time = make_interp_spline(times, radii, k=3)



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

    

ds, center = read_hdf5(location + "hdf5/")
# ds = ds[0:5]

punctures = np.loadtxt(location + "data/punctures.dat")

fontsize = 18
fontsizeBig = 21

minC = 10**10
maxC = 0
minCabs = 10**10
maxCabs = 0
for i in range(0, len(ds), jump):
    file = ds[i]
    time = file.current_time
    puncture = punctures[punctures[:,0] == time][0][1:4]
    if z_symmetry:
        puncture[2] = 0 # force numeric 0
    radius = np.linalg.norm(puncture-center)
    point_min = center + (puncture-center) * (radius - width) / radius
    point_max = center + (puncture-center) * (radius + width) / radius
    ray = file.r[point_min:point_max]
    minC = min(minC, min(np.min(ray["C"]), np.min(ray["Cphys"])))
    maxC = max(maxC, max(np.max(ray["C"]), np.max(ray["Cphys"])))
    minCabs = min(minCabs, min(np.min(np.abs(ray["C"])), np.min(np.abs(ray["Cphys"]))))
    maxCabs = max(maxCabs, max(np.max(np.abs(ray["C"])), np.max(np.abs(ray["Cphys"]))))

print("Min/Max C/Cphys = (%f,%f)" % (minC, maxC))
print("Min/Max Abs C/Cphys = (%f,%f)" % (minCabs, maxCabs))

def do_plot(use_log):
    minCabsLocal = minCabs
    if use_log:
        minCabsLocal = max(1e-4, minCabs) # might be 0

    name_start = "log_" if use_log else ""

    for i in range(0, len(ds), jump):
        file = ds[i]
        time = file.current_time
        puncture = punctures[punctures[:,0] == time][0][1:4]
        if z_symmetry:
            puncture[2] = 0 # force numeric 0
        radius = np.linalg.norm(puncture-center)
        point_min = center + (puncture-center) * (radius - width) / radius
        point_max = center + (puncture-center) * (radius + width) / radius
        ray = file.r[point_min:point_max]

        srt = np.argsort(ray["radius"]) # ray does not have elements ordered

        unit = file.length_unit.to('code_length')
        points = np.transpose([ray["index", "x"][srt] - unit * center[0], ray["index", "y"][srt] - unit * center[1], ray["index", "z"][srt] - unit * center[2]])
        radii = np.array([np.linalg.norm(point) - radius for point in points])

        fig, ax = plt.subplots(figsize=(12, 9))
        ax.plot(radii, np.array(ray["C"][srt]),     label=r"$C$")
        ax.plot(radii, np.array(ray["Cphys"][srt]), label=r"$C_{phys}$")

        ax.plot(radii, 10000*(radii - AH_radius_vs_time(time)),'k')
        ax.plot(radii, 10000*(radii + AH_radius_vs_time(time)),'k')

        ax.legend(fontsize=fontsize)

        if use_log:
            ax.set_ylim([minCabsLocal * 0.9 , maxCabs * 1.1])
        else:
            ax.set_ylim([minC * (1.1 if minC < 0 else 0.9), maxC * (0.9 if maxC < 0 else 1.1)])

        plt.xlabel("Radius (M)", fontsize=fontsize)

        if use_log:
            plt.yscale("log")

        ax.tick_params(axis='both', which='major', labelsize=fontsize)

        fig.suptitle(r'$C$ vs $C_{phys}$ ( t = %.2fM )' % time, fontsize=fontsizeBig, position=(0.5,0.93))

        plt.savefig(name_start + ("CvsCphys_%04d" % i) + ".png", bbox_inches = 'tight')
        plt.close()

    print ("Making a movie...")
    os.system('ffmpeg -f image2 -framerate 1 -i ' + name_start + 'CvsCphys_%04d.png ' + name_start + 'CvsCphys.avi -y')
    print ("I've finished!")

do_plot(False)
do_plot(True)


# old attempts with ProfilePlot and LinePlot

# for i in range(0, len(ds), jump):
#     file = ds[i]
#     sphere = file.sphere(center, 10)

#     # profiles = [yt.create_profile(sphere, "radius", fields="C", weight_field=None),
#                 # yt.create_profile(sphere, "radius", fields=["C","Cphys"], weight_field=None)]
#     profiles = []
#     profiles.append(yt.create_profile(sphere, "radius", ["C","Cphys"], weight_field=None))
#     profiles.append(yt.create_profile(sphere, "radius", ["Cphys","C"], weight_field=None))
#     labels = ["C", "Cphys"]
#     plot = yt.ProfilePlot.from_profiles(profiles, labels=labels)
#     # plot = yt.ProfilePlot(sphere, "radius", [("chombo", "C"), ("chombo", "Cphys")], weight_field=None, label=[r"$C$", r"$C_{phys}$"], x_log=False)
#     # plot.set_log("radius", False)
#     # plot.set_log("C", False)
#     # plot.set_xlabel("Radius")
#     plot.save("CvsCphys_%04d_log.png" % i)


# for i in range(0, len(ds), jump):
# # for file in ds:
#     file = ds[i]
#     plot = yt.LinePlot(
#         file,
#         ["C", "Cphys"],
#         [0, 0, rmin],
#         [0, 0, rmax],
#         100,
#         field_labels={"C": r"$C$", "Cphys": r"$C_{phys}$"}
#     )
#     plot.annotate_legend("C")
#     plot.set_xlabel("Radius")
#     plot.set_ylabel("")
#     plot.save("CvsCphys_%04d_log.png" % i, mpl_kwargs={"xbounds" : (0.001, 0.1)})
