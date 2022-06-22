# on a cluster with paralellism, run as
# 'mpiexec -n 8 -hosts=i01r01c01s01 python C_slices.py'
# host is the node running the script, 8 is the number of cores
# I realized making this '8' bigger does not help that much

import numpy as np
import yt
import matplotlib.pyplot as plt
import os

#################################################################
# USER DATA

width = 100.
location = '../' # base folder with /data and /hdf5 subfolders
jump = 10 # plot every 'jump' files
z_symmetry = True

#################################################################

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

    

ds, center = read_hdf5(location + "hdf5/")
# ds = ds[0:5]

punctures = np.loadtxt(location + "data/punctures.dat")

fontsize = 18
fontsizeBig = 21

minP = 10**10
maxP = 0
minPabs = 10**10
maxPabs = 0
sigma_asymp =1.0
sigma = 0.05
sigma_length = 0.92
sigma_width = 0.025 
for i in range(0, len(ds), jump):
    file = ds[i]
    time = file.current_time
    puncture = punctures[punctures[:,0] == time][0][1:4]
    if z_symmetry:
        puncture[2] = 0 # force numeric 0
    radius = np.linalg.norm(puncture-center)
    point_min = puncture * (radius - width) / radius
    point_max = puncture * (radius + width) / radius
    ray = file.r[point_min:point_max]
    chi = ray["chi"] 
    sigma_p = (sigma_asymp-sigma)/(1.0 +np.exp(-(chi/sigma_length -1.0)/sigma_width)) +sigma
    minP = min(minP,np.min(sigma_p))
    maxP = max(maxP,np.max(sigma_p))
    minPabs = min(minPabs, np.min(np.abs(sigma_p)))
    maxPabs = max(maxPabs, np.max(np.abs(sigma_p)))
    

def do_plot(use_log):
    minPabsLocal = minPabs

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
        chi = ray["chi"]
        sigma_p = (sigma_asymp-sigma)/(1.0 +np.exp(-(chi/sigma_length -1.0)/sigma_width)) + sigma        
        ax.plot(radii, np.array(sigma_p[srt]), label=r"$sigma$")
        ax.legend(fontsize=fontsize)

        ax.set_ylim([minP * (1.1 if minP < 0 else 0.9), maxP * (0.9 if maxP < 0 else 1.1)])

        plt.xlabel("Radius (M)", fontsize=fontsize)

        if use_log:
            plt.yscale("log")

        ax.tick_params(axis='both', which='major', labelsize=fontsize)

        fig.suptitle(r'$sigma$ ( t = %.2fM )' % time, fontsize=fontsizeBig, position=(0.5,0.93))

        plt.savefig(name_start + ("sigma_%04d" % i) + ".png", bbox_inches = 'tight')
        plt.close()

    print ("Making a movie...")
    os.system('ffmpeg -f image2 -framerate 1 -i ' + name_start + 'sigma_%04d.png ' + name_start + 'sigma.avi -y')
    print ("I've finished!")

do_plot(False)

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
