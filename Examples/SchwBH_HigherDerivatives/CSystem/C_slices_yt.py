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

rmin = 0.
rmax = 5.
location = './HigherDerivatives'
is_corner = True
jump = 1 # plot every 'jump' files

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

ds, center = read_hdf5(location, is_corner)

assert((center == [0,0,0]).all()) # fix ray domain below if not

# ds = ds[0:5]

fontsize = 18
fontsizeBig = 21

minC = 10**10
maxC = 0
minCabs = 10**10
maxCabs = 0
for i in range(0, len(ds), jump):
    file = ds[i]
    ray = file.r[[0, 0, rmin]:[0, 0, rmax]]
    minC = min(minC, min(np.min(ray["C"]), np.min(ray["Cphys"])))
    maxC = max(maxC, max(np.max(ray["C"]), np.max(ray["Cphys"])))
    minCabs = min(minCabs, min(np.min(np.abs(ray["C"])), np.min(np.abs(ray["Cphys"]))))
    maxCabs = max(maxCabs, max(np.max(np.abs(ray["C"])), np.max(np.abs(ray["Cphys"]))))

print("Min/Max C/Cphys = (%f,%f)" % (minC, maxC))
print("Min/Max Abs C/Cphys = (%f,%f)" % (minCabs, maxCabs))

def do_plot(use_log):
    name_end = "_log" if use_log else ""

    for i in range(0, len(ds), jump):
        file = ds[i]
        ray = file.r[[0, 0, rmin]:[0, 0, rmax]]

        srt = np.argsort(ray["radius"]) # ray does not have elements ordered

        fig, ax = plt.subplots(figsize=(12, 9))
        ax.plot(np.array(ray["index", "z"][srt]), np.array(ray["C"][srt]),     label=r"$C$")
        ax.plot(np.array(ray["index", "z"][srt]), np.array(ray["Cphys"][srt]), label=r"$C_{phys}$")
        ax.legend(fontsize=fontsize)

        if use_log:
            ax.set_ylim([minCabs * 0.9 , maxCabs * 1.1])
        else:
            ax.set_ylim([minC * (1.1 if minC < 0 else 0.9), maxC * (0.9 if maxC < 0 else 1.1)])

        plt.xlabel("Radius (M)", fontsize=fontsize)

        if use_log:
            plt.yscale("log")

        ax.tick_params(axis='both', which='major', labelsize=fontsize)

        time = file.current_time
        fig.suptitle(r'$C$ vs $C_{phys}$ ( t = %.2fM )' % time, fontsize=fontsizeBig, position=(0.5,0.93))

        plt.savefig(("CvsCphys_%04d" % i) + name_end + ".png", bbox_inches = 'tight')
        plt.close()

    print ("Making a movie...")
    os.system('ffmpeg -f image2 -framerate 1 -i CvsCphys_%04d' + name_end + '.png CvsCphys' + name_end + '.avi -y')
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
