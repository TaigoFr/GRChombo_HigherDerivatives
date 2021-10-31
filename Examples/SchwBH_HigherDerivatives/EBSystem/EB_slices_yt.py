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

minE11 = 10**10
maxE11 = 0
minE11abs = 10**10
maxE11abs = 0
for i in range(0, len(ds), jump):
    file = ds[i]
    ray = file.r[[0, 0, rmin]:[0, 0, rmax]]
    minE11 = min(minE11, min(np.min(ray["E11"]), np.min(ray["Ephys11"])))
    maxE11 = max(maxE11, max(np.max(ray["E11"]), np.max(ray["Ephys11"])))
    minE11abs = min(minE11abs, min(np.min(np.abs(ray["E11"])), np.min(np.abs(ray["Ephys11"]))))
    maxE11abs = max(maxE11abs, max(np.max(np.abs(ray["E11"])), np.max(np.abs(ray["Ephys11"]))))

print("Min/Max E11/Ephys11 = (%f,%f)" % (minE11, maxE11))
print("Min/Max Abs E11/Ephys11 = (%f,%f)" % (minE11abs, maxE11abs))

def do_plot(use_log):
    name_end = "_log" if use_log else ""

    for i in range(0, len(ds), jump):
        file = ds[i]
        ray = file.r[[0, 0, rmin]:[0, 0, rmax]]

        srt = np.argsort(ray["radius"]) # ray does not have elements ordered

        fig, ax = plt.subplots(figsize=(12, 9))
        ax.plot(np.array(ray["index", "z"][srt]), np.abs(np.array(ray["E11"][srt])),     label=r"$|E11|$")
        ax.plot(np.array(ray["index", "z"][srt]), np.abs(np.array(ray["Ephys11"][srt])), label=r"$|E11_{phys}|$")
        ax.legend(fontsize=fontsize)

        # if use_log:
        #     ax.set_ylim([minE11abs * 0.9 , maxE11abs * 1.1])
        # else:
        #     ax.set_ylim([minE11 * (1.1 if minE11 < 0 else 0.9), maxE11 * (0.9 if maxE11 < 0 else 1.1)])
        ax.set_ylim([0 , maxE11abs * 1.05])

        plt.xlabel("Radius (M)", fontsize=fontsize)

        if use_log:
            plt.yscale("log")

        ax.tick_params(axis='both', which='major', labelsize=fontsize)

        time = file.current_time
        fig.suptitle(r'$E11$ vs $E11_{phys}$ ( t = %.2fM )' % time, fontsize=fontsizeBig, position=(0.5,0.93))

        plt.savefig(("E11vsEphys11_%04d" % i) + name_end + ".png", bbox_inches = 'tight')
        plt.close()

    print ("Making a movie...")
    os.system('ffmpeg -f image2 -framerate 2 -i E11vsEphys11_%04d' + name_end + '.png E11vsEphys11' + name_end + '.avi -y')
    print ("I've finished!")

do_plot(False)
do_plot(True)
