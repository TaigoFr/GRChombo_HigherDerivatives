# this file is intended to plot the summary of the result of an array job from the 'job_array_summary.py'

# use as 'python job_array_plot.py summary_file.txt'

import sys # to get command line arguments
import json
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

if len(sys.argv)>1:
    summary_filename = sys.argv[1]

    summary_file = open(summary_filename, 'r')

    lines = summary_file.readlines()

    content = []
    for line in lines:
        if line[0]!="#" and len(line) > 3:
            content.append(json.loads(line))

    params_to_vary = content[0]
    params = np.array(content[1])
    success = np.array(content[2])
    final_times = np.array(content[3])

    # PLOT OF FINAL TIME

    y0 = params[0][1]
    ny = 1
    for i in range(1,len(params)):
        if params[i][1] == y0:
            ny += 1

    nx = len(params)//ny

    xs = np.log2(params[:,0])
    ys = np.log2(params[:,1])

    xs_split = np.array(np.split(xs, ny))
    ys_split = np.array(np.split(ys, nx))
    final_times_split = np.array(np.split(final_times, ny))

    fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    ax = Axes3D(fig)
    ax.plot_surface(xs_split, ys_split, final_times_split)

    ax.set_xlabel("log2(" + params_to_vary[0] + ")")
    ax.set_ylabel("log2(" + params_to_vary[1] + ")")
    ax.set_zlabel("Duration of run (M)")

    plt.draw()
    out_name = 'job_duration.png'
    plt.savefig(out_name, bbox_inches = 'tight')

    plt.show()
    plt.close()

    # PLOT OF AH FINAL MASS
    if len(content) >= 6:
        ah_masses = np.array(content[5])

        valid_x = xs[success]
        valid_y = ys[success]
        valid_ah_masses = ah_masses[success]

        if len(valid_ah_masses) > 0:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(valid_x, valid_y, valid_ah_masses)

            ax.set_xlabel("log2(" + params_to_vary[0] + ")")
            ax.set_ylabel("log2(" + params_to_vary[1] + ")")
            ax.set_zlabel("AH final mass (M)")

            fig.tight_layout()
            plt.draw()
            out_name = 'AH_final_mass.png'
            plt.savefig(out_name, bbox_inches = 'tight')

            plt.show()
            plt.close()

else:
    print("Run as 'python job_array_plot.py summary_file.txt'")