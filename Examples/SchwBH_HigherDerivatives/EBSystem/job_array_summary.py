# this file is intended to generate a summary of the result of an array job from the 'job_array_setup.py'

# use as 'python job_array_summary.py total'

import sys # to get command line arguments
import json
from os import path

if len(sys.argv)>1:
    total = int(sys.argv[1])

    folders = ['search_%d' % i for i in range(1,total+1)]

    # assume params to vary is common for all folders
    params_to_vary = []
    specific_params = open(folders[0] + "/specific_params.txt")
    lines = specific_params.readlines()
    specific_params.close()

    for line in lines:
        if line[0] != "#" and len(line) > 3:
            params_to_vary.append(line.split()[0])

    # now read all values
    params = []
    success = []
    final_times = []
    ah_areas = []
    ah_masses = []
    for i in range(total):
        params.append([])
        for line in lines:
            if line[0] != "#" and len(line) > 3:
                params[i].append(float(line.split()[-1]))

        pout = open(folders[i] + "/pout/pout.0")
        pout_lines = pout.readlines()
        pout.close()

        for l in range(len(pout_lines)):
            line = pout_lines[len(pout_lines) - 1 - l]
            if 'GRAMRLevel' in line:
                final_times.append(float(line.split()[5]))
                break

        success.append(pout_lines[-1] == 'GRChombo finished.\n')

        if path.exists(folders[i] + "/data/stats_AH1.dat"):
            AH_stats = open(folders[i] + "/data/stats_AH1.dat")
            stats_lines = AH_stats.readlines()
            AH_stats.close()
            lastline = stats_lines[-1].split()
            ah_areas.append(float(lastline[2]))
            ah_masses.append(float(lastline[3]))
            # TODO write E_diff, B_diff, etc

    print(params_to_vary)
    print(params)
    print(success)
    print(final_times)
    if len(ah_areas) > 0:
        print(ah_areas)
        print(ah_masses)

    output = open("job_array_summary.txt", 'w')

    output.write("# params_to_vary:\n")
    json.dump(params_to_vary, output)
    output.write("\n\n# params:\n")
    json.dump(params, output)
    output.write("\n\n# success:\n")
    json.dump(success, output)
    output.write("\n\n# final_times:\n")
    json.dump(final_times, output)

    if len(ah_areas) > 0:
        output.write("\n\n# ah areas:\n")
        json.dump(ah_areas, output)
        output.write("\n\n# ah masses:\n")
        json.dump(ah_masses, output)

    output.close()

else:
    print("Run as 'python job_array_summary.py total'")
