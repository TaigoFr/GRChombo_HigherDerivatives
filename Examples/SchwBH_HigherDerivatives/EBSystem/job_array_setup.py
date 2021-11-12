# this file is intended to generate an auxiliary text file with the specific parameters you want to vary.
# (e.g. tau and sigma)

# it will create a file 'specific_params.txt' for the index 'i' of an array-job
# Run as: 'python job_array_setup.py index', where index is the index in the array (total set below)
# Then, run GRChombo as 'MainExample.ex params.txt FILE=specific_params.txt'

import sys # to get command line arguments
import math # pow

# eb_version = 2
# use_last_index_raised = 1

float_params_start   = 10.
float_params_ratio   = 1./3.
float_params_iterations = 7
float_params_list = [float_params_start * pow(float_params_ratio, i) for i in range(float_params_iterations)]

params_to_vary   = ["rescale_tau_by_lapse", "rescale_sigma_by_lapse", ["advection_type", "advection_coeff"], "tau" , "eb_sigma"]
params_iterations = [[0,1], [0,1,2], [[0, 0], [1, 0.5], [1, 1.], [2, 0.5], [2, 1.]], float_params_list, float_params_list]

print("List of floats =", float_params_list)
num_params = len(params_to_vary)

total = 1
for p in range(num_params):
    total *= len(params_iterations[p])

print("Total =", total) # 1470 = 2 * 3 * 5 * 7 * 7

if len(sys.argv)>1:
    index = int(sys.argv[1])-1 # indices in job arrays start in 1

    assert(index >= 0 and index < total)

    file = open('specific_params.txt', 'w')
    file.write("# file generated automatically for index %d / %d\n\n" % (index, total))

    indices = [0]*num_params
    params = [0]*num_params
    for p in range(num_params):
        size_of_batch = total / len(params_iterations[p])
        indices[p] = int(index // size_of_batch)
        index = index % size_of_batch
        total = size_of_batch
        params[p] = params_iterations[p][indices[p]]
        file.write("# Parameter %d:\n" % p)
        file.write("# Index = %d\n" % indices[p])
        if isinstance(params_to_vary[p], list):
            for p2 in range(len(params_to_vary[p])):
                file.write(params_to_vary[p][p2] + " = " + str(params[p][p2]) + "\n")
            file.write("\n")
        else:
            file.write(params_to_vary[p] + " = " + str(params[p]) + "\n\n")

    print("Indices    for this job:", indices)
    print("Parameters for this job:", params)
    file.close()

else:
    print("Run as 'python job_array_setup.py index'")
