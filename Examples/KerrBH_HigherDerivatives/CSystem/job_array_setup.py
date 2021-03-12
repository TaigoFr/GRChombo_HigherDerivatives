# this file is intended to generate an auxiliary text file with the specific parameters you want to vary.
# (e.g. tau and sigma)

# it will create a file 'specific_params.txt' for the index 'i' of an array-job
# Run as: 'python job_array_setup.py index total', where index is the index in the array, of a total of 'total' jobs
# Then, run GRChombo as 'MainExample.ex params.txt FILE=specific_params.txt'


import sys # to get command line arguments
import math # pow

params_to_vary = ["tau" , "c_sigma"]
params_max     = [1. , 1.]
params_ratio   = [0.5 , 0.5] # search range is [1, 0.5, 0.25, 0.125, ...]

if len(sys.argv)>2:
    index = int(sys.argv[1]) # indices in job arrays start in 1
    total = int(sys.argv[2])

    assert(index > 0 and index <= total)

    num_params = len(params_to_vary)

    # round to nearest perfect "square"
    # add 1e-7 just to make sure 'int(' doesn't suffer from float-point precision errors
    total_per_params = int(math.pow(total, 1. / num_params) + 1e-7)
    actual_total = math.pow(total_per_params, num_params)

    if actual_total < total:
        print("Total %d is not a perfect square." % total)
    print("Using %d jobs for each of the %d parameters, for a total of %d." % (total_per_params, num_params, actual_total))

    if index > actual_total:
        print('Not running index %d bigger than the actual total to run %d' % (index, actual_total))
        exit()

    sub_total_index = index - 1
    indices = [0]*num_params
    params = [0]*num_params
    file = open('specific_params.txt', 'w')
    file.write("# file generated automatically for index %d / %d\n\n" % (index, actual_total))
    for p in range(num_params):
        indices[p] = sub_total_index % total_per_params
        sub_total_index = int(sub_total_index / total_per_params)
        params[p] = params_max[p] * math.pow(params_ratio[p], indices[p])
        file.write("# Parameter %d:\n" % p)
        file.write("# Index = %d\n" % indices[p])
        file.write("# Using param_max = %f and param_ratio = %f\n" % (params_max[p], params_ratio[p]))
        file.write(params_to_vary[p] + " = " + str(params[p]) + "\n\n")

    print("Indices    for this job:", indices)
    print("Parameters for this job:", params)
    file.close()

else:
    print("Run as 'python job_array_setup.py index total'")