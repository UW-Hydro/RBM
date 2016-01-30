#!/usr/local/anaconda/bin/python

# This script generates:
#   - A routing station file (for the next step to generate flow and energy file)
#   - RBM control file (with running period and flow and energy file missing and to be subsitute)
# Note:
#   - For different basin or different Mohseni parameters, this script needs to be rerun

import numpy as np
import sys
import subprocess
import my_functions

cfg = my_functions.read_config(sys.argv[1])  # Read config file

#==========================================================#
# Read config file parameters
#==========================================================#
# [INPUT]
# Flow direction file, arcGIS ascii format
# 1-8 for direction, 9 for basin outlet (only one outlet), -1 for inactive grid cells
flowdir_asc = cfg['INPUT']['flowdir_asc']
# Template control file for RBM parameter preparation
# The following options will be filled in:
#   - <OUTPUT_DIR> # Directory for output routing station files
#   - <BASIN_CODE> # basin_code
#   - <TOPOLOGY_FILE>
#   - <NETWORK_FILE>  # RBM control file to be generated
#   - <MOHSENI_FILE>
# The following options will NOT be filled in (but will be filled in directly in the RBM control file when preparing flow and energy input files for RBM):
#   - <START_DATE>
#   - <END_DATE>
#   - <OUTPUT_FLOW_FILE>
#   - <OUTPUT_ENERGY_FILE>
control_template = cfg['INPUT']['control_template']

# [TOOLS]
# Perl script for building topology file
topo_pl = cfg['TOOLS']['topo_pl']
# Perl script for generating RBM control file & routing station file
build_input_pl = cfg['TOOLS']['build_input_pl']

# [MOHSENI]
# Mohseni parameters (currently spatially constants
# alpha is actually 'alpha-mu'
mohseni_param = {}
mohseni_param['alpha'] = cfg['MOHSENI']['alpha']
mohseni_param['beta'] = cfg['MOHSENI']['beta']
mohseni_param['gamma'] = cfg['MOHSENI']['gamma']
mohseni_param['mu'] = cfg['MOHSENI']['mu']
mohseni_param['timelag'] = cfg['MOHSENI']['timelag']

# [OUTPUT]
# Output directory for all temporary files in the process
output_tmp_dir = cfg['OUTPUT']['output_tmp_dir'] 
# Basin code, will be used as basename for topology and Mohseni parameter files
basin_code = cfg['OUTPUT']['basin_code']

#==========================================================#
# Generate topology file
#==========================================================#
subprocess.call('perl {} {} {}/{}.Topology'\
                    .format(topo_pl, flowdir_asc, output_tmp_dir, basin_code), \
                shell=True)

#==========================================================#
# Generate Mohseni parameter files
# (Currntly, Mohseni parameters are set to spatially constants)
#==========================================================#
#=== Load flow direction file header ===#
# Read header
header = ''
f = open(flowdir_asc, 'r')
for i in range(6):
    line = f.readline()
    header = header + line
f.close()
# Extract number of rows and columns
ncols = int(header.split()[1])
nrows = int(header.split()[3])

#=== Write Mohseni parameter files ===#
for param in mohseni_param.keys():
    # Create Mohseni parameter array
    param_array = np.ones([nrows, ncols]) * mohseni_param[param]
    f = open('{}/{}.Mohseni.{}'.format(output_tmp_dir, basin_code, param), 'w')
    f.write(header)
    np.savetxt(f, param_array, fmt='%.2f')
    f.close()

#==========================================================#
# Prepare control file (for RBM input preparation)
#==========================================================#
#=== Read in template file ===#
f = open(control_template, 'r')
content = f.read()
f.close()
#=== Replace options ===#
content = content.replace('<OUTPUT_DIR>', \
                          '{}'.format(output_tmp_dir))
content = content.replace('<BASIN_CODE>', \
                          '{}'.format(basin_code))
content = content.replace('<TOPOLOGY_FILE>', \
                          '{}/{}.Topology'.format(output_tmp_dir, basin_code))
content = content.replace('<NETWORK_FILE>', \
                          '{}/{}_Network'.format(output_tmp_dir, basin_code))
content = content.replace('<MOHSENI_FILE>', \
                          '{}/{}.Mohseni'.format(output_tmp_dir, basin_code))
#=== Write new control file ===#
f = open('{}/{}.RBM_prep.Control'.format(output_tmp_dir, basin_code), 'w')
f.write(content)
f.close()

#==========================================================#
# Generate RBM control file & routing station file
#==========================================================#
subprocess.call("{} {}/{}.RBM_prep"\
                    .format(build_input_pl, output_tmp_dir, basin_code), \
                shell=True)  # ".Control" is appended 



