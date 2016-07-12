#!/usr/local/anaconda/bin/python

# ------------------------------------------------------------------------------------
#
#           Code to read in _Network file and add reservoir information 
#
# ------------------------------------------------------------------------------------

# path new _Network file will be saved
save_path = '/raid3/rniemeyr/RBM/two_layer_model/RBM_Yixin/RIPS/model_run/source/VIC_RBM/RBM_processing/Tennessee_8th.Mohseni_v1/'
# RUN file by typing: python_code network_file reservoir_file 

import sys
import csv
import pandas as pd
import os.path

# ------------------------------------------------------------------------------------
#        read in reservoir information, write all nodes with a reservoir
# ------------------------------------------------------------------------------------

# reservoir_file = 'reservoir_to_model.csv'
reservoir_file = sys.argv[2]
with open(reservoir_file, 'rU') as csvfile:
	reader = csv.reader(csvfile, delimiter=',', dialect=csv.excel_tab)
	headers = csvfile.next() # read in headers
	#headers = headers[0]
	headers = headers.split(',')
	headers[len(headers)-1] = headers[len(headers)-1].rstrip('\n')  #strip \n
	reservoir_nodes = []
	for row in reader:
		# print ', '.join(row)
		x = ', '.join(row)
		x = x.split(',')
		reservoir_nodes.append(x)

reservoir_nodes = pd.DataFrame(reservoir_nodes)	
reservoir_nodes.columns = headers

# -------- sequence of numbers -----------
reservoir_tot_nodes = []
reservoir_tot_nodes_index = []
for i in range(0, len(reservoir_nodes['start_node']) ):
	x = range(int(reservoir_nodes['start_node'][i]), int(reservoir_nodes['end_node'][i])+1 )
	reservoir_tot_nodes.extend(x)
	x = repeat(i+1, len(x) )  # want index to start with 1, not 0
	reservoir_tot_nodes_index.extend(x)

reservoir_tot_nodes2 = pd.DataFrame(reservoir_tot_nodes)

# ------------------------------------------------------------------------------------
#        read in reservoir information, write all nodes with a reservoir
# ------------------------------------t ------------------------------------------------

# network_file = "Holston_8th_Network"
network_file = sys.argv[1]
with open(network_file) as f:
    network =  f.read().splitlines() 

# -------- get dimensions of river network ----------    
net_str = network[4]
dim = [int(s) for s in net_str.split() if s.isdigit()]  # nreach, flow_cells, heat_cells, source

# ------------- reform the network file ------------
network2 = network[0:3]

x = save_path + reservoir_file  # insert path for reservoir information file
network2.append(x)

network2.append(network[3])   # start and end date

nresx = len(reservoir_nodes)
x = network[4] + '     ' + str(nresx) + '     TRUE'  # a boolean for if reservoirs are simulated with two-layer models
network2.append(x)
				
# ------------- loop to get the nodes each file --------
for i in range(5, len(network)-1):
	x = network[i]
	if x[0:4] == 'Node':
		x2 = [int(s) for s in x.split() if s.isdigit()]
		
		# if node in NETWORK file has a reservoir
		if x2[0] in reservoir_tot_nodes:
			
			nodesx = int(np.where(reservoir_tot_nodes2 == x2[0])[0])
			nodesx2 = reservoir_tot_nodes_index[nodesx]
			x4 = nodesx2.astype('|S10')
			x4 = x4.rjust(5-len(x4))
			x3 = network[i] + '  ' +  x4
			network2.append(x3)
			
		# if node in NETWORK file does not have reservoir
		else:
			x3 = network[i] + '  ' +  '   0'
			network2.append(x3)
	else:
		network2.append(x)

# ------------- write the file --------
network_file = network_file + "_2"  # REMOVE once you have a network file you like (i.e. overwrite network file)
completeName = os.path.join(save_path, network_file)         

with open (network_file, 'w') as fo:
   for d in network2:
     fo.write(d + '\n')
     
    
