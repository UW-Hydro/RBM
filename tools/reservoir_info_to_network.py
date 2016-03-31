#!/usr/local/anaconda/bin/python

# ------------------------------------------------------------------------------------
#
#           Code to read in _Network file and add reservoir information
#
# ------------------------------------------------------------------------------------

import sys

# ------------------------------------------------------------------------------------
#        read in reservoir information, write all nodes with a reservoir
# ------------------------------------------------------------------------------------

f = open(sys.argv[1])  # run the script like this: ./reservoir_info_to_network.py reservoir_to_model.csv
csv_f = csv.reader(f)

# ----- loop to read in first and last nodes, sequence across those, and get
#        a long array with all the nodes with reservoirs
reservoir_nodes = []
for row in csv_f:
        x = row[3:5]  # get the first and last node
        x = range(int(x[0]),(int(x[1])+1))
        reservoir_nodes.extend(x)

# ------------------------------------------------------------------------------------
#        read in reservoir information, write all nodes with a reservoir
# ------------------------------------------------------------------------------------

with open("Salmon_0.5_Network") as f:
    network =  f.read().splitlines()

# -------- get dimensions of river network ----------
str = network[4]
dim = [int(s) for s in str.split() if s.isdigit()]  # nreach, flow_cells, heat_cells, source

# ------------- reform the network file ------------
network2 = network[0:4]
x = network[4] + '     TRUE'  # a boolean for if reservoirs are simulated with two-layer models
network2.append(x)

# ------------- loop to get the nodes each file --------
for i in range(5, len(network)-1):
        x = network[i]
        if x[0:4] == 'Node':
                x2 = [int(s) for s in x.split() if s.isdigit()]

                # if node in NETWORK file has a reservoir
                if x2[0] in reservoir_nodes:
                        x3 = network[i] + '     TRUE'
                        network2.append(x3)

                # if node in NETWORK file does not have reservoir
                else:
                        x3 = network[i] + '     FALSE'
                        network2.append(x3)
        else:
                network2.append(x)

# ------------- write the file --------
with open ('myfile', 'w') as fo:
   for d in network2:
     fo.write(d + '\n')



