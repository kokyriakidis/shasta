#!/usr/bin/python3

import shasta
import argparse
import GetConfig

# Read the config file.
config = GetConfig.getConfig()

# Parse the command line arguments.
parser = argparse.ArgumentParser(description=
    'Rerun Mode 3 assembly.')
parser.add_argument('threadCount0', type=int, 
    help='Number of threads for high level parallelization')
parser.add_argument('threadCount1', type=int, 
    help='Number of threads for low level parallelization')
        
arguments = parser.parse_args() 
if arguments.threadCount0 <= 0 or arguments.threadCount1 <=0:
    raise Exception("Numbers of threads must be positive.")    


# Create the Assembler object and access what we need.
a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges(True)
a.accessMarkerGraphReverseComplementEdge()
a.accessMarkerGraphConsensus()
a.accessMarkerGraphPrimaryJourneys()

# Redo Mode 3 assembly.
shasta.openPerformanceLog('Mode3Reassembly.log')
a.flagPrimaryMarkerGraphEdges(
    int(config['MarkerGraph']['minPrimaryEdgeCoverage']), 
    int(config['MarkerGraph']['maxPrimaryEdgeCoverage']), 
    0)
a.createMarkerGraphPrimaryJourneys(0)
a.findMode3bPaths(arguments.threadCount0, arguments.threadCount1)
 
