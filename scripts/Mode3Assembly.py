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
    
"""
parser.add_argument(
    "--recomputePrimaryJourneys",
    dest="recomputePrimaryJourneys",
    action="store_true",
    help="If this boolean flag is specified, the primary journeys "
         "are recomputed. This is only necessary if "
         "--MarkerGraph.minPrimaryEdgeCoverage or "
         "--MarkerGraph.maxPrimaryEdgeCoverage have changed."
) 
"""  
        
parser.add_argument(
    "--debug",
    dest="debug",
    action="store_true",
)    
        
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

# Open a performance log.
shasta.openPerformanceLog('Mode3Assembly.log')

"""
# If requested, recompute the primary journeys.
if arguments.recomputePrimaryJourneys:
    a.flagPrimaryMarkerGraphEdges(
        int(config['MarkerGraph']['minPrimaryEdgeCoverage']), 
        int(config['MarkerGraph']['maxPrimaryEdgeCoverage']), 
        0)
    a.createMarkerGraphPrimaryJourneys(0)
"""

# Run Mode 3 assembly.
# a.mode3Assembly0(arguments.threadCount0, arguments.threadCount1)
a.mode3Assembly1(arguments.threadCount0, arguments.debug)
 
