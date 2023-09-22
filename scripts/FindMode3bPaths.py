#!/usr/bin/python3

import shasta
import argparse

parser = argparse.ArgumentParser(description=
    'Find and assemble paths in the complete marker graph.')
parser.add_argument('threadCount0', type=int, 
    help='Number of threads for high level parallelization')
parser.add_argument('threadCount1', type=int, 
    help='Number of threads for low level parallelization')
        
arguments = parser.parse_args() 
if arguments.threadCount0 <= 0 or arguments.threadCount1 <=0:
    raise Exception("Numbers of threads must be positive.")    

a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.accessMarkerGraphConsensus()
shasta.openPerformanceLog('FindMode3Paths.log')
a.findMode3bPaths(arguments.threadCount0, arguments.threadCount1)
 
