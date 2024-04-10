#!/usr/bin/python3

import shasta
import argparse

parser = argparse.ArgumentParser(description=
    'Load a mode3::AssemblyGraph representing a connected component of the primary graph and assemble it.')
parser.add_argument('fileName', type=str, 
    help='Filename ')
parser.add_argument('threadCount', type=int, 
    help='Number of threads.')
parser.add_argument(
    "--debug",
    dest="debug",
    action="store_true",
)  
        
arguments = parser.parse_args() 
if arguments.threadCount <= 0:
    raise Exception("Numbers of threads must be positive.")    

a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.accessMarkerGraphReverseComplementEdge()
a.accessMarkerGraphConsensus()
shasta.openPerformanceLog('Mode3AssembleComponent.log')
a.mode3AssembleComponent(arguments.fileName, arguments.threadCount, arguments.debug)
 
