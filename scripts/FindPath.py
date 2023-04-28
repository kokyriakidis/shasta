#!/usr/bin/python3

import shasta
import argparse

parser = argparse.ArgumentParser(description=
    'Find and assemble a path in the complete marker graph.')
    
parser.add_argument('startEdgeId', type=int,
    help='The id of the marker graph edge to start from.')
parser.add_argument('direction', type=int, choices=range(2),
    help='The direction to move in while creating the path: 0=forward, 1=backward.')
    
arguments = parser.parse_args()   

a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.findCompleteMarkerGraphPath(arguments.startEdgeId, arguments.direction)
 
