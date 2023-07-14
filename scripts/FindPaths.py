#!/usr/bin/python3

import shasta
import argparse

parser = argparse.ArgumentParser(description=
    'Find and assemble paths in the complete marker graph.')
        
arguments = parser.parse_args()   

a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.accessMarkerGraphConsensus()
a.findCompleteMarkerGraphPaths()
 
