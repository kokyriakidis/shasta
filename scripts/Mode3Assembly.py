#!/usr/bin/python3

import shasta
import argparse
import GetConfig

# Read the config file.
config = GetConfig.getConfig()

# Parse the command line arguments.
parser = argparse.ArgumentParser(description=
    'Run Mode 3 assembly starting from the marker graph.')
            
parser.add_argument(
    "--debug",
    dest="debug",
    action="store_true",
)    
        
arguments = parser.parse_args() 

# Open a performance log.
shasta.openPerformanceLog('Mode3Assembly.log')

# Create the Assembler object and access what we need.
options = shasta.AssemblerOptions('shasta.conf')
a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdgeMarkerIntervals()
a.accessMarkerGraphConsensus()
a.accessDisjointSetsHistogram()

# Run Mode 3 assembly.
a.mode3Assembly(0, options.assemblyOptions.mode3Options, arguments.debug)
 
