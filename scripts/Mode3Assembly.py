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



# Create the Assembler object and access what we need.
options = shasta.AssemblerOptions('shasta.conf')
a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges(True)
a.accessMarkerGraphReverseComplementEdge()
a.accessMarkerGraphConsensus()
a.accessDisjointSetsHistogram()

# Open a performance log.
shasta.openPerformanceLog('Mode3Assembly.log')

# Flag primary marker graph edges.
a.flagPrimaryMarkerGraphEdges(
    int(config['MarkerGraph']['minPrimaryEdgeCoverage']), 
    int(config['MarkerGraph']['maxPrimaryEdgeCoverage']), 
    0)

# Run Mode 3 assembly.
a.mode3Assembly(0, options.assemblyOptions.mode3Options, arguments.debug)
 
