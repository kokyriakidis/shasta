#!/usr/bin/python3

import shasta
import argparse

parser = argparse.ArgumentParser(description=
    'Load a mode3::AssemblyGraph representing a connected component of the primary graph and assemble it.')
    
parser.add_argument('component', type=int, 
    help='The connected component to assemble.')
    
parser.add_argument(
    "--no-assemble-sequence",
    dest="dontAssembleSequence",
    action="store_true",
)  
        
parser.add_argument(
    "--debug",
    dest="debug",
    action="store_true",
)  
        
arguments = parser.parse_args() 



options = shasta.AssemblerOptions('shasta.conf')
a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphEdgeMarkerIntervals()
a.accessMarkerGraphConsensus()
shasta.openPerformanceLog('Mode3AssembleComponent.log')
fileName = 'AssemblyGraph-' + str(arguments.component) + '.data'
a.mode3AssembleComponent(fileName, 0, 
	options.assemblyOptions.mode3Options, not arguments.dontAssembleSequence, arguments.debug)
 
