#!/usr/bin/python3

import shasta
import argparse

# Parse the arguments.
parser = argparse.ArgumentParser(description=
    'Assembly graph detangling using path following.')    
parser.add_argument('stage', type=str)
parser.add_argument('component', type=int)
arguments = parser.parse_args()
stage = arguments.stage
component = arguments.component

# Get Shasta options.
options = shasta.AssemblerOptions('shasta.conf')

# Create the Assembler and access what we need.
assembler = shasta.Assembler()
assembler.accessMarkers()
assembler.accessMode3Assembler()

# Create the mode3::AssemblyGraph for this assembly stage and component.
assemblyGraph = shasta.Mode3AssemblyGraph(stage, component, assembler, options.assemblyOptions.mode3Options)

# Detangle with path following.
assemblyGraph.detangle2()
