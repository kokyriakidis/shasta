#!/usr/bin/python3

import shasta
import GetConfig

# Read the config file.
config = GetConfig.getConfig()
        
a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges(True)
a.accessDisjointSetsHistogram()
a.flagPrimaryMarkerGraphEdges(
    int(config['MarkerGraph']['minPrimaryEdgeCoverage']), 
    int(config['MarkerGraph']['maxPrimaryEdgeCoverage']), 
    0)
 
