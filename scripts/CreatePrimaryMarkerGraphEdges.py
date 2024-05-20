#!/usr/bin/python3

import shasta
import GetConfig

config = GetConfig.getConfig()


a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessDisjointSetsHistogram()
a.createPrimaryMarkerGraphEdges(
    minPrimaryCoverage = int(config['Assembly']['mode3.minPrimaryCoverage']),
    maxPrimaryCoverage = int(config['Assembly']['mode3.maxPrimaryCoverage']))

