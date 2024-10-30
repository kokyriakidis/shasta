#!/usr/bin/python3

import shasta

shasta.openPerformanceLog('CountKmers.log')

a = shasta.Assembler()
a.accessMarkers()
a.countKmers()

