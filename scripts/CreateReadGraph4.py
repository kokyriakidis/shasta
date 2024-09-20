#!/usr/bin/python3

import os
import shasta

# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessMarkers()
a.accessAlignmentDataReadWrite()
a.accessCompressedAlignments()

a.createReadGraph4()


