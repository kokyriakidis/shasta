#!/usr/bin/python3

import shasta
import sys

helpMessage = """
Invoke with two arguments, the id and strand of the read.
"""

if not len(sys.argv) == 3:
    print(helpMessage)
    exit(1)

readId = int(sys.argv[1])
strand = int(sys.argv[2])
if not (strand==0 or strand==1):
    print('Invalid strand')
    print(helpMessage)
    exit(1)


fileName = 'Markers-' + str(readId) + '-' + str(strand) + '.csv'

a = shasta.Assembler()
a.accessKmers()
a.accessMarkers()
a.writeMarkers(readId=readId, strand=strand, fileName=fileName)
