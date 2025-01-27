
#include "ConfigurationTable.hpp"

namespace shasta {
   const std::vector< pair<string, string> > configurationTable = {
    {"Nanopore-Dec2019", R"zzz(# This file contains Shasta options that, as of December 2019,
# are known to work with Oxford Nanopore reads under the following 
# circumstances:

# - Human genome assembly.
# - Coverage between 40x and 80x. If you have more coverage than that,
#   you can use option "--Reads.minReadLength" to adjust coverage as desired.
# - Reads from Guppy 3.0.5 base caller. Also known to work with
#   reads from other Guppy releases 3.0.x and 3.1.x.

# To use this configuration file, specify Shasta option "--config PathToThisFile". 
# If you specify any conflicting values on the command line,
# the values specified on the command line take precedence.

# In most cases, for best performance on a large assembly 
# you will usually also want to use the following options, which 
# cannot be specified in a configuration file:
# --memoryMode filesystem
# --memoryBacking 2M
# Using these options require root access.

# For detailed information on all available options see here:
# https://paoloshasta.github.io/shasta/CommandLineOptions.html

# For information on running a small assembly for which 
# performance is not essential see here:
# https://paoloshasta.github.io/shasta/QuickStart.html

# For more information on running an assembly see here:
# https://paoloshasta.github.io/shasta/Running.html

# For information on optimizing assembly performance see here:
# https://paoloshasta.github.io/shasta/Performance.html



[Reads]
# If you have extra coverage, use this option to adjust coverage.
minReadLength = 10000

[MinHash]
minBucketSize = 5
maxBucketSize = 30
minFrequency = 5

[Align]
minAlignedFraction = 0.4

[Assembly]
consensusCaller = Bayesian:guppy-3.0.5-a

)zzz"},
    {"Nanopore-UL-Dec2019", R"zzz(# This file contains Shasta options that, as of December 2019,
# are known to work with Ultra-Long (UL) Oxford Nanopore reads 
# under the following circumstances:

# - Human genome assembly.
# - Coverage between 40x and 80x. If you have more coverage than that,
#   you can use option "--Reads.minReadLength" to adjust coverage as desired.
# - Reads from Guppy 3.0.5 base caller. Also known to work with
#   reads from other Guppy releases 3.0.x and 3.1.x.

# To use this configuration file, specify Shasta option "--config PathToThisFile". 
# If you specify any conflicting values on the command line,
# the values specified on the command line take precedence.

# In most cases, for best performance on a large assembly 
# you will usually also want to use the following options, which 
# cannot be specified in a configuration file:
# --memoryMode filesystem
# --memoryBacking 2M
# Using these options require root access.

# For detailed information on all available options see here:
# https://paoloshasta.github.io/shasta/CommandLineOptions.html

# For information on running a small assembly for which 
# performance is not essential see here:
# https://paoloshasta.github.io/shasta/QuickStart.html

# For more information on running an assembly see here:
# https://paoloshasta.github.io/shasta/Running.html

# For information on optimizing assembly performance see here:
# https://paoloshasta.github.io/shasta/Performance.html



[Reads]
# If you have extra coverage, use this option to adjust coverage.
minReadLength = 50000

[MinHash]
minBucketSize = 5
maxBucketSize = 40
minFrequency = 10

[Align]
maxSkip = 60
maxDrift = 60
minAlignedMarkerCount = 400

[Assembly]
consensusCaller = Bayesian:guppy-3.0.5-a

)zzz"},
    {"Nanopore-Sep2020", R"zzz(# DO NOT USE THIS FILE IF YOU HAVE READS CREATED BY A
# GUPPY VERSION OLDER THAN 3.6.0.

# This file contains Shasta options which attempt to partially automate
# parameter selection. It is based on an earlier config, which, as of Jun 2020,
# was known to work with Oxford Nanopore reads under the following circumstances:

# - Human genome assembly.
# - Coverage between 40x and 80x. If you have more coverage than that,
#   you can use option "--Reads.minReadLength" to adjust coverage as desired.
# - Reads from Guppy 3.6.0 base caller. If you have reads
#   from an older version of Guppy, use configuration file
#   Nanopore-Dec2019.conf instead.

# The automation provided by this config is particularly applicable to
# low coverage or non-human samples. It also matches or exceeds continuity
# in human samples, relative to the appropriately chosen 3.6.0 or 3.6.0-UL conf.
# Automation can also be activated with parameters designed for earlier basecallers,
# if needed, but updating to guppy 3.6.0 or higher will greatly improve assembly
# quality and is therefore strongly recommended.

# To use this configuration file, specify Shasta option 
# "--config AbsolutePathToThisFile". 
# If you specify any conflicting values on the command line,
# the values specified on the command line take precedence.

# In most cases, for best performance on a large assembly 
# you will usually also want to use the following options, which 
# cannot be specified in a configuration file:
# --memoryMode filesystem
# --memoryBacking 2M
# Using these options requires root access.

# For detailed information on all available options see here:
# https://paoloshasta.github.io/shasta/CommandLineOptions.html

# For information on running a small assembly for which 
# performance is not essential see here:
# https://paoloshasta.github.io/shasta/QuickStart.html

# For more information on running an assembly see here:
# https://paoloshasta.github.io/shasta/Running.html

# For information on optimizing assembly performance see here:
# https://paoloshasta.github.io/shasta/Performance.html



[Reads]
# If you have extra coverage, use this option to adjust coverage.
minReadLength = 10000
noCache = True

[Kmers]
# Due to the higher accuracy of Guppy 3.6.0 we use longer
# markers than usual.
k = 14

[MinHash]
minBucketSize = 5
maxBucketSize = 30
minFrequency = 5

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30

# The following Align parameters are set to very permissive values to allow the majority of alignments
# to be assessed during the initial stage of automatic alignment parameter selection
maxSkip = 100
maxDrift = 100
maxTrim = 100
minAlignedMarkerCount = 10
minAlignedFraction = 0.1

[ReadGraph]
# This method uses the observed distribution of alignment stats to choose a cutoff for
# maxSkip, maxDrift, maxTrim, minAlignedMarkerCount, and minAlignedFraction
creationMethod = 2

[MarkerGraph]
simplifyMaxLength = 10,100,1000,10000,100000
crossEdgeCoverageThreshold = 3

# Automatically determine this using PeakFinder
minCoverage = 0

[Assembly]
consensusCaller = Bayesian:guppy-3.6.0-a
detangleMethod = 2


)zzz"},
    {"Nanopore-UL-Sep2020", R"zzz(# DO NOT USE THIS FILE IF YOU HAVE READS CREATED BY A
# GUPPY VERSION OLDER THAN 3.6.0.

# This file contains Shasta options which attempt to partially automate
# parameter selection. It is based on an earlier config, which, as of Jun 2020,
# was known to work with Oxford Nanopore reads under the following circumstances:

# - Human genome assembly.
# - Coverage between 40x and 80x. If you have more coverage than that,
#   you can use option "--Reads.minReadLength" to adjust coverage as desired.
# - Reads from Guppy 3.6.0 base caller. If you have reads
#   from an older version of Guppy, use configuration file
#   Nanopore-UL-Dec2019.conf instead.

# The automation provided by this config is particularly applicable to
# low coverage or non-human samples. It also matches or exceeds continuity
# in human samples, relative to the appropriately chosen 3.6.0 or 3.6.0-UL conf.
# Automation can be activated with parameters designed for earlier basecallers,
# if needed, but updating to guppy 3.6.0 or higher will greatly improve assembly
# quality and is therefore strongly recommended.

# In most cases, for best performance on a large assembly 
# you will usually also want to use the following options, which 
# cannot be specified in a configuration file:
# --memoryMode filesystem
# --memoryBacking 2M
# Using these options requires root access.

# For detailed information on all available options see here:
# https://paoloshasta.github.io/shasta/CommandLineOptions.html

# For information on running a small assembly for which 
# performance is not essential see here:
# https://paoloshasta.github.io/shasta/QuickStart.html

# For more information on running an assembly see here:
# https://paoloshasta.github.io/shasta/Running.html

# For information on optimizing assembly performance see here:
# https://paoloshasta.github.io/shasta/Performance.html



[Reads]
# If you have extra coverage, use this option to adjust coverage.
minReadLength = 50000
noCache = True

[Kmers]
# Due to the higher accuracy of Guppy 3.6.0 we use longer
# markers than usual.
k = 14

[MinHash]
minBucketSize = 10
maxBucketSize = 30
minFrequency = 5

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30

# The following Align parameters are set to very permissive values to allow the majority of alignments
# to be assessed during the initial stage of automatic alignment parameter selection
maxSkip = 100
maxDrift = 100
maxTrim = 100
minAlignedMarkerCount = 10
minAlignedFraction = 0.1

[ReadGraph]
# This method uses the observed distribution of alignment stats to choose a cutoff for
# maxSkip, maxDrift, maxTrim, minAlignedMarkerCount, and minAlignedFraction
creationMethod = 2

[MarkerGraph]
simplifyMaxLength = 10,100,1000,10000,100000
crossEdgeCoverageThreshold = 3

# Automatically determine this using PeakFinder
minCoverage = 0

[Assembly]
consensusCaller = Bayesian:guppy-3.6.0-a
detangleMethod = 2


)zzz"},
    {"Nanopore-UL-iterative-Sep2020", R"zzz(# This configuration file is EXPERIMENTAL and should only
# be used under the following conditions:
# - Nanopore reads created by Guppy 3.6.0 or newer.
# - Ultra-Long (UL) reads with typical N50 80Kb or more.
# - High coverage 80X.
# Iterative assembly results in some separation of haplotypes
# and does better at resolving segmental duplications.

[Reads]
minReadLength = 30000 
noCache = True

[Kmers]
k = 10 

[MinHash]
minBucketSize = 10 
maxBucketSize = 40 
minFrequency = 5 

[Align]
alignMethod = 3 
matchScore = 6
gapScore = -3 
downsamplingFactor = 0.05 
maxSkip = 100
maxDrift = 100
maxTrim = 100
minAlignedMarkerCount = 10
minAlignedFraction = 0.1
sameChannelReadAlignment.suppressDeltaThreshold = 30 

[ReadGraph]
maxAlignmentCount = 12
creationMethod = 2

[MarkerGraph]
minCoveragePerStrand = 3
simplifyMaxLength = 10,100
crossEdgeCoverageThreshold = 3 

[Assembly]
detangleMethod = 2 
consensusCaller = Bayesian:guppy-3.6.0-a
iterative = True

)zzz"},
    {"Nanopore-OldGuppy-Sep2020", R"zzz(# This file contains Shasta options which attempt to partially automate
# parameter selection. It is based on an earlier config, which, as of Jun 2020,
# was known to work with Oxford Nanopore reads under the following circumstances:

# - Human genome assembly.
# - Coverage between 40x and 80x. If you have more coverage than that,
#   you can use option "--Reads.minReadLength" to adjust coverage as desired.
# - Reads from Guppy 3.0.5 base caller. Also known to work with
#   reads from other Guppy releases 3.0.x and 3.1.x.

# The automation provided by this config is particularly applicable to
# low coverage or non-human samples. It also matches or exceeds continuity
# in human samples, relative to the appropriately chosen config file.
# Updating to guppy 3.6.0 or higher will greatly improve assembly
# quality and is therefore strongly recommended.

# To use this configuration file, specify Shasta option
# "--config AbsolutePathToThisFile".
# If you specify any conflicting values on the command line,
# the values specified on the command line take precedence.

# In most cases, for best performance on a large assembly
# you will usually also want to use the following options, which 
# cannot be specified in a configuration file:
# --memoryMode filesystem
# --memoryBacking 2M
# Using these options requires root access.

# For detailed information on all available options see here:
# https://paoloshasta.github.io/shasta/CommandLineOptions.html

# For information on running a small assembly for which 
# performance is not essential see here:
# https://paoloshasta.github.io/shasta/QuickStart.html

# For more information on running an assembly see here:
# https://paoloshasta.github.io/shasta/Running.html

# For information on optimizing assembly performance see here:
# https://paoloshasta.github.io/shasta/Performance.html



[Reads]
# If you have extra coverage, use this option to adjust coverage.
minReadLength = 10000
noCache = True

[MinHash]
minBucketSize = 5
maxBucketSize = 30
minFrequency = 5

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30

# The following Align parameters are set to very permissive values to allow the majority of alignments
# to be assessed during the initial stage of automatic alignment parameter selection
maxSkip = 100
maxDrift = 100
maxTrim = 100
minAlignedMarkerCount = 10
minAlignedFraction = 0.1

[ReadGraph]
# This method uses the observed distribution of alignment stats to choose a cutoff for
# maxSkip, maxDrift, maxTrim, minAlignedMarkerCount, and minAlignedFraction
creationMethod = 2

[MarkerGraph]
simplifyMaxLength = 10,100,1000,10000,100000
crossEdgeCoverageThreshold = 3

# Automatically determine this using PeakFinder
minCoverage = 0

[Assembly]
consensusCaller = Bayesian:guppy-3.0.5-a
detangleMethod = 2

)zzz"},
    {"Nanopore-Plants-Apr2021", R"zzz(# This was used in April 2021 for an assembly of Lolium Perenne
# under the following conditions:

# - Oxford Nanopore reads, basecalled with Guppy 4.0.14.
# - Coverage: 30x.
# - Reads N50: 62 Kb. 

# The assembly was run:

# - Using Shasta at commit 266d4c8e65ff408db2dfd1381bb14648be83dad3.
# - On AWS Graviton2 (ARM) r6g.16xlarge instance type, 512 GB, 64 vCPUs.
# - Using memory options --memoryMode filesystem --memoryBacking 2M.

# Assembly results:

# - Assembled 2.185 Mb of sequence. Estimated genome size 
#   is 2.46 Gb from k-mer analysis, 2.72 Gb from flow cytometry.
# - Assembly N50: 5.5 Mb.
# - Elapsed time for assembly: 95 minutes.

# Many thanks to Dario Copetti (Molecular Plant Breeding, ETH Zurich, Switzerland)
# for providing access to the reads in advance of publication. 

# Also see Shasta issue #200 for some discussion
# https://github.com/paoloshasta/shasta/issues/200



[Reads]
noCache = True

[MinHash]
minBucketSize = 5
maxBucketSize = 30
minHashIterationCount = 50
minFrequency = 5

[Align]
downsamplingFactor = 0.05
sameChannelReadAlignment.suppressDeltaThreshold = 30
maxSkip = 60
maxDrift = 20
maxTrim = 60
minAlignedMarkerCount = 200
minAlignedFraction = 0.3

[MarkerGraph]
simplifyMaxLength = 10,100,1000,10000,100000
crossEdgeCoverageThreshold = 3
minCoverage = 0

[Assembly]
consensusCaller = Bayesian:guppy-3.6.0-a
detangleMethod = 2


)zzz"},
    {"Nanopore-Oct2021", R"zzz(# This is known to work at least under the following conditions:
# - Oxford Nanopore reads.
# - Guppy 5 base caller.
# - Human genome.
# - Coverage 30x to 80x.

[Reads]
minReadLength = 10000
noCache = True

[Kmers]
k = 14

[MinHash]
minBucketSize = 5
maxBucketSize = 30
minFrequency = 5

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30

# The following Align parameters are set to very permissive values to allow the majority of alignments
# to be assessed during the initial stage of automatic alignment parameter selection
# (ReadGraph.creationMethod 2).
maxSkip = 100
maxDrift = 100
maxTrim = 100
minAlignedMarkerCount = 10
minAlignedFraction = 0.1

[ReadGraph]
# This uses the observed distribution of alignment statistics to choose thresholds for
# maxSkip, maxDrift, maxTrim, minAlignedMarkerCount, and minAlignedFraction
creationMethod = 2

[MarkerGraph]
simplifyMaxLength = 10,100,1000,10000,100000
crossEdgeCoverageThreshold = 3

# Adaptive estimation of coverage threshold to generate marker graph vertices.
minCoverage = 0

[Assembly]
consensusCaller = Bayesian:guppy-5.0.7-a
detangleMethod = 2


)zzz"},
    {"Nanopore-UL-Oct2021", R"zzz(# This is known to work at least under the following conditions:
# - Oxford Nanopore Ultra-Long (UL) reads.
# - Guppy 5 base caller.
# - Human genome.
# - Coverage 30x to 80x.



[Reads]
minReadLength = 50000
noCache = True

[Kmers]
k = 14

[MinHash]
minBucketSize = 10
maxBucketSize = 30
minFrequency = 5

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30

# The following Align parameters are set to very permissive values to allow the majority of alignments
# to be assessed during the initial stage of automatic alignment parameter selection
# (ReadGraph.creationMethod 2).
maxSkip = 100
maxDrift = 100
maxTrim = 100
minAlignedMarkerCount = 10
minAlignedFraction = 0.1

[ReadGraph]
# This uses the observed distribution of alignment statistics to choose thresholds for
# maxSkip, maxDrift, maxTrim, minAlignedMarkerCount, and minAlignedFraction
creationMethod = 2

[MarkerGraph]
simplifyMaxLength = 10,100,1000,10000,100000
crossEdgeCoverageThreshold = 3

# Adaptive estimation of coverage threshold to generate marker graph vertices.
minCoverage = 0

[Assembly]
consensusCaller = Bayesian:guppy-5.0.7-a
detangleMethod = 2


)zzz"},
    {"HiFi-Oct2021", R"zzz(# This was tested on the  "HG002 Data Freeze (v1.0) 
# Recommended downsampled data mix", a 34x HiFi dataset
# for HG002. See here for more information:
# https://github.com/human-pangenomics/HG002_Data_Freeze_v1.0#hg002-data-freeze-v10-recommended-downsampled-data-mix
# It assembled 3078 Mb of sequence with an N50 of 35 Mb.
# The quality of assembled sequence was estimated
# at Q = 46 dB via k-mer based methods. It is limited
# by the fact that this creates a haploid assembly
# in which one of the alleles of each heterozygous locus
# is removed.



[Reads]
minReadLength = 8000
noCache = True

[Kmers]
k = 14

[MinHash]
hashFraction = 0.05
minHashIterationCount = 100
minFrequency = 3
minBucketSize = 10
maxBucketSize = 60

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
minAlignedFraction = 0.97
minAlignedMarkerCount = 200
maxSkip = 6
maxDrift = 4
maxTrim = 2

[ReadGraph]
maxAlignmentCount = 30
maxChimericReadDistance = 2

[MarkerGraph]
minCoverage = 6
simplifyMaxLength = 10,100,1000,10000,100000
crossEdgeCoverageThreshold = 3

[Assembly]
consensusCaller = Modal
detangleMethod = 2


)zzz"},
    {"Nanopore-UL-Jan2022", R"zzz(# This is known to work at least under the following conditions:
# - Oxford Nanopore Ultra-Long (UL) reads with read N50 50 Kb or more.
# - Guppy 5 base caller with "super" accuracy.
# - Human genome.
# - Coverage 60x to 80x. If you have more coverage, 
#   adjust --Reads.minReadLength or --Reads.desiredCoverage
#   to bring coverage down to this range.
#   For a human genome you can set --Reads.desiredCoverage 200G. 

# Under the above conditions, this should give an assembly with
# N50 in the tens of Mb.



[Reads]
minReadLength = 50000
noCache = True

[Kmers]
k = 14

[MinHash]
minBucketSize = 10
maxBucketSize = 50
minFrequency = 5

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30

# The following Align parameters are set to very permissive values to allow the majority of alignments
# to be assessed during the initial stage of automatic alignment parameter selection
# (ReadGraph.creationMethod 2).
maxSkip = 100
maxDrift = 100
maxTrim = 100
minAlignedMarkerCount = 10
minAlignedFraction = 0.1

[ReadGraph]
# This uses the observed distribution of alignment statistics to choose thresholds for
# maxSkip, maxDrift, maxTrim, minAlignedMarkerCount, and minAlignedFraction
creationMethod = 2
maxAlignmentCount = 12
strandSeparationMethod = 2

[MarkerGraph]
simplifyMaxLength = 10,100,1000,10000,100000
crossEdgeCoverageThreshold = 3

# Adaptive estimation of coverage threshold to generate marker graph vertices.
minCoverage = 0

[Assembly]
consensusCaller = Bayesian:guppy-5.0.7-a
detangleMethod = 2


)zzz"},
    {"Nanopore-Phased-Jan2022", R"zzz(
# Configuration file that uses assembly mode 2
# to create a phased diploid assembly using 
# Nanopore reads.

# This is known to work at least under the following conditions:
# - Oxford Nanopore reads.
# - Guppy 5 base caller with "super" accuracy.
# - Human genome.
# - Coverage 40x to 80x. If you have more coverage, 
#   adjust --Reads.minReadLength or --Reads.desiredCoverage
#   to bring coverage down to this range.



[Reads]

# Don't use the RLE representation of the reads.
representation = 0

# Read length. Adjust as necessary,
# or use desiredCoverage option to get coverage 
# around 40x to 80x.
minReadLength = 10000

noCache = True



[Kmers]

# Because we are not using the RLE representation,
# use shorter k-mers.
k = 8

# Compensate for the non-RLE representation
# to get about the same number of markers
# we would get with the RLE representation and probability = 0.1.
probability = 0.07



[MinHash]
minBucketSize = 5
maxBucketSize = 30
minFrequency = 5



[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30

# Permissive alignment criteria as required for read graph creation method 2.
maxSkip = 100
maxDrift = 100
maxTrim = 100
minAlignedMarkerCount = 10
minAlignedFraction = 0.1



[ReadGraph]

# Automatic adjustment of alignment criteria.
creationMethod = 2

# Strict strand separation is required for Mode 2 (phased) assembly.
strandSeparationMethod = 2

maxAlignmentCount = 6



[MarkerGraph]
minCoverage = 6
minCoveragePerStrand = 1
minEdgeCoverage = 6
minEdgeCoveragePerStrand = 1



[Assembly]
mode = 2
consensusCaller = Bayesian:guppy-5.0.7-a
pruneLength = 100
mode2.bubbleRemoval.minConcordantReadCount = 2





)zzz"},
    {"Nanopore-UL-Phased-Jan2022", R"zzz(
# Configuration file that uses assembly mode 2
# to create a phased diploid assembly using 
# Nanopore Ultra-Long (UL) reads.

# This is known to work at least under the following conditions:
# - Oxford Nanopore Ultra-Long (UL) reads with read N50 50 Kb or more.
# - Guppy 5 base caller with "super" accuracy.
# - Human genome.
# - Coverage 60x to 80x. If you have more coverage, 
#   adjust --Reads.minReadLength or --Reads.desiredCoverage
#   to bring coverage down to this range.
#   For a human genome you can set --Reads.desiredCoverage 200G. 



[Reads]

# Don't use the RLE representation of the reads.
representation = 0

# Read length cutoff for UL reads, adjust as necessary,
# or use desiredCoverage option to get coverage 
# around 40x to 80x.
minReadLength = 50000

noCache = True



[Kmers]

# Because we are not using the RLE representation,
# use shorter k-mers.
k = 8

# Compensate for the non-RLE representation
# to get about the same number of markers
# we would get with the RLE representation and probability = 0.1.
probability = 0.07



[MinHash]
minBucketSize = 10
maxBucketSize = 50
minFrequency = 5



[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30

# Permissive alignment criteria as required for read graph creation method 2.
maxSkip = 100
maxDrift = 100
maxTrim = 100
minAlignedMarkerCount = 10
minAlignedFraction = 0.1



[ReadGraph]

# Automatic adjustment of alignment criteria.
creationMethod = 2

# Strict strand separation is required for Mode 2 (phased) assembly.
strandSeparationMethod = 2

maxAlignmentCount = 12



[MarkerGraph]
minCoverage = 6
minCoveragePerStrand = 1
minEdgeCoverage = 6
minEdgeCoveragePerStrand = 1



[Assembly]
mode = 2
consensusCaller = Bayesian:guppy-5.0.7-a
pruneLength = 100
mode2.bubbleRemoval.minConcordantReadCount = 2





)zzz"},
    {"Nanopore-May2022", R"zzz(# This is known to work at least under the following conditions:
# - Oxford Nanopore Ultra-Long reads.
# - Human genomes.
# - Guppy 5 base caller or newer with "super" accuracy.
# - Coverage 40x to 80x. If you have more coverage, 
#   adjust --Reads.minReadLength or --Reads.desiredCoverage
#   to bring coverage down to this range.

# Under the above conditions, this should give an assembly with
# N50 in the tens of Mb.


[Reads]
minReadLength = 10000
noCache = True

[Kmers]
k = 14

[MinHash]
minBucketSize = 5
maxBucketSize = 30
minFrequency = 5

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30

# The following Align parameters are set to very permissive values to allow the majority of alignments
# to be assessed during the initial stage of automatic alignment parameter selection
# (ReadGraph.creationMethod 2).
maxSkip = 100
maxDrift = 100
maxTrim = 100
minAlignedMarkerCount = 10
minAlignedFraction = 0.1

[ReadGraph]
# This uses the observed distribution of alignment statistics to choose thresholds for
# maxSkip, maxDrift, maxTrim, minAlignedMarkerCount, and minAlignedFraction
creationMethod = 2

[MarkerGraph]
simplifyMaxLength = 10,100,1000,10000,100000
crossEdgeCoverageThreshold = 3

# Adaptive estimation of coverage threshold to generate marker graph vertices.
minCoverage = 0

[Assembly]
consensusCaller = Bayesian:guppy-5.0.7-b
detangleMethod = 2


)zzz"},
    {"Nanopore-Phased-May2022", R"zzz(

# This is known to work at least under the following conditions:
# - Phased assembly.
# - Oxford Nanopore Ultra-Long reads.
# - Human genomes.
# - Guppy 5 base caller or newer with "super" accuracy.
# - Coverage 40x to 80x. If you have more coverage, 
#   adjust --Reads.minReadLength or --Reads.desiredCoverage
#   to bring coverage down to this range.



[Reads]
# Read length. Adjust as necessary,
# or use desiredCoverage option to get coverage 
# around 40x to 80x.
minReadLength = 10000

noCache = True



[MinHash]
minBucketSize = 5
maxBucketSize = 30
minFrequency = 5



[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30

# Permissive alignment criteria as required for read graph creation method 2.
maxSkip = 100
maxDrift = 100
maxTrim = 100
minAlignedMarkerCount = 10
minAlignedFraction = 0.1



[ReadGraph]

# Automatic adjustment of alignment criteria.
creationMethod = 2

# Strict strand separation is required for Mode 2 (phased) assembly.
strandSeparationMethod = 2

maxAlignmentCount = 6



[MarkerGraph]
minCoverage = 6
minCoveragePerStrand = 1
minEdgeCoverage = 6
minEdgeCoveragePerStrand = 1



[Assembly]
mode = 2
consensusCaller = Bayesian:guppy-5.0.7-b
pruneLength = 100
mode2.bubbleRemoval.minConcordantReadCount = 2





)zzz"},
    {"Nanopore-UL-May2022", R"zzz(# This is known to work at least under the following conditions:
# - Oxford Nanopore Ultra-Long (UL) reads with read N50 50 Kb or more.
# - Guppy 5 base caller or newer with "super" accuracy.
# - Human genome.
# - Coverage 60x to 80x. If you have more coverage, 
#   adjust --Reads.minReadLength or --Reads.desiredCoverage
#   to bring coverage down to this range.
#   For a human genome you can set --Reads.desiredCoverage 200G. 

# Under the above conditions, this should give an assembly with
# N50 in the tens of Mb.



[Reads]
minReadLength = 50000
noCache = True

[Kmers]
k = 14

[MinHash]
minBucketSize = 10
maxBucketSize = 50
minFrequency = 5

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30

# The following Align parameters are set to very permissive values to allow the majority of alignments
# to be assessed during the initial stage of automatic alignment parameter selection
# (ReadGraph.creationMethod 2).
maxSkip = 100
maxDrift = 100
maxTrim = 100
minAlignedMarkerCount = 10
minAlignedFraction = 0.1

[ReadGraph]
# This uses the observed distribution of alignment statistics to choose thresholds for
# maxSkip, maxDrift, maxTrim, minAlignedMarkerCount, and minAlignedFraction
creationMethod = 2
maxAlignmentCount = 12
strandSeparationMethod = 2

[MarkerGraph]
simplifyMaxLength = 10,100,1000,10000,100000
crossEdgeCoverageThreshold = 3

# Adaptive estimation of coverage threshold to generate marker graph vertices.
minCoverage = 0

[Assembly]
consensusCaller = Bayesian:guppy-5.0.7-b
detangleMethod = 2


)zzz"},
    {"Nanopore-UL-Phased-May2022", R"zzz(
# Configuration file that uses assembly mode 2
# to create a phased diploid assembly using 
# Nanopore Ultra-Long (UL) reads.

# This is known to work at least under the following conditions:
# - Oxford Nanopore Ultra-Long (UL) reads with read N50 50 Kb or more.
# - Guppy 5 base caller with "super" accuracy.
# - Human genome.
# - Coverage 60x to 80x. If you have more coverage, 
#   adjust --Reads.minReadLength or --Reads.desiredCoverage
#   to bring coverage down to this range.
#   For a human genome you can set --Reads.desiredCoverage 200G. 



[Reads]
# Read length cutoff for UL reads, adjust as necessary,
# or use desiredCoverage option to get coverage 
# around 40x to 80x.
minReadLength = 50000
noCache = True



[MinHash]
minBucketSize = 10
maxBucketSize = 50
minFrequency = 5



[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30

# Permissive alignment criteria as required for read graph creation method 2.
maxSkip = 100
maxDrift = 100
maxTrim = 100
minAlignedMarkerCount = 10
minAlignedFraction = 0.1



[ReadGraph]

# Automatic adjustment of alignment criteria.
creationMethod = 2

# Strict strand separation is required for Mode 2 (phased) assembly.
strandSeparationMethod = 2

maxAlignmentCount = 12



[MarkerGraph]
minCoverage = 6
minCoveragePerStrand = 1
minEdgeCoverage = 6
minEdgeCoveragePerStrand = 1



[Assembly]
mode = 2
consensusCaller = Bayesian:guppy-5.0.7-b
pruneLength = 100
mode2.bubbleRemoval.minConcordantReadCount = 2





)zzz"},
    {"Nanopore-Human-SingleFlowcell-May2022", R"zzz(# Shasta assembly configuration specialized for human assemblies
# at low coverage using a single flowcell of nanopore reads.
# It is meant to be used under the following conditions:
#
# - Human genome.
# - Nanopore reads, single flowcell (100-120 Gb of coverage or around 35x).
# - Guppy 5.0.7 base caller or newer, with "super" accuracy.
# - Read N50 around 30 Kb.
# - Haploid assembly.



[Reads]
minReadLength = 10000
noCache = True

[Kmers]
k = 14

[MinHash]
minBucketSize = 5
maxBucketSize = 30
minHashIterationCount = 100
minFrequency = 5

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30
maxSkip = 30
maxDrift = 15
maxTrim = 30
minAlignedMarkerCount = 200
minAlignedFraction = 0.6

[ReadGraph]
creationMethod = 0
maxAlignmentCount = 12

[MarkerGraph]
simplifyMaxLength = 10,100,1000,10000,100000
crossEdgeCoverageThreshold = 3

# Adaptive estimation of coverage threshold to generate marker graph vertices.
minCoverage = 0

[Assembly]
consensusCaller = Bayesian:guppy-5.0.7-b
detangleMethod = 2

)zzz"},
    {"Nanopore-Human-SingleFlowcell-Phased-May2022", R"zzz(# Shasta assembly configuration specialized for phased diploid 
# human assemblies at low coverage using
# a single flowcell of nanopore reads.
# Because of low coverage, it assembles relatively smalll diploid bubbles.
# It is meant to be used under the following conditions:
#
# - Human genome.
# - Nanopore reads, single flowcell (100-120 Gb of coverage or around 35x).
# - Guppy 5.0.7 base caller or newer, with "super" accuracy.
# - Read N50 around 30 Kb.
# - Phased diploid assembly.



[Reads]
noCache = True

[MinHash]
minBucketSize = 5
maxBucketSize = 30
minHashIterationCount = 100
minFrequency = 5

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30
maxSkip = 30
maxDrift = 15
maxTrim = 30
minAlignedMarkerCount = 200
minAlignedFraction = 0.6

[ReadGraph]
strandSeparationMethod = 2
maxAlignmentCount = 12

[MarkerGraph]
minCoverage = 6
minCoveragePerStrand = 1
minEdgeCoverage = 6
minEdgeCoveragePerStrand = 1

[Assembly]
mode = 2
consensusCaller = Bayesian:guppy-5.0.7-b
pruneLength = 100
mode2.bubbleRemoval.minConcordantReadCount = 2





)zzz"},
    {"Nanopore-UL-Phased-Nov2022", R"zzz(
# Configuration file that uses assembly mode 2
# to create a phased diploid assembly using 
# Nanopore Ultra-Long (UL) reads.

# This is known to work at least under the following conditions:
# - Oxford Nanopore Ultra-Long (UL) reads with read N50 50 Kb or more.
# - Guppy 5 or 6 base caller with "super" accuracy.
# - Human genome.
# - Coverage 60x to 80x. If you have more coverage, 
#   adjust --Reads.minReadLength or --Reads.desiredCoverage
#   to bring coverage down to this range.
#   For a human genome you can set --Reads.desiredCoverage 200G. 
# Compared to the Nanopore-UL-Phased-May2022 configuration,
# this provides better phasing accuracy without significantly decreasing
# the size of phased bubbles.



[Reads]
# Read length cutoff for UL reads, adjust as necessary,
# or use desiredCoverage option to get coverage 
# around 40x to 80x.
minReadLength = 50000
noCache = True



[MinHash]
minBucketSize = 10
maxBucketSize = 50
minFrequency = 5



[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30

# Permissive alignment criteria as required for read graph creation method 2.
maxSkip = 100
maxDrift = 100
maxTrim = 100
minAlignedMarkerCount = 10
minAlignedFraction = 0.1



[ReadGraph]

# Automatic adjustment of alignment criteria.
creationMethod = 2

# Strict strand separation is required for Mode 2 (phased) assembly.
strandSeparationMethod = 2

maxAlignmentCount = 12



[MarkerGraph]
minCoverage = 6
minCoveragePerStrand = 1
minEdgeCoverage = 6
minEdgeCoveragePerStrand = 1



[Assembly]
mode = 2
consensusCaller = Bayesian:guppy-5.0.7-b
pruneLength = 100
mode2.bubbleRemoval.minConcordantReadCount = 2
mode2.phasing.minLogP = 50





)zzz"},
    {"Nanopore-R10-Fast-Nov2022", R"zzz(# This Shasta assembly configuration was tested under the following conditions:
#  - ONT R10 chemistry, fast mode, chimera rate around 2%.
#  - Guppy 6 basecaller with "super" accuracy.
#  - Human genome at low coverage (one flowcell, coverage around 35x).
#  - Haploid assembly.
# Under these conditions it produced a haploid assembly with N50 around 20 Mb
# and with base level accuracy limited mostly by heterozygosity 
# of the human genome.  
# This configuration might also be usable for other genomes and at higher coverage,
# but that was not tested.

[Reads]
representation = 0
minReadLength = 10000
noCache = True

[Kmers]
k = 14

[MinHash]
minHashIterationCount = 100
minBucketSize = 10
maxBucketSize = 40
minFrequency = 5

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30
minAlignedMarkerCount = 1000
minAlignedFraction = 0.85
maxSkip = 20
maxDrift = 10
maxTrim = 20

[ReadGraph]
creationMethod = 0
maxAlignmentCount = 15

[MarkerGraph]
simplifyMaxLength = 10,100,1000,10000,100000
crossEdgeCoverageThreshold = 3

# Adaptive estimation of coverage threshold to generate marker graph vertices.
minCoverage = 0

[Assembly]
consensusCaller = Modal
detangleMethod = 2


)zzz"},
    {"Nanopore-R10-Slow-Nov2022", R"zzz(# This Shasta assembly configuration was tested under the following conditions:
#  - ONT R10 chemistry, slow mode, chimera rate around 2%.
#  - Guppy 6 basecaller with "super" accuracy.
#  - Human genome at medium coverage (two flowcells, coverage around 45x).
#  - Haploid assembly.
# Under these conditions it produced a haploid assembly with N50 around 30 Mb
# and with base level accuracy limited mostly by heterozygosity 
# of the human genome.  
# This configuration might also be usable for other genomes 
# and at higher or lower  coverage, but that was not tested.

[Reads]
representation = 0
minReadLength = 10000
noCache = True

[Kmers]
k = 14

[MinHash]
minHashIterationCount = 100
minBucketSize = 20
maxBucketSize = 60
minFrequency = 5

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30
minAlignedMarkerCount = 1200
minAlignedFraction = 0.9
maxSkip = 12
maxDrift = 8
maxTrim = 10

[ReadGraph]
creationMethod = 0
maxAlignmentCount = 15

[MarkerGraph]
simplifyMaxLength = 10,100,1000,10000,100000
crossEdgeCoverageThreshold = 3

# Adaptive estimation of coverage threshold to generate marker graph vertices.
minCoverage = 0

[Assembly]
consensusCaller = Modal
detangleMethod = 2


)zzz"},
    {"Nanopore-Phased-R10-Fast-Nov2022", R"zzz(# This Shasta assembly configuration for phased diploid assembly
# was tested under the following conditions:

#  - ONT R10 chemistry, fast mode, chimera rate around 2%.
#  - Guppy 6 basecaller with "super" accuracy.
#  - Human genome at low coverage (HG002, one flowcell, coverage around 35x).
#  - Phased diploid assembly.

# Under these conditions, a test run for HG002 produced an assembly 
# consisting of about 2.7 Gb in bubble chains and 0.3 Gb outside bubble chains.
# Of the 2.7 Gb in bubble chains, about 2.0 Gb were assembled diploid and phased. 

# The N50 for bubble chains was about 5 Mb, and the N50 for phased bubbles was about 0.5 Mb.
# This is longer than most genes, which means that for many genes 
# the assembly contains both haplotypes, entirely phased.

# When mapping each branch of a phased bubble against the correct reference haplotype
# for HG002, base level quality was around Q = 47 dB for mismatches,
# and around Q = 40 dB for indels.

# The fraction assembled diploid and phased can be improved with the use of 
# Ultra-Long reads (e. g. 2.7 Gb were assembled diploid and phased with R9, 
# Ultra-Long reads at high coverage). A separate assembly configuration
# for phased diploid assembly using R10 Ultra-Long reads will
# be provided when possible. 

# This configuration might also be usable under different conditions,
# but that was not tested.



[Reads]
representation = 0
minReadLength = 10000
noCache = True

[Kmers]
k = 14

[MinHash]
minHashIterationCount = 100
minBucketSize = 10
maxBucketSize = 40
minFrequency = 5

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30
minAlignedMarkerCount = 1000
minAlignedFraction = 0.85
maxSkip = 20
maxDrift = 10
maxTrim = 20

[ReadGraph]
creationMethod = 0
maxAlignmentCount = 15
strandSeparationMethod = 2

[MarkerGraph]
minCoverage = 6
minCoveragePerStrand = 1
minEdgeCoverage = 6
minEdgeCoveragePerStrand = 1

[Assembly]
mode = 2
consensusCaller = Modal
pruneLength = 100
mode2.bubbleRemoval.minConcordantReadCount = 2




)zzz"},
    {"Nanopore-Phased-R10-Slow-Nov2022", R"zzz(# This Shasta assembly configuration for phased diploid assembly
# was tested under the following conditions:

#  - ONT R10 chemistry, slow mode, chimera rate around 2%.
#  - Guppy 6 basecaller with "super" accuracy.
#  - Human genome at medium coverage (HG002, two flowcells, coverage around 45x).
#  - Phased diploid assembly.

# Under these conditions, a test run for HG002 produced an assembly 
# consisting of about 2.7 Gb in bubble chains and 0.3 Gb outside bubble chains.
# Of the 2.7 Gb in bubble chains, about 2.2 Gb were assembled diploid and phased. 

# The N50 for bubble chains and the N50 for phased bubbles was about 10 Mb.


# When mapping each branch of a phased bubble against the correct reference haplotype
# for HG002, base level quality was around Q = 50 dB for mismatches,
# and around Q = 40 dB for indels.

# The fraction assembled diploid and phased can be improved with the use of 
# Ultra-Long reads (e. g. 2.7 Gb were assembled diploid and phased with R9, 
# Ultra-Long reads at high coverage). A separate assembly configuration
# for phased diploid assembly using R10 Ultra-Long reads will
# be provided when possible. 

# This configuration might also be usable under different conditions,
# but that was not tested.



[Reads]
representation = 0
minReadLength = 10000
noCache = True

[Kmers]
k = 14

[MinHash]
minHashIterationCount = 100
minBucketSize = 20
maxBucketSize = 60
minFrequency = 5

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30
minAlignedMarkerCount = 1200
minAlignedFraction = 0.9
maxSkip = 12
maxDrift = 8
maxTrim = 10

[ReadGraph]
creationMethod = 0
strandSeparationMethod = 2
maxAlignmentCount = 15

[MarkerGraph]
minCoverage = 6
minCoveragePerStrand = 1
minEdgeCoverage = 6
minEdgeCoveragePerStrand = 1


[Assembly]
mode = 2
consensusCaller = Modal
pruneLength = 100
mode2.bubbleRemoval.minConcordantReadCount = 2



)zzz"},
    {"Nanopore-ncm23-May2024", R"zzz(# This assembly configuration is for nanopore reads generated using the 
# "Experimental extremely high-accuracy, ultra-long sequencing kit"
# from the ONT December 2023 data release:
# https://labs.epi2me.io/gm24385_ncm23_preview/

# It uses Mode 3 assembly to create a phased assembly.
# It was only tested for a human genome at coverage 40x to 60x,
# but it should work at lower or higher coverage,
# within reasonable limits, because it includes some
# provisions for coverage adaptivity.

[Reads]
representation = 0
minReadLength = 10000
noCache = True
palindromicReads.deltaThreshold = 300

[Kmers]
k = 30
probability = 0.05

[MinHash]
minHashIterationCount = 50
minBucketSize = 0
maxBucketSize = 0
minFrequency = 5

[Align]
alignMethod = 5
sameChannelReadAlignment.suppressDeltaThreshold = 30
minAlignedMarkerCount = 1000
minAlignedFraction = 0.9
maxSkip = 20
maxDrift = 10
maxTrim = 20

[ReadGraph]
maxAlignmentCount = 20
strandSeparationMethod = 2

[Assembly]
mode = 3


)zzz"},
    {"Nanopore-r10.4.1_e8.2-400bps_sup-Herro-Sep2024", R"zzz(# Shasta assembly configuration for ONT r10.4.1_e8.2-400bps_sup reads,
# error-corrected with HERRO, as announced by Oxford Nanopore in May 2024:
# see https://labs.epi2me.io/lc2024_t2t/

[Reads]
representation = 0
minReadLength = 10000
noCache = True
palindromicReads.deltaThreshold = 300

[Kmers]
k = 30
probability = 0.05

[MinHash]
minHashIterationCount = 50
minBucketSize = 0
maxBucketSize = 0
minFrequency = 5

[Align]
alignMethod = 5
sameChannelReadAlignment.suppressDeltaThreshold = 30
minAlignedMarkerCount = 700
minAlignedFraction = 0.99
maxSkip = 20
maxDrift = 10
maxTrim = 20
align5.driftRateTolerance = 0.05
align5.minBandExtend = 10

[ReadGraph]
maxAlignmentCount = 40
strandSeparationMethod = 2

[Assembly]
mode = 3
mode3.primaryGraph.maxLoss = 0.5

)zzz"},
    {"Nanopore-r10.4.1_e8.2-400bps_sup-Raw-Sep2024", R"zzz(# Shasta assembly configuration for ONT r10.4.1_e8.2-400bps_sup raw reads,
# without error correction, as announced by Oxford Nanopore in May 2024:
# see https://labs.epi2me.io/lc2024_t2t/

[Reads]
representation = 0
minReadLength = 10000
noCache = True
palindromicReads.deltaThreshold = 300

[Kmers]
k = 30
probability = 0.05

[MinHash]
m = 10
minHashIterationCount = 50
minBucketSize = 0
maxBucketSize = 0
minFrequency = 70

[Align]
alignMethod = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30
minAlignedMarkerCount = 500
minAlignedFraction = 0.83
maxSkip = 10
maxDrift = 14
maxTrim = 10
align6.maxInBandCount = 1000000
align6.maxGlobalFrequencyMultiplier = 0.6
align6.minLowFrequencyCount = 6

[ReadGraph]
maxAlignmentCount = 20
strandSeparationMethod = 2
preferAlignedFraction = True

[Assembly]
mode = 3
mode3.primaryGraph.maxLoss = 0.5
mode3.assemblyGraph.phaseErrorThreshold = 0.2
mode3.assemblyGraph.minLogP = 10


)zzz"},
    {"Nanopore-r10.4.1_e8.2-400bps_sup-Herro-Jan2025", R"zzz(# Shasta assembly configuration for ONT r10.4.1_e8.2-400bps_sup reads,
# error-corrected with HERRO, as announced by Oxford Nanopore in May 2024:
# see https://labs.epi2me.io/lc2024_t2t/

[Reads]
representation = 0
minReadLength = 10000
noCache = True
palindromicReads.deltaThreshold = 300

[Kmers]
k = 50
probability = 0.05

[MinHash]
minHashIterationCount = 50
minBucketSize = 0
maxBucketSize = 0
minFrequency = 5

[Align]
alignMethod = 5
sameChannelReadAlignment.suppressDeltaThreshold = 30
minAlignedMarkerCount = 500
minAlignedFraction = 0.98
maxSkip = 50
maxDrift = 50
maxTrim = 50
align5.driftRateTolerance = 0.05
align5.minBandExtend = 10

[ReadGraph]
creationMethod = 4
maxAlignmentCount = 10
strandSeparationMethod = 0
maxChimericReadDistance = 0
epsilon = 1e-4
delta = 5e-4
WThreshold = 1e-8
WThresholdForBreaks = 1e+15

[Assembly]
mode = 3
mode3.primaryGraph.maxLoss = 0.5
mode3.minAnchorCoverage = 0
mode3.maxAnchorCoverage = 0
mode3.maxAnchorCoverageMultiplier = 1.18
mode3.assemblyGraph.pruneLength = 150000

)zzz"},
    {"Nanopore-r10.4.1_e8.2-400bps_sup-Raw-Jan2025", R"zzz(# Shasta assembly configuration for ONT r10.4.1_e8.2-400bps_sup raw reads,
# without error correction, as announced by Oxford Nanopore in May 2024:
# see https://labs.epi2me.io/lc2024_t2t/

[Reads]
representation = 0
minReadLength = 10000
noCache = True
palindromicReads.deltaThreshold = 300

[Kmers]
k = 26
probability = 0.04

[MinHash]
m = 10
minHashIterationCount = 50
minBucketSize = 0
maxBucketSize = 0
minFrequency = 70

[Align]
alignMethod = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30
minAlignedMarkerCount = 500
minAlignedFraction = 0.83
maxSkip = 10
maxDrift = 14
maxTrim = 10
align6.maxInBandCount = 1000000
align6.maxGlobalFrequencyMultiplier = 0.6
align6.minLowFrequencyCount = 6

[ReadGraph]
maxAlignmentCount = 20
strandSeparationMethod = 2
preferAlignedFraction = True

[Assembly]
mode = 3
mode3.primaryGraph.maxLoss = 0.5
mode3.assemblyGraph.phaseErrorThreshold = 0.2
mode3.assemblyGraph.minLogP = 10


)zzz"}
    };
}
