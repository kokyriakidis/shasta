#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "Reads.hpp"
#include "seqan.hpp"
using namespace shasta;

void Assembler::alignOrientedReads5(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    int matchScore,
    int mismatchScore,
    int gapScore,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    ostream& html)
{

    // EXPOSE WHEN CODE STABILIZES *******
    const uint64_t maxMarkerFrequency = 1;
    const double driftRateTolerance = 0.02;

    // Get the marker KmerIds for the two oriented reads.
    array<span<KmerId>, 2> allMarkerKmerIds;
    array<vector<KmerId>, 2> allMarkerKmerIdsVectors;
    if(markerKmerIds.isOpen()) {
        allMarkerKmerIds[0] = markerKmerIds[orientedReadId0.getValue()];
        allMarkerKmerIds[1] = markerKmerIds[orientedReadId1.getValue()];
    } else {
        // This is slower and will happen if markerKmerIds is not available.
        // Resize the vectors and make the spans point to the vectors.
        // Then call getOrientedReadMarkerKmerIds to fill them in.
        allMarkerKmerIdsVectors[0].resize(markers.size(orientedReadId0.getValue()));
        allMarkerKmerIdsVectors[1].resize(markers.size(orientedReadId1.getValue()));
        allMarkerKmerIds[0] = span<KmerId>(allMarkerKmerIdsVectors[0]);
        allMarkerKmerIds[1] = span<KmerId>(allMarkerKmerIdsVectors[1]);
        getOrientedReadMarkerKmerIds(orientedReadId0, allMarkerKmerIds[0]);
        getOrientedReadMarkerKmerIds(orientedReadId1, allMarkerKmerIds[1]);
    }


    // Get the low frequency markers in the two oriented reads, sorted by KmerId.
    array< span<uint32_t>, 2> lowFrequencyOrdinals;
    array< vector<uint32_t>, 2> lowFrequencyOrdinalsVectors;
    if(lowFrequencyMarkers.isOpen()) {
        // Use the stored copy.
        lowFrequencyOrdinals[0] = lowFrequencyMarkers[orientedReadId0.getValue()];
        lowFrequencyOrdinals[1] = lowFrequencyMarkers[orientedReadId1.getValue()];
    }
    else {
        // Compute them and store in the local vectors, then have the spans point to them.
        for(uint64_t i=0; i<2; i++) {
            computeLowFrequencyMarkers(allMarkerKmerIds[i], maxMarkerFrequency, lowFrequencyOrdinalsVectors[i]);
            lowFrequencyOrdinals[i] = span<uint32_t>(lowFrequencyOrdinalsVectors[i]);
        }
    }



    if(html) {
        for(uint64_t i=0; i<2; i++) {
            html << "<br>" << (i==0 ? orientedReadId0 : orientedReadId1) << " has " << allMarkerKmerIds[i].size() <<
                " markers of which " << lowFrequencyOrdinals[i].size() << " appear up to " <<
                maxMarkerFrequency << " times." << endl;
        }
    }



    // Find pairs of ordinals in the two oriented reads that correspond to
    // the same low frequency k-mers.
    class CommonKmerInfo {
    public:
        uint32_t ordinal0;
        uint32_t ordinal1;
        KmerId kmerId;
        uint64_t rank0 = invalid<uint64_t>;
        uint64_t rank1 = invalid<uint64_t>;
        uint64_t ordinalSum() const
        {
            return ordinal0 + ordinal1;
        }
        int64_t ordinalOffset() const
        {
            return int64_t(ordinal0) - int64_t(ordinal1);
        }
    };
    vector<CommonKmerInfo> commonKmerInfos;

    // Joint loop over the ordinals corresponding to low frequency markers.
    // They are both sorted by KmerId.
    const auto begin0 = lowFrequencyOrdinals[0].begin();
    const auto begin1 = lowFrequencyOrdinals[1].begin();
    const auto end0 = lowFrequencyOrdinals[0].end();
    const auto end1 = lowFrequencyOrdinals[1].end();
    auto it0 = begin0;
    auto it1 = begin1;
    while((it0 != end0) and (it1 != end1)) {
        const uint32_t ordinal0 = *it0;
        const uint32_t ordinal1 = *it1;
        const KmerId kmerId0 = allMarkerKmerIds[0][ordinal0];
        const KmerId kmerId1 = allMarkerKmerIds[1][ordinal1];

        if(kmerId0 < kmerId1) {

            // Go past the streak with this KmerId in lowFrequencyOrdinals[0].
            while(it0 != end0 and allMarkerKmerIds[0][*it0] == kmerId0) {
                ++it0;
            }


        } else if(kmerId1 < kmerId0) {

            // Go past the streak with this KmerId in lowFrequencyOrdinals[1].
            while(it1 != end1 and allMarkerKmerIds[1][*it1] == kmerId1) {
                ++it1;
            }

        } else {

            // We found a common low frequency marker k-mer.
            SHASTA_ASSERT(kmerId0 == kmerId1);
            const KmerId kmerId = kmerId0;

            // Look for the streak with this KmerId in lowFrequencyOrdinals[0].
            auto streakBegin0 = it0;
            auto streakEnd0 = it0 + 1;
            while(streakEnd0 != end0 and allMarkerKmerIds[0][*streakEnd0] == kmerId) {
                ++streakEnd0;
            }

            // Look for the streak with this KmerId in lowFrequencyOrdinals[1].
            auto streakBegin1 = it1;
            auto streakEnd1 = it1 + 1;
            while(streakEnd1 != end1 and allMarkerKmerIds[1][*streakEnd1] == kmerId) {
                ++streakEnd1;
            }

            // Look over pairs of markers in these streaks.
            for(auto jt0=streakBegin0; jt0!=streakEnd0; jt0++) {
                for(auto jt1=streakBegin1; jt1!=streakEnd1; jt1++) {
                    commonKmerInfos.push_back({*jt0, *jt1, kmerId});
                }
            }

            // Point to the next marker in lowFrequencyOrdinals[0] and lowFrequencyOrdinals[1].
            it0 = streakEnd0;
            it1 = streakEnd1;
        }
    }



    // Write the common unique markers.
    if(html) {
        html << "<h3>Common unique markers</h3>";
        html << "There are " << commonKmerInfos.size() << " common unique markers." << endl;
        html << "<p><table>"
            "<tr><th>Ordinal0<th>Ordinal1<th>Ordinal<br>offset<th>Ordinal<br>sum<th>KmerId<th>Kmer";
        const uint64_t k = assemblerInfo->k;
        for(const CommonKmerInfo& commonKmerInfo: commonKmerInfos) {
            const Kmer kmer(commonKmerInfo.kmerId, k);
            html << "<tr>"
                "<td class=centered>" << commonKmerInfo.ordinal0 <<
                "<td class=centered>" << commonKmerInfo.ordinal1 <<
                "<td class=centered>" << commonKmerInfo.ordinalOffset() <<
                "<td class=centered>" << commonKmerInfo.ordinalSum();

            // Write the KmerId in hex with the appropriate number of digits.
            const char oldFill = html.fill('0');
            html << "<td class=centered style='font-family:monospace'>" <<
                std::hex << std::setw(int(k/2)) << commonKmerInfo.kmerId << std::dec;
            html.fill(oldFill);

            // Write the Kmer.
            html << "<td class=centered style='font-family:monospace'>";
            kmer.write(html, k);
        }
        html << "</table>";
    }



    // Create a histogram of ordinal offsets for the common unique markers.
    std::map<int64_t, uint64_t> histogramMap;
    for(const CommonKmerInfo& commonKmerInfo: commonKmerInfos) {
        const int64_t offset = commonKmerInfo.ordinalOffset();
        auto it = histogramMap.find(offset);
        if(it == histogramMap.end()) {
            histogramMap.insert({offset, 1});
        } else {
            ++it->second;
        }
    }
    vector< pair<int64_t, uint64_t> > histogram;
    copy(histogramMap.begin(), histogramMap.end(), back_inserter(histogram));
    if(html) {
        html << "<h3>Histogram of ordinal offsets for the common unique markers</h3>"
            "<table>"
            "<tr><th>Ordinal<br>offset<th>Frequency";
        for(const auto& p: histogram) {
            html << "<tr>"
                "<td class=centered>" << p.first <<
                "<td class=centered>" << p.second;
        }
        html << "</table>";
    }



    // Find clusters of ordinal offsets.
    class Cluster {
    public:
        int64_t firstOffset;
        int64_t lastOffset;
        uint64_t uniqueMarkerCount;
    };
    vector<Cluster> clusters;
    const uint64_t minMarkerCount = min(allMarkerKmerIds[0].size(), allMarkerKmerIds[1].size());
    const int64_t offsetDeltaTolerance = int64_t(std::round(driftRateTolerance * double(minMarkerCount)));
    for(uint64_t i=0; i<histogram.size(); /* Increment later */) {
        Cluster cluster;
        const uint64_t firstOffsetIndexInHistogram = i;
        cluster.firstOffset = histogram[firstOffsetIndexInHistogram].first;
        for(++i; i < histogram.size(); ++i) {
            if(histogram[i].first > histogram[i-1].first + offsetDeltaTolerance) {
                break;
            }
        }
        const uint64_t lastOffsetIndexInHistogram = i-1;
        cluster.lastOffset = histogram[lastOffsetIndexInHistogram].first;
        cluster.uniqueMarkerCount = 0;
        for(uint64_t j=firstOffsetIndexInHistogram; j<=lastOffsetIndexInHistogram; j++) {
            cluster.uniqueMarkerCount += histogram[j].second;
        }
        clusters.push_back(cluster);
    }

    // Find the largest cluster.
    uint64_t largestClusterIndex = invalid<uint64_t>;
    uint64_t largestClusterSize = 0;
    for(uint64_t i=0; i<clusters.size(); i++) {
        const uint64_t clusterSize = clusters[i].uniqueMarkerCount;
        if(clusterSize > largestClusterSize) {
            largestClusterSize = clusterSize;
            largestClusterIndex = i;
        }
    }
    const Cluster& largestCluster = clusters[largestClusterIndex];

    // Write the clusters.
    if(html) {
        html << "<h3>Ordinal offset clusters</h3>";
        html << "<p>Ordinal offset clusters were computed using offset tolerance " << offsetDeltaTolerance;
        html << "<table><tr><th>First<br>offset<th>Last<br>offset<th>Size";
        for(uint64_t i=0; i<clusters.size(); i++) {
            const Cluster& cluster = clusters[i];
            html << "<tr";
            if(i == largestClusterIndex) {
                html << " style='background-color:pink'";
            }
            html << ">"
                "<td class=centered>" << cluster.firstOffset <<
                "<td class=centered>" << cluster.lastOffset <<
                "<td class=centered>" << cluster.uniqueMarkerCount;
        }
        html << "</table>";
    }



    // Only keep the common unique markers on the largest cluster.
    {
        vector<CommonKmerInfo> newCommonKmerInfos;
        for(const CommonKmerInfo& commonKmerInfo: commonKmerInfos) {
            const int64_t offset = commonKmerInfo.ordinalOffset();
            if(offset >= largestCluster.firstOffset and offset <= largestCluster.lastOffset) {
                newCommonKmerInfos.push_back(commonKmerInfo);
            }
        }
        commonKmerInfos.swap(newCommonKmerInfos);
    }



    // Fill in the ordinal ranks.
    std::ranges::sort(commonKmerInfos, std::ranges::less(), &CommonKmerInfo::ordinal0);
    for(uint64_t rank=0; rank<commonKmerInfos.size(); rank++) {
        commonKmerInfos[rank].rank0 = rank;
    }
    std::ranges::sort(commonKmerInfos, std::ranges::less(), &CommonKmerInfo::ordinal1);
    for(uint64_t rank=0; rank<commonKmerInfos.size(); rank++) {
        commonKmerInfos[rank].rank1 = rank;
    }


    // If there are any markers that don't have the same rank, remove them.
    {
        vector<CommonKmerInfo> newCommonKmerInfos;
        for(const CommonKmerInfo& commonKmerInfo: commonKmerInfos) {
            if(commonKmerInfo.rank0 == commonKmerInfo.rank1) {
                newCommonKmerInfos.push_back(commonKmerInfo);
            }
        }
        commonKmerInfos.swap(newCommonKmerInfos);

    }



    // Sort them by ordinalSum.
    class OrderByOrdinalSum {
    public:
        bool operator()(const CommonKmerInfo& x, const CommonKmerInfo& y) const
        {
            return x.ordinalSum() < y.ordinalSum();
        }
    };
    sort(commonKmerInfos.begin(), commonKmerInfos.end(), OrderByOrdinalSum());



    // Write the common unique markers we kept.
    if(html) {
        html << "<h3>Active common unique markers</h3>";
        html << "There are " << commonKmerInfos.size() << " active common unique markers, "
            "shown in the table sorted by ordinal sum."
            "<p><table>"
            "<tr><th>Ordinal0<th>Ordinal1<th>Ordinal<br>offset<th>Ordinal<br>sum<th>Rank0<th>Rank1<th>KmerId<th>Kmer";
        const uint64_t k = assemblerInfo->k;
        for(const CommonKmerInfo& commonKmerInfo: commonKmerInfos) {
            const Kmer kmer(commonKmerInfo.kmerId, k);
            html << "<tr>"
                "<td class=centered>" << commonKmerInfo.ordinal0 <<
                "<td class=centered>" << commonKmerInfo.ordinal1 <<
                "<td class=centered>" << commonKmerInfo.ordinalOffset() <<
                "<td class=centered>" << commonKmerInfo.ordinalSum() <<
                "<td class=centered>" << commonKmerInfo.rank0 <<
                "<td class=centered>" << commonKmerInfo.rank1;

            // Write the KmerId in hex with the appropriate number of digits.
            const char oldFill = html.fill('0');
            html << "<td class=centered style='font-family:monospace'>" <<
                std::hex << std::setw(int(k/2)) << commonKmerInfo.kmerId << std::dec;
            html.fill(oldFill);

            // Write the Kmer.
            html << "<td class=centered style='font-family:monospace'>";
            kmer.write(html, k);
        }
        html << "</table>";
    }



    // We should remove common unique markers that have a different rank
    // in the two oriented reads. This does not happen frequently and
    // for now just check for them.
    for(const CommonKmerInfo& commonKmerInfo: commonKmerInfos) {
        SHASTA_ASSERT(commonKmerInfo.rank0 == commonKmerInfo.rank1);
    }


    if(commonKmerInfos.size() < 2) {
        alignment.clear();
        alignmentInfo.create(alignment, uint32_t(allMarkerKmerIds[0].size()), uint32_t(allMarkerKmerIds[1].size()));
        return;
    }



    // Create the alignment by stitching together alignments computed
    // between each pair of consecutive unique k-mers that survived
    // the above process (the "active" markers).
    alignment.clear();
    SHASTA_ASSERT(commonKmerInfos.size() > 1);

    // First, do an alignment between the beginning and the
    // first active unique marker.
    // This alignment is constrained on the right only.
    // We could do a bit better than this for performance.
    {
        const CommonKmerInfo& firstCommonKmerInfo = commonKmerInfos.front();
        const uint32_t ordinalB0 = firstCommonKmerInfo.ordinal0;
        const uint32_t ordinalB1 = firstCommonKmerInfo.ordinal1;
        if(ordinalB0 > 0 and ordinalB1 > 0) {
            const span<const KmerId> kmerIds0(&allMarkerKmerIds[0][0], &allMarkerKmerIds[0][ordinalB0]);
            const span<const KmerId> kmerIds1(&allMarkerKmerIds[1][0], &allMarkerKmerIds[1][ordinalB1]);
            vector< pair<bool, bool> > seqanAlignment;
            seqanAlign(
                kmerIds0.begin(), kmerIds0.end(),
                kmerIds1.begin(), kmerIds1.end(),
                matchScore, mismatchScore, gapScore,
                true, false,    // Free on left
                seqanAlignment);
            uint32_t ordinal0 = 0;
            uint32_t ordinal1 = 0;
            for(const auto& p: seqanAlignment) {
                if(p.first and p.second and allMarkerKmerIds[0][ordinal0] == allMarkerKmerIds[1][ordinal1]) {
                    alignment.ordinals.push_back({ordinal0, ordinal1});
                }
                if(p.first) {
                    ++ordinal0;
                }
                if(p.second) {
                    ++ordinal1;
                }
            }
            SHASTA_ASSERT(ordinal0 == ordinalB0);
            SHASTA_ASSERT(ordinal1 == ordinalB1);
        }
    }


    for(uint64_t step=1; step<commonKmerInfos.size(); step++) {
        const CommonKmerInfo& commonKmerInfoA = commonKmerInfos[step-1];
        const CommonKmerInfo& commonKmerInfoB = commonKmerInfos[step];
        SHASTA_ASSERT(commonKmerInfoB.rank0 > commonKmerInfoA.rank0);
        SHASTA_ASSERT(commonKmerInfoB.rank1 > commonKmerInfoA.rank1);

        const uint32_t ordinalA0 = commonKmerInfoA.ordinal0;
        const uint32_t ordinalA1 = commonKmerInfoA.ordinal1;
        const uint32_t ordinalB0 = commonKmerInfoB.ordinal0;
        const uint32_t ordinalB1 = commonKmerInfoB.ordinal1;

        // Get the KmerIds between A and B for the two reads.
        // These are the Kmers that we will align in this step.
        const span<const KmerId> kmerIds0(&allMarkerKmerIds[0][ordinalA0 + 1], &allMarkerKmerIds[0][ordinalB0]);
        const span<const KmerId> kmerIds1(&allMarkerKmerIds[1][ordinalA1 +1 ], &allMarkerKmerIds[1][ordinalB1]);
        if(html) {
            html << "<br>Step " << step << " alignment lengths " << kmerIds0.size() << " " << kmerIds1.size();
        }

        // Add to the alignment the first marker of this step.
        alignment.ordinals.push_back({commonKmerInfoA.ordinal0, commonKmerInfoA.ordinal1});

        // If there is nothing to align, we are done for this step,
        if(kmerIds0.empty() or kmerIds1.empty()) {
            continue;
        }

        // Use seqan to compute the alignment for this step.
        // This alignment is constrained on both sides.
        vector< pair<bool, bool> > seqanAlignment;
        const int64_t alignmentScore = seqanAlign(
            kmerIds0.begin(), kmerIds0.end(),
            kmerIds1.begin(), kmerIds1.end(),
            matchScore, mismatchScore, gapScore,
            false, false,
            seqanAlignment);
        if(html) {
            html << "<br>Alignment score " << alignmentScore;
        }

        // Add to the alignment the ordinals of matching alignment positions.
        uint32_t ordinal0 = ordinalA0 + 1;
        uint32_t ordinal1 = ordinalA1 + 1;
        for(const auto& p: seqanAlignment) {
            if(p.first and p.second and allMarkerKmerIds[0][ordinal0] == allMarkerKmerIds[1][ordinal1]) {
                alignment.ordinals.push_back({ordinal0, ordinal1});
            }
            if(p.first) {
                ++ordinal0;
            }
            if(p.second) {
                ++ordinal1;
            }
        }
        SHASTA_ASSERT(ordinal0 == ordinalB0);
        SHASTA_ASSERT(ordinal1 == ordinalB1);
    }

    // Add the last active marker.
    const CommonKmerInfo& lastCommonKmerInfo = commonKmerInfos.back();
    alignment.ordinals.push_back({lastCommonKmerInfo.ordinal0, lastCommonKmerInfo.ordinal1});



    // Do an alignment between the last active unique marker and the end.
    // This alignment is constrained on the left only.
    // We could do a bit better than this for performance.
    {
        const CommonKmerInfo& lastCommonKmerInfo = commonKmerInfos.back();
        const uint32_t ordinalA0 = lastCommonKmerInfo.ordinal0 + 1;
        const uint32_t ordinalA1 = lastCommonKmerInfo.ordinal1 + 1;
        const uint32_t ordinalB0 = uint32_t(allMarkerKmerIds[0].size());
        const uint32_t ordinalB1 = uint32_t(allMarkerKmerIds[1].size());
        if( ordinalA0 < ordinalB0 and ordinalA1 < ordinalB1) {
            const span<const KmerId> kmerIds0(&allMarkerKmerIds[0][ordinalA0], &allMarkerKmerIds[0][ordinalB0]);
            const span<const KmerId> kmerIds1(&allMarkerKmerIds[1][ordinalA1], &allMarkerKmerIds[1][ordinalB1]);
            vector< pair<bool, bool> > seqanAlignment;
            seqanAlign(
                kmerIds0.begin(), kmerIds0.end(),
                kmerIds1.begin(), kmerIds1.end(),
                matchScore, mismatchScore, gapScore,
                false, true,    // Free on right
                seqanAlignment);
            uint32_t ordinal0 = ordinalA0;
            uint32_t ordinal1 = ordinalA1;
            for(const auto& p: seqanAlignment) {
                if(p.first and p.second and allMarkerKmerIds[0][ordinal0] == allMarkerKmerIds[1][ordinal1]) {
                    alignment.ordinals.push_back({ordinal0, ordinal1});
                }
                if(p.first) {
                    ++ordinal0;
                }
                if(p.second) {
                    ++ordinal1;
                }
            }
            SHASTA_ASSERT(ordinal0 == ordinalB0);
            SHASTA_ASSERT(ordinal1 == ordinalB1);
        }
    }

    // Store the alignment info.
    alignmentInfo.create(alignment, uint32_t(allMarkerKmerIds[0].size()), uint32_t(allMarkerKmerIds[1].size()));

}



void Assembler::computeLowFrequencyMarkers(
    uint64_t maxMarkerFrequency,
    uint64_t threadCount)
{
    // Check that we have what we need.
    SHASTA_ASSERT(markerKmerIds.isOpen());

    // Get the number of reads.
    const uint64_t readCount = getReads().readCount();

    // Adjust the number of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Store the maxMarkerFrequency so all threads can see it.
    computeLowFrequencyMarkersData.maxMarkerFrequency = maxMarkerFrequency;

    // Initialize the low frequency markers.
    lowFrequencyMarkers.createNew(largeDataName("LowFrequencyMarkers"), largeDataPageSize);

    // Pass 1 just counts the number of low frequency markers for each oriented read.
    const uint64_t batchSize = 1;
    lowFrequencyMarkers.beginPass1(2 * readCount);
    setupLoadBalancing(readCount, batchSize);
    runThreads(&Assembler::computeLowFrequencyMarkersThreadFunctionPass1, threadCount);

    // Pass 2 stores the low frequency markers for each oriented read.
    setupLoadBalancing(getReads().readCount(), batchSize);
    lowFrequencyMarkers.beginPass2();
    runThreads(&Assembler::computeLowFrequencyMarkersThreadFunctionPass2, threadCount);
    lowFrequencyMarkers.endPass2(false, true);
}



void Assembler::computeLowFrequencyMarkersThreadFunctionPass1(uint64_t threadId)
{
    computeLowFrequencyMarkersThreadFunctionPass12(1);
}
void Assembler::computeLowFrequencyMarkersThreadFunctionPass2(uint64_t threadId)
{
    computeLowFrequencyMarkersThreadFunctionPass12(2);
}
void Assembler::computeLowFrequencyMarkersThreadFunctionPass12(uint64_t pass)
{
    const uint64_t maxMarkerFrequency = computeLowFrequencyMarkersData.maxMarkerFrequency;
    vector<uint32_t> lowFrequencyOrdinals;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over oriented reads in this batch.
        for(uint32_t readId=ReadId(begin); readId!=ReadId(end); ++readId) {
            for(uint32_t strand=0; strand<2; strand++) {
                const OrientedReadId orientedReadId(readId, strand);

                // Compute the low frequency markers.
                computeLowFrequencyMarkers(
                    markerKmerIds[orientedReadId.getValue()],
                    maxMarkerFrequency,
                    lowFrequencyOrdinals);

                if(pass == 1) {
                    // Just make space for them.
                    lowFrequencyMarkers.incrementCountMultithreaded(
                        orientedReadId.getValue(),
                        lowFrequencyOrdinals.size());
                } else {
                    // Store them.
                    copy(lowFrequencyOrdinals.begin(), lowFrequencyOrdinals.end(),
                        lowFrequencyMarkers.begin(orientedReadId.getValue()));
                }
            }
        }
    }
}



// Compute low frequency markers for a single oriented read.
// On return, the lowFrequencyOrdinals vector contains the ordinals corresponding
// to low frequency markers, sorted by KmerId.
// Low frequency markers are the ones that occur up to maxMarkerFrequency
// times on the oriented read.
void Assembler::computeLowFrequencyMarkers(
    const span<const KmerId>& kmerIds,          // The marker KmerIds for the oriented reads, sorted by ordinal
    uint64_t maxMarkerFrequency,
    vector<uint32_t>& lowFrequencyOrdinals)     // The ordinals of the low frequency markers, sorted by KmerId
{

    // Create a vector of ordinals, sorted by ordinal.
    const uint64_t markerCount = kmerIds.size();
    vector<uint32_t> allOrdinals(markerCount);
    std::iota(allOrdinals.begin(), allOrdinals.end(), uint32_t(0));

    // Now sort them by KmerId.
    class SortHelper {
    public:
        SortHelper(const span<const KmerId>& kmerIds) : kmerIds(kmerIds) {}
        bool operator()(uint32_t ordinal0, uint32_t ordinal1) const
        {
            return kmerIds[ordinal0] < kmerIds[ordinal1];
        }
    private:
        const span<const KmerId>& kmerIds;
    };
    sort(allOrdinals.begin(), allOrdinals.end(), SortHelper(kmerIds));



    // Loop over streaks with the same KmerId.
    lowFrequencyOrdinals.clear();
    for(uint64_t streakBegin=0; streakBegin<markerCount; /* Increment later */) {
        const KmerId kmerId = kmerIds[allOrdinals[streakBegin]];

        // Find the streak with this KmerId.
        uint64_t streakEnd = streakBegin + 1;
        while(true) {
            if(streakEnd == markerCount or kmerIds[allOrdinals[streakEnd]] != kmerId) {
                break;
            }
            ++streakEnd;
        }
        const uint64_t streakLength = streakEnd - streakBegin;

        // If short enough, copy to the low frequency ordinals.
        if(streakLength <= maxMarkerFrequency) {
            copy(allOrdinals.begin() + streakBegin, allOrdinals.begin() + streakEnd,
                back_inserter(lowFrequencyOrdinals));
        }

        // Prepare to process the next streak.
        streakBegin = streakEnd;
    }
}
