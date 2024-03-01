#include "Assembler.hpp"
#include "deduplicate.hpp"
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



    // Gather the unique markers in the two oriented reads.
    // For performance, this should done for all oriented reads
    // before alignment computation begins.
    class KmerInfo {
    public:
        uint32_t ordinal;
        KmerId kmerId;
        // Comparison by KmerId only.
        // This is used below to find unique markers.
        bool operator<(const KmerInfo& that) const
        {
            return kmerId < that.kmerId;
        }
        bool operator==(const KmerInfo& that) const
        {
            return kmerId == that.kmerId;
        }
    };
    array<vector<KmerInfo>, 2> kmerInfos;

    // Loop over the two oriented reads.
    for(uint64_t i=0; i<2; i++) {

        // Gather the markers for this oriented read.
        const uint64_t markerCount = allMarkerKmerIds[i].size();
        kmerInfos[i].resize(markerCount);
        for(uint32_t ordinal=0; ordinal<markerCount; ordinal++) {
            kmerInfos[i][ordinal] = {ordinal, allMarkerKmerIds[i][ordinal]};
        }

        // Deduplicate and count occurrences.
        vector<uint64_t> frequency;
        deduplicateAndCount(kmerInfos[i], frequency);

        // Only keep the unique ones.
        vector<KmerInfo> uniqueKmerInfos;
        for(uint64_t j=0; j<frequency.size(); j++) {
            if(frequency[j] == 1) {
                uniqueKmerInfos.push_back(kmerInfos[i][j]);
            }
        }
        kmerInfos[i].swap(uniqueKmerInfos);

        if(html) {
            html << "<br>" << (i==0 ? orientedReadId0 : orientedReadId1) << " has " << markerCount <<
                " markers of which " << kmerInfos[i].size() << " unique." << endl;
        }
    }



    // Find markers with a KmerId that is unique in both oriented reads
    // and occurs in both.
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
    const auto begin0 = kmerInfos[0].begin();
    const auto begin1 = kmerInfos[1].begin();
    const auto end0 = kmerInfos[0].end();
    const auto end1 = kmerInfos[1].end();
    auto it0 = begin0;
    auto it1 = begin1;
    while((it0 != end0) and (it1 != end1)) {
        const KmerId kmerId0 = it0->kmerId;
        const KmerId kmerId1 = it1->kmerId;
        if(kmerId0 < kmerId1) {
            ++it0;
        } else if(kmerId1 < kmerId0) {
            ++it1;
        } else {
            SHASTA_ASSERT(kmerId0 == kmerId1);
            commonKmerInfos.push_back({it0->ordinal, it1->ordinal, kmerId0});
            ++it0;
            ++it1;
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


#if 0
    // Recreate the kmerInfos vectors, this time including only the unique markers.
    for(uint64_t i=0; i<2; i++) {
        kmerInfos[i].clear();
    }
    for(const CommonKmerInfo& commonKmerInfo: commonKmerInfos) {
        kmerInfos[0].push_back({commonKmerInfo.ordinal0, commonKmerInfo.kmerId});
        kmerInfos[1].push_back({commonKmerInfo.ordinal1, commonKmerInfo.kmerId});
    }

    // Sort them by ordinal.
    class OrderKmerInfosByOrdinal {
    public:
        bool operator()(const KmerInfo& x, const KmerInfo& y) const
        {
            return x.ordinal < y.ordinal;
        }
    };
    for(uint64_t i=0; i<2; i++) {
        sort(kmerInfos[i].begin(), kmerInfos[i].end(), OrderKmerInfosByOrdinal());
    }

    if(html) {
        html << "<table>";
        for(uint64_t j0=0; j0<kmerInfos[0].size(); j0++) {
            const KmerInfo& kmerInfo0 = kmerInfos[0][j0];
            for(uint64_t j1=0; j1<kmerInfos[1].size(); j1++) {
                const KmerInfo& kmerInfo1 = kmerInfos[1][j1];
                if(kmerInfo0.kmerId == kmerInfo1.kmerId) {
                    html << "<tr>"
                        "<td class=centered>" << j0 <<
                        "<td class=centered>" << j1;
                }
            }
        }
        html << "</table>";
    }



    // We will align the sequences of these unique k-mers.
    array<vector<KmerId>, 2> compressedSequences;
    for(uint64_t i=0; i<2; i++) {
        for(const KmerInfo& kmerInfo: kmerInfos[i]) {
            compressedSequences[i].push_back(kmerInfo.kmerId);
        }
    }
    vector< pair<bool, bool> > booleanAlignment;
    seqanAlign(
        compressedSequences[0].begin(), compressedSequences[0].end(),
        compressedSequences[1].begin(), compressedSequences[1].end(),
        matchScore, mismatchScore, gapScore,
        true, true, booleanAlignment);

    if(html) {
        html << "<br>Alignment of common unique markers follows:<pre>";
        uint64_t position0 = 0;
        uint64_t position1 = 0;
        for(const auto& p: booleanAlignment) {
            if(p.first and p.second) {
                if(compressedSequences[0][position0] == compressedSequences[1][position1]) {
                    html << ".";
                } else {
                    html << "X";
                }
            } else if (p.first){
                SHASTA_ASSERT(not p.second);
                html << ".";
            } else {
                SHASTA_ASSERT(p.second);
                html << "-";
            }
            if(p.first) {
                ++position0;
            }
            if(p.second) {
                ++position1;
            }
        }
        html << "\n";
        position0 = 0;
        position1 = 0;
        for(const auto& p: booleanAlignment) {
            if(p.first and p.second) {
                if(compressedSequences[0][position0] == compressedSequences[1][position1]) {
                    html << ".";
                } else {
                    html << "X";
                }
            } else if (p.second){
                SHASTA_ASSERT(not p.first);
                html << ".";
            } else {
                SHASTA_ASSERT(p.first);
                html << "-";
            }
            if(p.first) {
                ++position0;
            }
            if(p.second) {
                ++position1;
            }
        }
        html << "</pre>";
    }
#endif
}
