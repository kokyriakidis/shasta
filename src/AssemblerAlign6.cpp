#include "Assembler.hpp"
#include "KmerCounter.hpp"
#include "orderPairs.hpp"
using namespace shasta;



// Version that uses banded alignments.
void Assembler::alignOrientedReads6(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    int matchScore,
    int mismatchScore,
    int gapScore,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    ostream& html)
{
    // EXPOSE WHEN CODE STABILIZES
    const uint64_t maxLocalFrequency = 1000000000;
    const uint64_t minGlobalFrequency = 10;
    const uint64_t maxGlobalFrequency = 50;

    SHASTA_ASSERT(kmerCounter and kmerCounter->isAvailable());

    // Get the marker KmerIds for the two oriented reads.
    array<span<KmerId>, 2> orientedReadKmerIds;
    array<vector<KmerId>, 2> orientedReadKmerIdsVectors;
    if(markerKmerIds.isOpen()) {
        if(html) {
            html << "<br>Using stored marker KmerIds." << endl;
        }

        orientedReadKmerIds[0] = markerKmerIds[orientedReadId0.getValue()];
        orientedReadKmerIds[1] = markerKmerIds[orientedReadId1.getValue()];

    } else {
        if(html) {
            html << "<br>Stored sorted marker KmerIds are not available - computing them." << endl;
        }

        // This is slower and will happen if markerKmerIds is not available.
        // Resize the vectors and make the spans point to the vectors.
        // Then call getOrientedReadMarkerKmerIds to fill them in.
        orientedReadKmerIdsVectors[0].resize(markers.size(orientedReadId0.getValue()));
        orientedReadKmerIdsVectors[1].resize(markers.size(orientedReadId1.getValue()));
        orientedReadKmerIds[0] = span<KmerId>(orientedReadKmerIdsVectors[0]);
        orientedReadKmerIds[1] = span<KmerId>(orientedReadKmerIdsVectors[1]);
        getOrientedReadMarkerKmerIds(orientedReadId0, orientedReadKmerIds[0]);
        getOrientedReadMarkerKmerIds(orientedReadId1, orientedReadKmerIds[1]);
    }



    // Align6 needs markers sorted by KmerId.
    // Use the ones from sortedMarkers if available, or else compute them.
    array<span< pair<KmerId, uint32_t> >, 2> orientedReadSortedMarkersSpans;
    array<vector< pair<KmerId, uint32_t> >, 2> orientedReadSortedMarkers;
    if(sortedMarkers.isOpen()) {
        if(html) {
            html << "<br>Using stored sorted markers." << endl;
        }

        // Make the spans point to the stored sorted markers.
        orientedReadSortedMarkersSpans[0] = sortedMarkers[orientedReadId0.getValue()];
        orientedReadSortedMarkersSpans[1] = sortedMarkers[orientedReadId1.getValue()];

    } else {

        // We don't have the sorted markers. we have to compute them here.
        if(html) {
            html << "<br>Stored sorted markers are not available - computing them." << endl;
        }

        for(uint64_t i=0; i<2; i++) {

            // Unsorted markers for this oriented read.
            const span<const KmerId>& km = orientedReadKmerIds[i];

            // Sorted markers for this oriented read.
            vector<pair<KmerId, uint32_t> >& sm = orientedReadSortedMarkers[i];

            // Copy the unsorted markers.
            const uint64_t n = km.size();
            sm.resize(n);
            for(uint64_t ordinal=0; ordinal<n; ordinal++) {
                sm[ordinal] = make_pair(km[ordinal], uint32_t(ordinal));
            }

            // Sort them.
            sort(sm.begin(), sm.end(), OrderPairsByFirstOnly<KmerId, uint32_t>());

            // Make the span point to the data in the vector.
            orientedReadSortedMarkersSpans[i] = sm;
        }
    }



    // Find pairs of markers in the two oriented reads that have the same k-mer
    // and satisfy the frequency requirements.
    class MarkerPairInfo {
    public:
        uint32_t ordinal0;
        uint32_t ordinal1;
        KmerId kmerId;
        uint64_t globalFrequency;
        uint64_t localFrequency0;
        uint64_t localFrequency1;
        uint64_t ordinalSum() const
        {
            return ordinal0 + ordinal1;
        }
        int64_t ordinalOffset() const
        {
            return int64_t(ordinal0) - int64_t(ordinal1);
        }
    };
    vector<MarkerPairInfo> markerPairInfos;



    if(html) {
        html <<
            "<h3>Marker pairs</h2>"
            "<table><tr>"
            "<th>Ordinal0"
            "<th>Ordinal1"
            "<th>Ordinal<br>sum"
            "<th>Ordinal<br>offset"
            "<th>Kmer"
            "<th>KmerId"
            "<th>Local<br>frequency0"
            "<th>Local<br>frequency1"
            "<th>Global<br>frequency";
    }


    // Joint loop over markers sorted by KmerId to find the non-zero elements of
    // the alignment matrix. That is, we want to find pairs (ordinal0, ordinal1)
    // such that KmerId(orientedReadId0, ordinal0) == KmerId(orientedReadId1, ordinal1).
    // Joint loop over the sorted markers, looking for common markers.
    auto begin0 = orientedReadSortedMarkersSpans[0].begin();
    auto begin1 = orientedReadSortedMarkersSpans[1].begin();
    auto end0 = orientedReadSortedMarkersSpans[0].end();
    auto end1 = orientedReadSortedMarkersSpans[1].end();

    auto it0 = begin0;
    auto it1 = begin1;
    while(it0!=end0 && it1!=end1) {
        if(it0->first < it1->first) {
            ++it0;
        } else if(it1->first < it0->first) {
            ++it1;
        } else {

            // We found a common KmerId.
            const KmerId kmerId = it0->first;

            // This KmerId could appear more than once in each of the sequences,
            // so we need to find the streak of this KmerId.
            auto it0Begin = it0;
            auto it1Begin = it1;
            auto it0End = it0Begin;
            auto it1End = it1Begin;
            while(it0End!=end0 && it0End->first == kmerId) {
                ++it0End;
            }
            while(it1End!=end1 && it1End->first == kmerId) {
                ++it1End;
            }

            // If the local global frequencies are in the desired range,
            // loop over pairs in the streaks.
            const uint64_t localFrequency0 = it0End - it0Begin;
            const uint64_t localFrequency1 = it1End - it1Begin;
            const uint64_t globalFrequency = kmerCounter->getFrequency(kmerId);
            if(
                (localFrequency0 <= maxLocalFrequency) and
                (localFrequency1 <= maxLocalFrequency) and
                (globalFrequency >= minGlobalFrequency) and
                (globalFrequency <= maxGlobalFrequency)) {

                MarkerPairInfo markerPairInfo;
                markerPairInfo.kmerId = kmerId;
                markerPairInfo. globalFrequency = globalFrequency;
                markerPairInfo.localFrequency0 = localFrequency0;
                markerPairInfo.localFrequency1 = localFrequency1;

                for(auto jt0=it0Begin; jt0!=it0End; ++jt0) {
                    markerPairInfo. ordinal0 = jt0->second;
                    for(auto jt1=it1Begin; jt1!=it1End; ++jt1) {
                        markerPairInfo.ordinal1 = jt1->second;
                        markerPairInfos.push_back(markerPairInfo);

                        if(html) {
                            html <<
                                "<tr>"
                                "<td class=centered>" << markerPairInfo.ordinal0 <<
                                "<td class=centered>" << markerPairInfo.ordinal1 <<
                                "<td class=centered>" << markerPairInfo.ordinalSum() <<
                                "<td class=centered>" << markerPairInfo.ordinalOffset() <<
                                "<td class=centered style='font-family:courier'>";
                            Kmer(kmerId, assemblerInfo->k).write(html, assemblerInfo->k);

                            html <<
                                "<td class=centered>" << kmerId <<
                                "<td class=centered>" << localFrequency0 <<
                                "<td class=centered>" << localFrequency1 <<
                                "<td class=centered>" << globalFrequency;
                        }
                    }
                }
            }

            // Continue the joint loop over KmerId's.
            it0 = it0End;
            it1 = it1End;
        }
    }
    if(html) {
        html << "</table>";
    }



    // Create a histogram of ordinal offsets for the active common markers.
    std::map<int64_t, uint64_t> histogramMap;
    for(const MarkerPairInfo& markerPairInfo: markerPairInfos) {
        const int64_t offset = markerPairInfo.ordinalOffset();
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
        html << "<h3>Histogram of ordinal offsets for the marker pairs</h3>"
            "<table>"
            "<tr><th>Ordinal<br>offset<th>Gap to<br>previous<th>Frequency";
        int64_t previousOffset = invalid<int64_t>;
        for(const auto& p: histogram) {
            const int64_t offset = p.first;
            html << "<tr><td class=centered>" << offset << "<td class=centered>";
            if(previousOffset != invalid<int64_t>) {
                html << offset - previousOffset - 1;
            }
            html << "<td class=centered>" << p.second;
            previousOffset = offset;
        }
        html << "</table>";
    }


    SHASTA_ASSERT(0);
}
