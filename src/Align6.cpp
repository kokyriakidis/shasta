// Shasta.
#include "Align6.hpp"
#include "Alignment.hpp"
#include "KmerCounter.hpp"
#include "longestPath.hpp"
#include "orderPairs.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include "fstream.hpp"
#include "vector.hpp"



Align6::Align6(
    const array<span< pair<KmerId, uint32_t> >, 2>& orientedReadSortedMarkersSpans,
    uint64_t k,
    const KmerCounter& kmerCounter,
    uint64_t maxSkip,
    uint64_t maxDrift,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    ostream& html)
{
    // EXPOSE WHEN CODE STABILIZES
    const uint64_t maxLocalFrequency = 1000000000;
    const uint64_t minGlobalFrequency = 10;
    const uint64_t maxGlobalFrequency = 50;
    const double driftRateTolerance = 0.05;

    // Sanity check.
    SHASTA_ASSERT(kmerCounter.isAvailable());


    // Class to store information about a pair of markers in the
    // two oriented reads that have the same KmerId.
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

        // Sort by ordinalSum.
        bool operator<(const MarkerPairInfo& that) const
        {
            return ordinalSum() < that.ordinalSum();
        }
    };



    // Joint loop over markers sorted by KmerId to find the low frequency
    // marker pairs. That is, we want to find pairs (ordinal0, ordinal1)
    // such that KmerId(orientedReadId0, ordinal0) == KmerId(orientedReadId1, ordinal1).
    // Do a joint loop over the sorted markers, looking for common markers.
    vector<MarkerPairInfo> lowFrequencyMarkerPairInfos;
    const auto begin0 = orientedReadSortedMarkersSpans[0].begin();
    const auto begin1 = orientedReadSortedMarkersSpans[1].begin();
    const auto end0 = orientedReadSortedMarkersSpans[0].end();
    const auto end1 = orientedReadSortedMarkersSpans[1].end();

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
            const uint64_t globalFrequency = kmerCounter.getFrequency(kmerId);
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
                        lowFrequencyMarkerPairInfos.push_back(markerPairInfo);

                    }
                }
            }

            // Continue the joint loop over KmerId's.
            it0 = it0End;
            it1 = it1End;
        }
    }



    if(lowFrequencyMarkerPairInfos.empty()) {
        if(html) {
            html << "<br>No low frequency marker pairs found.";
        }
        return;
    }


    // Write out the low frequency marker pairs.
    if(html) {
        html <<
            "<h3>Low frequency marker pairs</h2>"
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

        for(const MarkerPairInfo& markerPairInfo: lowFrequencyMarkerPairInfos) {
            html <<
                "<tr>"
                "<td class=centered>" << markerPairInfo.ordinal0 <<
                "<td class=centered>" << markerPairInfo.ordinal1 <<
                "<td class=centered>" << markerPairInfo.ordinalSum() <<
                "<td class=centered>" << markerPairInfo.ordinalOffset() <<
                "<td class=centered style='font-family:courier'>";

            Kmer(markerPairInfo.kmerId, k).write(html, k);

            html <<
                "<td class=centered>" << markerPairInfo.kmerId <<
                "<td class=centered>" << markerPairInfo.localFrequency0 <<
                "<td class=centered>" << markerPairInfo.localFrequency1 <<
                "<td class=centered>" << markerPairInfo.globalFrequency;
        }

        html << "</table>";
    }



    // Create a histogram of ordinal offsets for the low frequency marker pairs.
    std::map<int64_t, uint64_t> histogramMap;
    for(const MarkerPairInfo& markerPairInfo: lowFrequencyMarkerPairInfos) {
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



    // Write out the histogram.
    if(html) {
        html << "<h3>Histogram of ordinal offsets for the low frequency marker pairs</h3>"
            "<table>"
            "<tr><th>Ordinal<br>offset<th>Frequency";
        for(const auto& p: histogram) {
            const int64_t offset = p.first;
            const uint64_t frequency = p.second;

            html << "<tr><td class=centered>" << offset;
            html << "<td class=centered>" << frequency;
        }
        html << "</table>";
    }



    // Create a convolution kernel that will be used to get a smoothed version
    // of the histogram.
    const uint64_t minLength = min(
        orientedReadSortedMarkersSpans[0].size(),
        orientedReadSortedMarkersSpans[1].size());
    const int64_t kernelWidth = int64_t(std::round(driftRateTolerance * double(minLength)));
    vector<double> kernel(2 * kernelWidth + 1);
    for(uint64_t i=0; i<kernel.size(); i++) {
        const double x = double(int64_t(i) - kernelWidth) / double(kernelWidth);
        kernel[i] = 1. - fabs(x);
    }



    // Write the kernel.
    if(html) {
        html <<
            "<h3>Convolution kernel</h3>"
            "<table>"
            "<tr><th>Delta<th>Kernel";
        for(uint64_t i=0; i<kernel.size(); i++) {
            html <<
                "<tr><td class=centered>" << int64_t(i) - kernelWidth <<
                "<td class=centered>" << kernel[i];
        }
        html << "</table>";
    }



    // Compute the convolution.
    const int64_t minHistogramOffset = histogram.front().first;
    const int64_t maxHistogramOffset = histogram.back().first;
    const int64_t minDistributionOffset = minHistogramOffset - kernelWidth;
    const int64_t maxDistributionOffset = maxHistogramOffset + kernelWidth;
    vector<double> distribution(maxDistributionOffset - minDistributionOffset + 1, 0.);
    for(const auto& p: histogram) {
        const int64_t offset = p.first;
        const uint64_t frequency = p.second;
        const double weight = double(frequency);
        for(uint64_t i=0; i<kernel.size(); i++) {
            const int64_t j = offset + int64_t(i) - kernelWidth - minDistributionOffset;
            SHASTA_ASSERT(j >= 0);
            SHASTA_ASSERT(j < int64_t(distribution.size()));
            distribution[j] += weight * kernel[i];
        }
    }



    // Write the smoothed offset distribution.
    if(html) {
        html <<
            "<h3>Smoothed offset distribution</h3>"
            "<table>"
            "<tr><th>Offset<th>Value";
        for(uint64_t i=0; i<distribution.size(); i++) {
            const int64_t offset = int64_t(i) + minDistributionOffset;
            html << "<tr><td class=centered>" << offset <<
                "<td class=centered>" << distribution[i];
        }
        html << "</table>";
    }



    // The tentative band is centered at the peak of the smoothed offset distribution.
    const int64_t tentativeBandCenter =
        int64_t(std::max_element(distribution.begin(), distribution.end()) - distribution.begin()) + minDistributionOffset;
    const int64_t tentativeBandLow = tentativeBandCenter - kernelWidth;
    const int64_t tentativeBandHigh = tentativeBandCenter + kernelWidth;
    if(html) {
        html << "<h3>Tentative band</h3>"
            "Tentative band includes offset " << tentativeBandLow << " through " << tentativeBandHigh <<
            "<br>Tentative band is centered at " << tentativeBandCenter;
    }



    // Gather all marker pairs within this tentative band.
    vector<MarkerPairInfo> inBandMarkerPairInfos;
    it0 = begin0;
    it1 = begin1;
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

            // Loop over pairs in the streaks.
            const uint64_t localFrequency0 = it0End - it0Begin;
            const uint64_t localFrequency1 = it1End - it1Begin;
            const uint64_t globalFrequency = kmerCounter.getFrequency(kmerId);

            MarkerPairInfo markerPairInfo;
            markerPairInfo.kmerId = kmerId;
            markerPairInfo. globalFrequency = globalFrequency;
            markerPairInfo.localFrequency0 = localFrequency0;
            markerPairInfo.localFrequency1 = localFrequency1;

            for(auto jt0=it0Begin; jt0!=it0End; ++jt0) {
                markerPairInfo. ordinal0 = jt0->second;
                for(auto jt1=it1Begin; jt1!=it1End; ++jt1) {
                    markerPairInfo.ordinal1 = jt1->second;
                    const int64_t offset = markerPairInfo.ordinalOffset();
                    if((offset >= tentativeBandLow) and (offset <= tentativeBandHigh)) {
                        inBandMarkerPairInfos.push_back(markerPairInfo);
                    }
                }
            }

            // Continue the joint loop over KmerId's.
            it0 = it0End;
            it1 = it1End;
        }
    }

    // Sort them by ordinalSum.
    sort(inBandMarkerPairInfos.begin(), inBandMarkerPairInfos.end());



    // Create connected components of a graph in which each vertex represents a MarkerPairInfo in
    // the inBandMarkerPairInfos vector, and a directed edge is created
    // if two MarkerPairInfo are compatible with maxSkip, maxDrift.
    const uint64_t maxOffsetSumDelta = 2 * maxSkip;
    vector<uint64_t> rank(inBandMarkerPairInfos.size());
    vector<uint64_t> parent(inBandMarkerPairInfos.size());
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t v=0; v<inBandMarkerPairInfos.size(); v++) {
        disjointSets.make_set(v);
    }
    for(uint64_t v0=0; v0<inBandMarkerPairInfos.size(); v0++) {
        const MarkerPairInfo& markerPairInfo0 = inBandMarkerPairInfos[v0];
        const uint64_t ordinalSum0 = markerPairInfo0.ordinalSum();
        for(uint64_t v1=v0+1; v1<inBandMarkerPairInfos.size(); v1++) {
            const MarkerPairInfo& markerPairInfo1 = inBandMarkerPairInfos[v1];
            const uint64_t ordinalSum1 = markerPairInfo1.ordinalSum();
            if(ordinalSum1 - ordinalSum0 > maxOffsetSumDelta) {
                break;
            }

            // Check skip.
            const uint64_t skip = max(
                abs(int32_t(markerPairInfo0.ordinal0) - int32_t(markerPairInfo1.ordinal0)),
                abs(int32_t(markerPairInfo0.ordinal1) - int32_t(markerPairInfo1.ordinal1)));
            if(skip > maxSkip) {
                continue;
            }

            // Check drift.
            const int64_t offset0 =  markerPairInfo0.ordinalOffset();
            const int64_t offset1 =  markerPairInfo1.ordinalOffset();
            const uint64_t drift = labs(offset0 - offset1);
            if(drift > maxDrift) {
                continue;
            }

            // Update the disjoint sets to take this edge into account.
            disjointSets.union_set(v0, v1);
        }
    }

    // Compute connected components.
    vector<uint64_t> component(inBandMarkerPairInfos.size());
    for(uint64_t v=0; v<inBandMarkerPairInfos.size(); v++) {
        component[v] = disjointSets.find_set(v);
    }



    // Write out the marker pairs in the tentative band.
    if(html) {
        html <<
            "<h3>Marker pairs in the tentative band</h3>"
            "<table><tr>"
            "<th>Id"
            "<th>Ordinal0"
            "<th>Ordinal1"
            "<th>Ordinal<br>sum"
            "<th>Ordinal<br>offset"
            "<th>Kmer"
            "<th>KmerId"
            "<th>Local<br>frequency0"
            "<th>Local<br>frequency1"
            "<th>Global<br>frequency"
            "<th>Component";

        for(uint64_t i=0; i<inBandMarkerPairInfos.size(); i++) {
            const MarkerPairInfo& markerPairInfo = inBandMarkerPairInfos[i];
            html <<
                "<tr>"
                "<td class=centered>" << i <<
                "<td class=centered>" << markerPairInfo.ordinal0 <<
                "<td class=centered>" << markerPairInfo.ordinal1 <<
                "<td class=centered>" << markerPairInfo.ordinalSum() <<
                "<td class=centered>" << markerPairInfo.ordinalOffset() <<
                "<td class=centered style='font-family:courier'>";

            Kmer(markerPairInfo.kmerId, k).write(html, k);

            html <<
                "<td class=centered>" << markerPairInfo.kmerId <<
                "<td class=centered>" << markerPairInfo.localFrequency0 <<
                "<td class=centered>" << markerPairInfo.localFrequency1 <<
                "<td class=centered>" << markerPairInfo.globalFrequency <<
                "<td class=centered>" << component[i];
        }

        html << "</table>";
    }



    // Count low frequency marker pairs in each connected component.
    vector<uint64_t> count(inBandMarkerPairInfos.size(), 0);
    for(uint64_t i=0; i<inBandMarkerPairInfos.size(); i++) {
        const MarkerPairInfo& markerPairInfo = inBandMarkerPairInfos[i];
        if(markerPairInfo.globalFrequency > maxGlobalFrequency) {
            continue;
        }
        if(markerPairInfo.globalFrequency < minGlobalFrequency) {
            continue;
        }
        if(markerPairInfo.localFrequency0 > maxLocalFrequency) {
            continue;
        }
        if(markerPairInfo.localFrequency1 > maxLocalFrequency) {
            continue;
        }
        const uint64_t c = component[i];
        ++count[c];
    }

    if(html) {
        html << "<h3>Number of low frequency marker pairs in each component</h3>"
            "<table><tr><th>Component<th>Low<br>frequency<br>markers";
        for(uint64_t c=0; c<inBandMarkerPairInfos.size(); c++) {
            if(count[c] > 0) {
                html << "<tr><td class=centered>" << c << "<td class=centered>" << count[c];
            }
        }
        html << "</table>";
    }


    // Find the connected component with the most low frequency markers.
    const uint64_t bestComponent = max_element(count.begin(), count.end()) - count.begin();
    if(html) {
        html << "<br>Will use marker pairs from component " << bestComponent;
    }


    // Gather the marker pairs in this component.
    vector<MarkerPairInfo> activeMarkerPairs;
    for(uint64_t v=0; v<inBandMarkerPairInfos.size(); v++) {
        if(component[v] == bestComponent) {
            activeMarkerPairs.push_back(inBandMarkerPairInfos[v]);
        }
    }



    // Write out the active marker pairs.
    if(html) {
        html <<
            "<h3>Active marker pairs</h3>"
            "<table><tr>"
            "<th>Id"
            "<th>Ordinal0"
            "<th>Ordinal1"
            "<th>Ordinal<br>sum"
            "<th>Ordinal<br>offset"
            "<th>Kmer"
            "<th>KmerId"
            "<th>Local<br>frequency0"
            "<th>Local<br>frequency1"
            "<th>Global<br>frequency";

        for(uint64_t i=0; i<activeMarkerPairs.size(); i++) {
            const MarkerPairInfo& markerPairInfo = activeMarkerPairs[i];
            html <<
                "<tr>"
                "<td class=centered>" << i <<
                "<td class=centered>" << markerPairInfo.ordinal0 <<
                "<td class=centered>" << markerPairInfo.ordinal1 <<
                "<td class=centered>" << markerPairInfo.ordinalSum() <<
                "<td class=centered>" << markerPairInfo.ordinalOffset() <<
                "<td class=centered style='font-family:courier'>";

            Kmer(markerPairInfo.kmerId, k).write(html, k);

            html <<
                "<td class=centered>" << markerPairInfo.kmerId <<
                "<td class=centered>" << markerPairInfo.localFrequency0 <<
                "<td class=centered>" << markerPairInfo.localFrequency1 <<
                "<td class=centered>" << markerPairInfo.globalFrequency;
        }

        html << "</table>";
    }



    // Use the active marker pairs to create a directed graph
    // whose edges correspond to MarkerPairInfos compatible with maxSkip, maxDrift.
    using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS>;
    Graph graph(activeMarkerPairs.size());
    for(uint64_t v0=0; v0<activeMarkerPairs.size(); v0++) {
        const MarkerPairInfo& markerPairInfo0 = activeMarkerPairs[v0];
        const uint64_t ordinalSum0 = markerPairInfo0.ordinalSum();
        for(uint64_t v1=v0+1; v1<activeMarkerPairs.size(); v1++) {
            const MarkerPairInfo& markerPairInfo1 = activeMarkerPairs[v1];
            const uint64_t ordinalSum1 = markerPairInfo1.ordinalSum();
            if(ordinalSum1 - ordinalSum0 > maxOffsetSumDelta) {
                break;
            }

            // Check skip.
            const uint64_t skip = max(
                abs(int32_t(markerPairInfo0.ordinal0) - int32_t(markerPairInfo1.ordinal0)),
                abs(int32_t(markerPairInfo0.ordinal1) - int32_t(markerPairInfo1.ordinal1)));
            if(skip > maxSkip) {
                continue;
            }

            // Check drift.
            const int64_t offset0 =  markerPairInfo0.ordinalOffset();
            const int64_t offset1 =  markerPairInfo1.ordinalOffset();
            const uint64_t drift = labs(offset0 - offset1);
            if(drift > maxDrift) {
                continue;
            }

            // Add an edge corresponding to these marker pairs.
            add_edge(v0, v1, graph);
        }
    }


    // Conpute the longest path. This gives us the alignment.
    vector<uint64_t> longestPath;
    shasta::longestPath(graph, longestPath);
    if(html) {
        html << "The longest path uses " << longestPath.size() <<
            " vertices out of " << activeMarkerPairs.size() << endl;
    }

    // Create the alignment from the longest path.
    alignment.clear();
    for(const uint64_t v:longestPath) {
        const MarkerPairInfo& markerPairInfo = activeMarkerPairs[v];
        alignment.ordinals.push_back({markerPairInfo.ordinal0, markerPairInfo.ordinal1});
    }


    // Store the alignment info.
    alignmentInfo.create(
        alignment,
        uint32_t(orientedReadSortedMarkersSpans[0].size()),
        uint32_t(orientedReadSortedMarkersSpans[1].size()));
    // alignmentInfo.uniquenessMetric = uniquenessMetric;


}
