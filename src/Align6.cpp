// Shasta.
#include "Align6.hpp"
#include "Align6Marker.hpp"
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



Align6::Align6(
    const array<span<Align6Marker>, 2>& orientedReadMarkers,
    uint64_t k,
    uint64_t maxSkip,
    uint64_t maxDrift,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    ostream& html) :

    orientedReadMarkers(orientedReadMarkers),
    k(k),
    maxSkip(maxSkip),
    maxDrift(maxDrift),
    alignment(alignment),
    alignmentInfo(alignmentInfo),
    html(html),
    maxOffsetSumDelta(2 * maxSkip)
{

    // Find ordinal offsets of the low frequency marker pairs.
    // If there are two few of them, store an empty alignment and return.
    computeLowFrequencyMarkerPairOffsets();
    if(lowFrequencyMarkerPairOffsets.size() < minLowFrequencyCount) {
        if(html) {
            html << "<br>Two few low frequency marker pairs found.";
        }
        storeEmptyAlignment();
        return;
    }

    // Compute a histogram of ordinal offsets of the low frequency marker pairs.
    computeOffsetHistogram();

    // Create a convolution kernel that will be used to get a smoothed version
    // of the histogram.
    createConvolutionKernel();

    // Compute the offset distribution, a smoothed version of the offset histogram.
    computeOffsetDistribution();

    // Use the offset distribution to compute a band.
    computeBand();

    // Gather all marker pairs in the band, regardless of frequency.
    gatherMarkerPairsInBand();

    // Compute connected components of the marker pairs in the band.
    // Two marker pairs belong to the same component if
    // their ordinals are compatible with maxSkip and maxDrift.
    computeComponents();

    // Write the marker pairs in the band and the component they belong to.
    writeMarkerPairsInBand();

    // Find the connected component with the most low frequency markers.
    const uint64_t bestLowFrequencyMarkerCount = findBestComponent();
    if(bestLowFrequencyMarkerCount < minLowFrequencyCount) {
        if(html) {
            html << "<br>Too few low frequency marker pairs found in the best component.";
        }
        storeEmptyAlignment();
        return;
    }

    // Gather the marker pairs in the best component.
    // These are the ones that will be used to compute the alignment.
    gatherActiveMarkerPairs();

    // Use the active marker pairs to compute the alignment.
    computeAlignment();
}



// Find low frequency marker pairs and store their offsets.
// The low frequency marker pairs consists of pairs (ordinal0, ordinal1)
// such that KmerId(orientedReadId0, ordinal0) == KmerId(orientedReadId1, ordinal1),
// and that satisfy the requirements on local and global frequency.
void Align6::computeLowFrequencyMarkerPairOffsets()
{
    lowFrequencyMarkerPairOffsets.clear();

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
    }

    const auto begin0 = orientedReadMarkers[0].begin();
    const auto begin1 = orientedReadMarkers[1].begin();
    const auto end0 = orientedReadMarkers[0].end();
    const auto end1 = orientedReadMarkers[1].end();

    auto it0 = begin0;
    auto it1 = begin1;
    while(it0!=end0 && it1!=end1) {
        if(it0->kmerId < it1->kmerId) {
            ++it0;
        } else if(it1->kmerId < it0->kmerId) {
            ++it1;
        } else {

            // We found a common KmerId.
            const KmerId kmerId = it0->kmerId;
            const uint64_t globalFrequency = it0->globalFrequency;

            // This KmerId could appear more than once in each of the sequences,
            // so we need to find the streak of this KmerId.
            auto it0Begin = it0;
            auto it1Begin = it1;
            auto it0End = it0Begin;
            auto it1End = it1Begin;
            while(it0End!=end0 && it0End->kmerId == kmerId) {
                ++it0End;
            }
            while(it1End!=end1 && it1End->kmerId == kmerId) {
                ++it1End;
            }

            // If the local global frequencies are in the desired range,
            // loop over pairs in the streaks.
            const uint64_t localFrequency0 = it0End - it0Begin;
            const uint64_t localFrequency1 = it1End - it1Begin;
            if(
                (localFrequency0 <= maxLocalFrequency) and
                (localFrequency1 <= maxLocalFrequency) and
                (globalFrequency >= minGlobalFrequency) and
                (globalFrequency <= maxGlobalFrequency)) {

                MarkerPair markerPair;
                markerPair.kmerId = kmerId;
                markerPair. globalFrequency = globalFrequency;
                markerPair.localFrequency0 = localFrequency0;
                markerPair.localFrequency1 = localFrequency1;

                for(auto jt0=it0Begin; jt0!=it0End; ++jt0) {
                    markerPair. ordinal0 = jt0->ordinal;
                    for(auto jt1=it1Begin; jt1!=it1End; ++jt1) {
                        markerPair.ordinal1 = jt1->ordinal;
                        lowFrequencyMarkerPairOffsets.push_back(markerPair.ordinalOffset());

                        if(html) {
                            html <<
                                "<tr>"
                                "<td class=centered>" << markerPair.ordinal0 <<
                                "<td class=centered>" << markerPair.ordinal1 <<
                                "<td class=centered>" << markerPair.ordinalSum() <<
                                "<td class=centered>" << markerPair.ordinalOffset() <<
                                "<td class=centered style='font-family:courier'>";

                            Kmer(markerPair.kmerId, k).write(html, k);

                            html <<
                                "<td class=centered>" << markerPair.kmerId <<
                                "<td class=centered>" << markerPair.localFrequency0 <<
                                "<td class=centered>" << markerPair.localFrequency1 <<
                                "<td class=centered>" << markerPair.globalFrequency;
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

}



// This is called when we give up.
// Stores an empty alignment and the corresponding AlignmentInfo.
void Align6::storeEmptyAlignment()
{
    alignment.clear();
    alignmentInfo.create(
        alignment,
        uint32_t(orientedReadMarkers[0].size()),
        uint32_t(orientedReadMarkers[1].size()));
    alignmentInfo.uniquenessMetric = 0.;
}



// Compute a histogram of ordinal offsets for the low frequency marker pairs.
void Align6::computeOffsetHistogram()
{
    std::map<int64_t, uint64_t> offsetHistogramMap;
    for(const int64_t offset: lowFrequencyMarkerPairOffsets) {
        auto it = offsetHistogramMap.find(offset);
        if(it == offsetHistogramMap.end()) {
            offsetHistogramMap.insert({offset, 1});
        } else {
            ++it->second;
        }
    }

    // Copy it to the offsetHistogram vector.
    offsetHistogram.clear();
    copy(offsetHistogramMap.begin(), offsetHistogramMap.end(), back_inserter(offsetHistogram));

    // Write out the histogram.
    if(html) {
        html << "<h3>Histogram of ordinal offsets for the low frequency marker pairs</h3>"
            "<table>"
            "<tr><th>Ordinal<br>offset<th>Frequency";
        for(const auto& p: offsetHistogram) {
            const int64_t offset = p.first;
            const uint64_t frequency = p.second;

            html << "<tr><td class=centered>" << offset;
            html << "<td class=centered>" << frequency;
        }
        html << "</table>";
    }
}



// Create a convolution kernel that will be used to get a smoothed version
// of the offset histogram.
void Align6::createConvolutionKernel()
{
    // Define the half width of the kernel.
    const uint64_t minLength = min(
        orientedReadMarkers[0].size(),
        orientedReadMarkers[1].size());
    const int64_t kernelHalfWidth = int64_t(std::round(driftRateTolerance * double(minLength)));

    // Fill the kernel with a linear function that is 1 at 0 and goes to 0 at the edge of the kernel.
    kernel.resize(2 * kernelHalfWidth + 1);
    for(uint64_t i=0; i<kernel.size(); i++) {
        const double x = double(int64_t(i) - kernelHalfWidth) / double(kernelHalfWidth);
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
                "<tr><td class=centered>" << int64_t(i) - kernelHalfWidth <<
                "<td class=centered>" << kernel[i];
        }
        html << "</table>";
    }
}



void Align6::computeOffsetDistribution()
{
    offsetDistribution.resize(maxDistributionOffset() - minDistributionOffset() + 1, 0.);
    for(const auto& p: offsetHistogram) {
        const int64_t offset = p.first;
        const uint64_t frequency = p.second;
        const double weight = double(frequency);
        for(uint64_t i=0; i<kernel.size(); i++) {
            const int64_t j = offset + int64_t(i) - kernelHalfWidth() - minDistributionOffset();
            SHASTA_ASSERT(j >= 0);
            SHASTA_ASSERT(j < int64_t(offsetDistribution.size()));
            offsetDistribution[j] += weight * kernel[i];
        }
    }

    // Write the smoothed offset distribution.
    if(html) {
        html <<
            "<h3>Smoothed offset distribution</h3>"
            "<table>"
            "<tr><th>Offset<th>Value";
        for(uint64_t i=0; i<offsetDistribution.size(); i++) {
            const int64_t offset = int64_t(i) + minDistributionOffset();
            html << "<tr><td class=centered>" << offset <<
                "<td class=centered>" << offsetDistribution[i];
        }
        html << "</table>";
    }
}



// Use the offset distribution to compute a band.
void Align6::computeBand()
{
    bandCenter =
        int64_t(std::max_element(offsetDistribution.begin(), offsetDistribution.end()) -
            offsetDistribution.begin()) + minDistributionOffset();
    bandLow = bandCenter - kernelHalfWidth();
    bandHigh = bandCenter + kernelHalfWidth();

    if(html) {
        html << "<h3>Band</h3>"
            "Band includes offset " << bandLow << " through " << bandHigh <<
            " and is centered at " << bandCenter;
    }
}



// Gather all marker pairs in the band, regardless of frequency.
void Align6::gatherMarkerPairsInBand()
{
    inBandMarkerPairs.clear();

    const auto begin0 = orientedReadMarkers[0].begin();
    const auto begin1 = orientedReadMarkers[1].begin();
    const auto end0 = orientedReadMarkers[0].end();
    const auto end1 = orientedReadMarkers[1].end();

    auto it0 = begin0;
    auto it1 = begin1;
    while(it0!=end0 && it1!=end1) {
        if(it0->kmerId < it1->kmerId) {
            ++it0;
        } else if(it1->kmerId < it0->kmerId) {
            ++it1;
        } else {

            // We found a common KmerId.
            const KmerId kmerId = it0->kmerId;
            const uint64_t globalFrequency = it0->globalFrequency;

            // This KmerId could appear more than once in each of the sequences,
            // so we need to find the streak of this KmerId.
            auto it0Begin = it0;
            auto it1Begin = it1;
            auto it0End = it0Begin;
            auto it1End = it1Begin;
            while(it0End!=end0 && it0End->kmerId == kmerId) {
                ++it0End;
            }
            while(it1End!=end1 && it1End->kmerId == kmerId) {
                ++it1End;
            }

            // Loop over pairs in the streaks.
            const uint64_t localFrequency0 = it0End - it0Begin;
            const uint64_t localFrequency1 = it1End - it1Begin;

            MarkerPair markerPair;
            markerPair.kmerId = kmerId;
            markerPair. globalFrequency = globalFrequency;
            markerPair.localFrequency0 = localFrequency0;
            markerPair.localFrequency1 = localFrequency1;

            for(auto jt0=it0Begin; jt0!=it0End; ++jt0) {
                markerPair. ordinal0 = jt0->ordinal;
                for(auto jt1=it1Begin; jt1!=it1End; ++jt1) {
                    markerPair.ordinal1 = jt1->ordinal;
                    const int64_t offset = markerPair.ordinalOffset();
                    if((offset >= bandLow) and (offset <= bandHigh)) {
                        inBandMarkerPairs.push_back(markerPair);
                    }
                }
            }

            // Continue the joint loop over KmerId's.
            it0 = it0End;
            it1 = it1End;
        }
    }

    // Sort them by ordinalSum.
    sort(inBandMarkerPairs.begin(), inBandMarkerPairs.end());
}



// Find out if two MarkerPairInfos are compatible with maxSkip and maxDrift.
bool Align6::canBeConnected(
    const MarkerPair& markerPairInfo0,
    const MarkerPair& markerPairInfo1) const
{

    // Check skip.
    const uint64_t skip = max(
        abs(int32_t(markerPairInfo0.ordinal0) - int32_t(markerPairInfo1.ordinal0)),
        abs(int32_t(markerPairInfo0.ordinal1) - int32_t(markerPairInfo1.ordinal1)));
    if(skip > maxSkip) {
        return false;
    }

    // Check drift.
    const int64_t offset0 =  markerPairInfo0.ordinalOffset();
    const int64_t offset1 =  markerPairInfo1.ordinalOffset();
    const uint64_t drift = labs(offset0 - offset1);
    if(drift > maxDrift) {
        return false;
    }

    return true;
}



// Compute connected components of the marker pairs in the band.
// Two marker pairs belong to the same component if
// their ordinals are compatible with maxSkip and maxDrift.
void Align6::computeComponents()
{

    // Initialize the disjoint sets data structure.
    vector<uint64_t> rank(inBandMarkerPairs.size());
    vector<uint64_t> parent(inBandMarkerPairs.size());
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t v=0; v<inBandMarkerPairs.size(); v++) {
        disjointSets.make_set(v);
    }

    // Loop over pairs of MarkerPairInfos that are possibly compatible with maxSkip and maxDrift.
    for(uint64_t v0=0; v0<inBandMarkerPairs.size(); v0++) {
        const MarkerPair& markerPairInfo0 = inBandMarkerPairs[v0];
        const uint64_t ordinalSum0 = markerPairInfo0.ordinalSum();
        for(uint64_t v1=v0+1; v1<inBandMarkerPairs.size(); v1++) {
            const MarkerPair& markerPairInfo1 = inBandMarkerPairs[v1];
            const uint64_t ordinalSum1 = markerPairInfo1.ordinalSum();

            if(ordinalSum1 - ordinalSum0 > maxOffsetSumDelta) {
                break;
            }

            // If they are compatible with maxSkip and maxDrift,
            // update the disjoint sets to take this edge into account.
            if(canBeConnected(markerPairInfo0, markerPairInfo1)) {
                disjointSets.union_set(v0, v1);
            }
        }
    }

    // Compute connected components.
    component.resize(inBandMarkerPairs.size());
    for(uint64_t v=0; v<inBandMarkerPairs.size(); v++) {
        component[v] = disjointSets.find_set(v);
    }
}



// Write the marker pairs in the band and the component they belong to.
void Align6::writeMarkerPairsInBand()
{
    // Write out the marker pairs in the band.
    if(html) {
        html <<
            "<h3>Marker pairs in the band</h3>"
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

        for(uint64_t i=0; i<inBandMarkerPairs.size(); i++) {
            const MarkerPair& markerPair = inBandMarkerPairs[i];
            html <<
                "<tr>"
                "<td class=centered>" << i <<
                "<td class=centered>" << markerPair.ordinal0 <<
                "<td class=centered>" << markerPair.ordinal1 <<
                "<td class=centered>" << markerPair.ordinalSum() <<
                "<td class=centered>" << markerPair.ordinalOffset() <<
                "<td class=centered style='font-family:courier'>";

            Kmer(markerPair.kmerId, k).write(html, k);

            html <<
                "<td class=centered>" << markerPair.kmerId <<
                "<td class=centered>" << markerPair.localFrequency0 <<
                "<td class=centered>" << markerPair.localFrequency1 <<
                "<td class=centered>" << markerPair.globalFrequency <<
                "<td class=centered>" << component[i];
        }

        html << "</table>";
    }
}



// The best component is the one with the most low frequency markers.
// This returns the number of low frequency markers in the best component.
uint64_t Align6::findBestComponent()
{

    // Count low frequency marker pairs in each connected component.
    vector<uint64_t> count(inBandMarkerPairs.size(), 0);
    for(uint64_t i=0; i<inBandMarkerPairs.size(); i++) {
        const MarkerPair& markerPair = inBandMarkerPairs[i];
        if(markerPair.globalFrequency > maxGlobalFrequency) {
            continue;
        }
        if(markerPair.globalFrequency < minGlobalFrequency) {
            continue;
        }
        if(markerPair.localFrequency0 > maxLocalFrequency) {
            continue;
        }
        if(markerPair.localFrequency1 > maxLocalFrequency) {
            continue;
        }
        const uint64_t c = component[i];
        ++count[c];
    }

    if(html) {
        html << "<h3>Number of low frequency marker pairs in each component</h3>"
            "<table><tr><th>Component<th>Low<br>frequency<br>markers";
        for(uint64_t c=0; c<inBandMarkerPairs.size(); c++) {
            if(count[c] > 0) {
                html << "<tr><td class=centered>" << c << "<td class=centered>" << count[c];
            }
        }
        html << "</table>";
    }


    // Find the connected component with the most low frequency markers.
    bestComponent = max_element(count.begin(), count.end()) - count.begin();
    if(html) {
        html << "<br>Will use marker pairs from component " << bestComponent <<
            " with " << count[bestComponent] << " low frequency markers.";
    }

    return count[bestComponent];
}



// The active marker pairs are the ones  in the best component.
void Align6::gatherActiveMarkerPairs()
{
    activeMarkerPairs.clear();

    for(uint64_t v=0; v<inBandMarkerPairs.size(); v++) {
        if(component[v] == bestComponent) {
            activeMarkerPairs.push_back(inBandMarkerPairs[v]);
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
            const MarkerPair& markerPair = activeMarkerPairs[i];
            html <<
                "<tr>"
                "<td class=centered>" << i <<
                "<td class=centered>" << markerPair.ordinal0 <<
                "<td class=centered>" << markerPair.ordinal1 <<
                "<td class=centered>" << markerPair.ordinalSum() <<
                "<td class=centered>" << markerPair.ordinalOffset() <<
                "<td class=centered style='font-family:courier'>";

            Kmer(markerPair.kmerId, k).write(html, k);

            html <<
                "<td class=centered>" << markerPair.kmerId <<
                "<td class=centered>" << markerPair.localFrequency0 <<
                "<td class=centered>" << markerPair.localFrequency1 <<
                "<td class=centered>" << markerPair.globalFrequency;
        }

        html << "</table>";
    }
}



// Use the active marker pairs to compute the alignment.
void Align6::computeAlignment()
{
    // Use the active marker pairs to create a directed graph
    // whose edges correspond to MarkerPairInfos compatible with maxSkip, maxDrift.
    using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS>;
    Graph graph(activeMarkerPairs.size());

    for(uint64_t v0=0; v0<activeMarkerPairs.size(); v0++) {
        const MarkerPair& markerPairInfo0 = activeMarkerPairs[v0];
        const uint64_t ordinalSum0 = markerPairInfo0.ordinalSum();

        for(uint64_t v1=v0+1; v1<activeMarkerPairs.size(); v1++) {
            const MarkerPair& markerPairInfo1 = activeMarkerPairs[v1];

            const uint64_t ordinalSum1 = markerPairInfo1.ordinalSum();
            if(ordinalSum1 - ordinalSum0 > maxOffsetSumDelta) {
                break;
            }

            // If these two MarkerPair are compatible with maxSkip, maxDrift,
            // Add an edge corresponding to these marker pairs.
            if(canBeConnected(markerPairInfo0, markerPairInfo1)) {
                add_edge(v0, v1, graph);
            }
        }
    }



    // Compute the longest path. This gives us the alignment.
    vector<uint64_t> longestPath;
    shasta::longestPath(graph, longestPath);
    if(html) {
        html << "The longest path uses " << longestPath.size() <<
            " vertices out of " << activeMarkerPairs.size() << endl;
    }

    // Create the alignment from the longest path.
    alignment.clear();
    for(const uint64_t v:longestPath) {
        const MarkerPair& markerPair = activeMarkerPairs[v];
        alignment.ordinals.push_back({markerPair.ordinal0, markerPair.ordinal1});
    }


    // Store the alignment info.
    alignmentInfo.create(
        alignment,
        uint32_t(orientedReadMarkers[0].size()),
        uint32_t(orientedReadMarkers[1].size()));
    // alignmentInfo.uniquenessMetric = uniquenessMetric;
}



void Align6Marker::setGlobalFrequency(uint64_t globalFrequencyLong)
{
    if(globalFrequencyLong > std::numeric_limits<uint32_t>::max()) {
        globalFrequency = std::numeric_limits<uint32_t>::max();
    } else {
        globalFrequency = uint32_t(globalFrequencyLong);
    }
}
