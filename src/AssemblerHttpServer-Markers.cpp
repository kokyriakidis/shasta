// Shasta.
#include "Assembler.hpp"
#include "KmerChecker.hpp"
#include "KmerCounter.hpp"
#include "MarkerKmers.hpp"
#include "Reads.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/algorithm/string.hpp>



void Assembler::exploreReadMarkers(const vector<string>& request, ostream& html)
{
    SHASTA_ASSERT(assemblerInfo->readRepresentation == 0);

    // Get the request parameters.
    ReadId readId = 0;
    const bool readIdIsPresent = getParameterValue(request, "readId", readId);
    Strand strand = 0;
    const bool strandIsPresent = getParameterValue(request, "strand", strand);


    // Write the form.
    html <<
        "<form>"
        "<table>"

        "<tr>"
        "<th class=left>Numeric read id"
        "<td><input type=text name=readId" <<
        (readIdIsPresent ? (" value=" + to_string(readId)) : "") <<
        " title='Enter a read id between 0 and " << reads->readCount()-1 << "'>"

        "<tr>"
        "<th class=left>Strand"
        "<td>";
    writeStrandSelection(html, "strand", strandIsPresent && strand==0, strandIsPresent && strand==1);

    html <<
        "</table>"
        "<input type=submit value='Display'>"
        "</form>";

    if(not readIdIsPresent) {
        html << "Specify a numeric read id.";
        return;
    }

    // If the strand is missing, stop here.
    if(not strandIsPresent) {
        return;
    }

    // Sanity checks.
    if(readId >= reads->readCount()) {
        html << "<p>Invalid read id.";
        return;
    }
    if(strand!=0 && strand!=1) {
        html << "<p>Invalid strand.";
        return;
    }


    // Access the read information we need.
    const OrientedReadId orientedReadId(readId, strand);
    const auto sequence = reads->getRead(readId);
    const span<const CompressedMarker> orientedReadMarkers = markers[orientedReadId.getValue()];




    // Page title.
    html << "<h2 title='Markers of read " << readId << " on strand " << strand;
    if(strand == 0) {
        html << " (input read without reverse complementing)";
    } else {
        html << " (reverse complement of input read)";
    }
    html << "'>Markers of oriented read " << orientedReadId << "</h2>";

    // Write a table with some summary information for the markers of this oriented read.
    const double readMarkerDensity = double(orientedReadMarkers.size()) / double(sequence.baseCount);
    const double assemblyMarkerDensity = double(markers.totalSize()) / double(2 * assemblerInfo->baseCount);
    const uint64_t expectedMarkerCount = uint64_t(std::round(assemblyMarkerDensity * double(sequence.baseCount)));
    html <<
        "<table>"
        "<tr><th class=left>Length in bases<td class=centered>" << sequence.baseCount <<
        "<tr><th class=left>Number of markers<td class=centered>" << orientedReadMarkers.size() <<
        "<tr><th class=left>Average marker density for this read<td class=centered>" <<
        readMarkerDensity <<
        "<tr><th class=left>Average marker density for this assembly<td class=centered>" <<
        assemblyMarkerDensity <<
        "<tr><th class=left width=400>"
        "Expected number of markers based on average marker density for this assembly<td class=centered>" <<
        expectedMarkerCount <<
        "<tr><th class=left>"
        "Deviation from expected number of markers<td class=centered>" <<
        int64_t(orientedReadMarkers.size()) - int64_t(expectedMarkerCount) <<
        "<tr><th class=left>"
        "Deviation from expected number of markers, relative to "
        "standard deviation of a Poisson distribution<td class=centered>" <<
        double(int64_t(orientedReadMarkers.size()) - int64_t(expectedMarkerCount)) /
        sqrt(double(expectedMarkerCount)) <<
        "</table>";


    // Count k-mers in this oriented read.
    // Reverse complemented k-mers are considered equivalent.
    const uint64_t k = assemblerInfo->k;
    std::map<Kmer, uint64_t> kmerFrequencyMap;
    for(uint64_t ordinal=0; ordinal<orientedReadMarkers.size(); ordinal++) {
        const Kmer kmer = getOrientedReadMarkerKmer(orientedReadId, ordinal);
        const KmerId kmerId = kmer.id(k);
        const Kmer rcKmer = kmer.reverseComplement(k);
        const KmerId rcKmerId = rcKmer.id(k);
        const Kmer canonicalKmer = (kmerId <= rcKmerId) ? kmer : rcKmer;
        const auto it = kmerFrequencyMap.find(canonicalKmer);
        if(it == kmerFrequencyMap.end()) {
            kmerFrequencyMap.insert({canonicalKmer, 1});
        } else {
            ++(it->second);
        }
    }



    // Begin the main table containing one row for each marker.
    html <<
        "<p><table>"
        "<tr>"
        "<th>Marker<br>ordinal"
        "<th>Begin<br>position"
        "<th>End<br>position"
        "<th>Kmer"
        "<th>Frequency<br>in this<br>oriented read"
        ;
    if(markerKmers and markerKmers->isOpen()) {
        html << "<th>Global<br>frequency";
    }

    if(kmerCounter and kmerCounter->isAvailable()) {
        html << "<th>Global<br>frequency";
    }



    // Write one row for each marker.
    for(uint64_t ordinal=0; ordinal<orientedReadMarkers.size(); ordinal++) {
        const uint64_t position = orientedReadMarkers[ordinal].position;
        const Kmer kmer = getOrientedReadMarkerKmer(orientedReadId, ordinal);
        const KmerId kmerId = kmer.id(k);
        const Kmer rcKmer = kmer.reverseComplement(k);
        const KmerId rcKmerId = rcKmer.id(k);
        const Kmer canonicalKmer = (kmerId <= rcKmerId) ? kmer : rcKmer;

        html <<
            "<tr>"
            "<td class=centered>" << ordinal <<
            "<td class=centered>" << position <<
            "<td class=centered>" << position + k <<
            "<td class=centered style='font-family:Courier New;'>";
        kmer.write(html, k);

        // Frequency of this Kmer in this read.
        html << "<td class=centered>" << kmerFrequencyMap[canonicalKmer];

        // Global frequency of this Kmer.
        if(kmerCounter and kmerCounter->isAvailable()) {
            html << "<td class=centered>" << kmerCounter->getFrequency(kmer);
        } else if(markerKmers and markerKmers->isOpen()) {
            html << "<td class=centered>" << markerKmers->getFrequency(kmer);
        }

    }

    html << "</table>";

}



void Assembler::exploreMarkerKmers(const vector<string>& request, ostream& html)
{
    SHASTA_ASSERT(assemblerInfo->readRepresentation == 0);
    SHASTA_ASSERT(markerKmers and markerKmers->isOpen());

    const uint64_t k = assemblerInfo->k;

    html << "<h2>Marker k-mer</h2>";

    // Get the request parameters.
    string kmerString;
    getParameterValue(request, "kmer", kmerString);
    boost::trim(kmerString);

    // Write the form.
    html <<
        "<p><form><table>"
        "<tr><th class=left>K-mer"
        "<td><input type=text name=kmer style='font-family:monospace' "
        "size=" << k + 10 << " "
        "value='" << kmerString << "'" <<
        " title='Enter a " << k << "-base k-mer.'>"
        "</table><input type=submit value='Get k-mer information'></form>";

    // If the k-mer string is empty, do nothing.
    if(kmerString.empty()) {
        return;
    }

    // Check the length.
    if(kmerString.size() != k) {
        html << "This k-mer is " << kmerString.size() << " bases long. "
            "This assembly uses " << k << "-base markers. "
            "Specify a " << k << "-mer." << endl;
        return;
    }

    // Construct the k-mer.
    Kmer kmer;
    for(uint64_t i=0; i<kmerString.size(); i++) {
        const char c = kmerString[i];
        const Base b = Base::fromCharacterNoException(c);
        if(not b.isValid()) {
            throw runtime_error("Invalid base character " + string(1, c) + " at k-mer position " + to_string(i));
        }
        kmer.set(i, b);
    }

    // Check if it is a marker.
    SHASTA_ASSERT(kmerChecker);
    const KmerId kmerId = KmerId(kmer.id(k));
    if(not kmerChecker->isMarker(kmerId)) {
        throw runtime_error("This assembly does not use this as a marker.");
    }



    // Summary table.
    const uint64_t coverage = markerKmers->getFrequency(kmer);
    const Kmer kmerRc = kmer.reverseComplement(k);
    html <<
        "<table><tr><th class=left>K-mer<td class=centered style='font-family:monospace'>";
    kmer.write(html, k);
    html <<
        "<tr><th class=left>Reverse complement K-mer<td class=centered style='font-family:monospace'>";
    kmerRc.write(html, k);
    html << "<tr><th class=left>Coverage<td class=centered>" << coverage <<
        "</table>";



    // Details table.
    vector<MarkerKmers::MarkerInfo> markerInfos;
    markerKmers->get(kmer, markerInfos);

    html <<
        "<p>"
        "<table>"
        "<tr><th>Oriented<br>read<th>Ordinal";
    for(const MarkerKmers::MarkerInfo& markerInfo: markerInfos) {
        html <<
            "<tr>"
            "<td class=centered>" << markerInfo.orientedReadId <<
            "<td class=centered>" << markerInfo.ordinal;
    }
    html << "</table>";


}
