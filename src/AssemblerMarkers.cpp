// Shasta.
#include "Assembler.hpp"
#include "extractKmer.hpp"
#include "findMarkerId.hpp"
#include "MarkerFinder.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
using namespace shasta;

// Standard library.
#include "fstream.hpp"


void Assembler::findMarkers(size_t threadCount)
{
    reads->checkReadsAreOpen();
    checkKmersAreOpen();

    markers.createNew(largeDataName("Markers"), largeDataPageSize);
    MarkerFinder markerFinder(
        assemblerInfo->k,
        kmerTable,
        getReads(),
        markers,
        threadCount);

}



void Assembler::accessMarkers()
{
    markers.accessExistingReadOnly(largeDataName("Markers"));
}

void Assembler::checkMarkersAreOpen() const
{
    if(!markers.isOpen()) {
        throw runtime_error("Markers are not accessible.");
    }
}


void Assembler::writeMarkers(ReadId readId, Strand strand, const string& fileName)
{
    // Check that we have what we need.
    checkKmersAreOpen();
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    reads->checkReadId(readId);

    // Get the markers.
    const OrientedReadId orientedReadId(readId, strand);
    const auto orientedReadMarkers = markers[orientedReadId.getValue()];

    // Write them out.
    ofstream csv(fileName);
    csv << "MarkerId,Ordinal,KmerId,Kmer,Position\n";
    for(uint32_t ordinal=0; ordinal<orientedReadMarkers.size(); ordinal++) {
        const CompressedMarker& marker = orientedReadMarkers[ordinal];
        const MarkerId markerId = getMarkerId(orientedReadId, ordinal);
        csv << markerId << ",";
        csv << ordinal << ",";
        csv << marker.kmerId << ",";
        csv << Kmer(marker.kmerId, assemblerInfo->k) << ",";
        csv << marker.position << "\n";
    }
}



vector<KmerId> Assembler::getMarkers(ReadId readId, Strand strand)
{
    const OrientedReadId orientedReadId(readId, strand);
    const auto orientedReadMarkers = markers[orientedReadId.getValue()];

    vector<KmerId> v;
    for(const CompressedMarker& marker: orientedReadMarkers) {
        v.push_back(marker.kmerId);
    }
    return v;
}


// Get markers sorted by KmerId for a given OrientedReadId.
void Assembler::getMarkersSortedByKmerId(
    OrientedReadId orientedReadId,
    vector<MarkerWithOrdinal>& markersSortedByKmerId) const
{
    const auto compressedMarkers = markers[orientedReadId.getValue()];
    markersSortedByKmerId.clear();
    markersSortedByKmerId.resize(compressedMarkers.size());

    for(uint32_t ordinal=0; ordinal<compressedMarkers.size(); ordinal++) {
        const CompressedMarker& compressedMarker = compressedMarkers[ordinal];
        markersSortedByKmerId[ordinal] = MarkerWithOrdinal(compressedMarker, ordinal);
    }

    // Sort by kmerId.
    sort(markersSortedByKmerId.begin(), markersSortedByKmerId.end());
}



// Given a marker by its OrientedReadId and ordinal,
// return the corresponding global marker id.
MarkerId Assembler::getMarkerId(
    OrientedReadId orientedReadId, uint32_t ordinal) const
{
    return
        (markers.begin(orientedReadId.getValue()) - markers.begin())
        + ordinal;
}

MarkerId Assembler::getReverseComplementMarkerId(
    OrientedReadId orientedReadId, uint32_t ordinal) const
{
    OrientedReadId orientedReadIdRc = orientedReadId;
    orientedReadIdRc.flipStrand();

    const uint32_t markerCount = uint32_t(markers.size(orientedReadId.getValue()));

    return getMarkerId(orientedReadIdRc, markerCount - 1 - ordinal);

}


// Inverse of the above: given a global marker id,
// return its OrientedReadId and ordinal.
// This requires a binary search in the markers toc.
// This could be avoided, at the cost of storing
// an additional 4 bytes per marker.
pair<OrientedReadId, uint32_t>
    Assembler::findMarkerId(MarkerId markerId) const
{
    return shasta::findMarkerId(markerId, markers);
}



// Given a MarkerId, compute the MarkerId of the
// reverse complemented marker.
MarkerId Assembler::findReverseComplement(MarkerId markerId) const
{
	// Find the oriented read id and marker ordinal.
	OrientedReadId orientedReadId;
	uint32_t ordinal;
	tie(orientedReadId, ordinal) = findMarkerId(markerId);

	// Reverse complement.
	ordinal = uint32_t(markers.size(orientedReadId.getValue()) - 1 - ordinal);
	orientedReadId.flipStrand();

	// Return the corresponding Markerid.
	return getMarkerId(orientedReadId, ordinal);
}



// Write the frequency of markers in oriented reads.
void Assembler::writeMarkerFrequency()
{
    const uint64_t k = assemblerInfo->k;
    const uint64_t kmerCount = 1ULL << (2ULL*k);
    SHASTA_ASSERT(markers.isOpen());
    vector<uint64_t> frequency(kmerCount, 0);

    const CompressedMarker* compressedMarker = markers.begin();
    const CompressedMarker* end = markers.end();
    for(; compressedMarker!=end; ++compressedMarker) {
        ++frequency[compressedMarker->kmerId];
    }

    ofstream csv("MarkerFrequency.csv");
    for(uint64_t kmerId=0; kmerId<kmerCount; kmerId++) {
        const uint64_t n = frequency[kmerId];
        if(n== 0) {
            continue;
        }
        const Kmer kmer(kmerId, k);
        kmer.write(csv, k);
        csv << "," << n << "\n";
    }
}



void Assembler::computeMarkerKmers(uint64_t threadCount)
{
    performanceLog << timestamp << "computeMarkerKmers begins." << endl;

    // Check that we have what we need.
    checkMarkersAreOpen();
    const uint64_t readCount = reads->readCount();

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Do it.
    // The layout is identical to that used by the markers.
    markerKmers.createNew(largeDataName("MarkerKmers"), largeDataPageSize);
    for(uint64_t readId=0; readId<readCount; readId++) {
        const OrientedReadId orientedReadId0(uint32_t(readId), 0);
        const OrientedReadId orientedReadId1(uint32_t(readId), 1);
        const uint64_t readMarkerCount = markers.size(orientedReadId0.getValue());
        SHASTA_ASSERT(markers.size(orientedReadId1.getValue()) == readMarkerCount);
        for(uint64_t strand=0; strand<2; strand++) {
            markerKmers.appendVector(readMarkerCount);
        }
    }
    markerKmers.unreserve();
    const uint64_t batchSize = 100;
    setupLoadBalancing(readCount, batchSize);
    runThreads(&Assembler::computeMarkerKmersThreadFunction, threadCount);

    performanceLog << timestamp << "computeMarkerKmers ends." << endl;
}



void Assembler::computeMarkerKmersThreadFunction(size_t threadId)
{
    const uint64_t k = assemblerInfo->k;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over reads in this batch.
        for(uint64_t readId=begin; readId!=end; ++readId) {
            auto read = reads->getRead(uint32_t(readId));

            // Access the two oriented reads for this read.
            const OrientedReadId orientedReadId0(uint32_t(readId), 0);
            const OrientedReadId orientedReadId1(uint32_t(readId), 1);


            // Access their markers.
            const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];
            const auto orientedReadMarkers1 = markers[orientedReadId1.getValue()];
            const uint64_t readMarkerCount = orientedReadMarkers0.size();
            SHASTA_ASSERT(orientedReadMarkers1.size() == readMarkerCount);

            // Access the marker kmers we are filling in.
            const auto orientedReadKmers0 = markerKmers[orientedReadId0.getValue()];
            const auto orientedReadKmers1 = markerKmers[orientedReadId1.getValue()];



            // Loop over markers of this read.
            for(uint64_t ordinal0=0; ordinal0<readMarkerCount; ordinal0++) {

                // Strand 0.
                const KmerId kmerId0FromMarker = orientedReadMarkers0[ordinal0].kmerId;
                Kmer kmer0;
                extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
                const KmerId kmerId0FromRead = KmerId(kmer0.id(k));
                SHASTA_ASSERT(kmerId0FromMarker == kmerId0FromRead);
                orientedReadKmers0[ordinal0] = kmerId0FromRead;

                // Strand 1.
                const Kmer kmer1 = kmer0.reverseComplement(k);
                const uint64_t ordinal1 = readMarkerCount - 1 - ordinal0;
                const KmerId kmerId1FromMarker = orientedReadMarkers1[ordinal1].kmerId;
                const KmerId kmerId1FromRead = KmerId(kmer1.id(k));
                SHASTA_ASSERT(kmerId1FromMarker == kmerId1FromRead);
                orientedReadKmers1[ordinal1] = kmerId1FromRead;
            }
        }
    }

}
