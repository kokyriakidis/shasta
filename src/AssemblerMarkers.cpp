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



void Assembler::computeMarkerKmerIds(uint64_t threadCount)
{
    performanceLog << timestamp << "computeMarkerKmerIds begins." << endl;

    // Check that we have what we need.
    checkMarkersAreOpen();
    const uint64_t readCount = reads->readCount();

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Do it.
    // The layout is identical to that used by the markers.
    markerKmerIds.createNew(largeDataName("MarkerKmerIds"), largeDataPageSize);
    for(uint64_t readId=0; readId<readCount; readId++) {
        const OrientedReadId orientedReadId0(uint32_t(readId), 0);
        const OrientedReadId orientedReadId1(uint32_t(readId), 1);
        const uint64_t readMarkerCount = markers.size(orientedReadId0.getValue());
        SHASTA_ASSERT(markers.size(orientedReadId1.getValue()) == readMarkerCount);
        for(uint64_t strand=0; strand<2; strand++) {
            markerKmerIds.appendVector(readMarkerCount);
        }
    }
    markerKmerIds.unreserve();
    const uint64_t batchSize = 100;
    setupLoadBalancing(readCount, batchSize);
    runThreads(&Assembler::computeMarkerKmerIdsThreadFunction, threadCount);



#if 1
    // Test the low level functions to extract Kmers/KmerIds.
    const uint64_t k = assemblerInfo->k;
    vector<Kmer> kmerVector;
    vector<KmerId> kmerIdVector;
    performanceLog << timestamp << "Testing." << endl;
    for(uint64_t readId=0; readId<readCount; readId++) {
        for(uint64_t strand=0; strand<2; strand++) {

            const OrientedReadId orientedReadId = OrientedReadId(ReadId(readId), Strand(strand));
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const auto orientedReadMarkerKmerIds = markerKmerIds[orientedReadId.getValue()];
            const uint64_t orientedReadMarkerCount = orientedReadMarkers.size();
            SHASTA_ASSERT(orientedReadMarkerKmerIds.size() == orientedReadMarkerCount);

            kmerVector.resize(orientedReadMarkerCount);
            kmerIdVector.resize(orientedReadMarkerCount);
            const span<Kmer> kmerSpan(kmerVector);
            const span<KmerId> kmerIdSpan(kmerIdVector);

            getOrientedReadMarkerKmers(orientedReadId, kmerSpan);
            getOrientedReadMarkerKmerIds(orientedReadId, kmerIdSpan);

            for(uint64_t ordinal=0; ordinal<orientedReadMarkerCount; ordinal++) {
                SHASTA_ASSERT(kmerVector[ordinal].id(k) == orientedReadMarkers[ordinal].kmerId);
                SHASTA_ASSERT(kmerIdVector[ordinal] == orientedReadMarkers[ordinal].kmerId);

                SHASTA_ASSERT(kmerVector[ordinal] == getOrientedReadMarkerKmer(orientedReadId, ordinal));
                SHASTA_ASSERT(kmerIdVector[ordinal] == getOrientedReadMarkerKmerId(orientedReadId, ordinal));
            }
        }
    }
#endif



    performanceLog << timestamp << "computeMarkerKmerIds ends." << endl;
}



void Assembler::computeMarkerKmerIdsThreadFunction(size_t threadId)
{

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over reads in this batch.
        for(uint64_t readId=begin; readId!=end; ++readId) {

            const OrientedReadId orientedReadId0(uint32_t(readId), 0);
            const OrientedReadId orientedReadId1(uint32_t(readId), 1);

            getReadMarkerKmerIds(
                ReadId(readId),
                markerKmerIds[orientedReadId0.getValue()],
                markerKmerIds[orientedReadId1.getValue()]);
        }
    }

}



// Get all marker Kmers for an oriented read.
void Assembler::getOrientedReadMarkerKmers(
    OrientedReadId orientedReadId,
    const span<Kmer>& kmers) const
{
    const ReadId readId = orientedReadId.getReadId();
    const Strand strand = orientedReadId.getStrand();

    if(strand == 0) {
        getOrientedReadMarkerKmersStrand0(readId, kmers);
    } else {
        getOrientedReadMarkerKmersStrand1(readId, kmers);
    }
}



void Assembler::getOrientedReadMarkerKmersStrand0(ReadId readId, const span<Kmer>& kmers0) const
{
    const uint64_t k = assemblerInfo->k;

    const auto read = reads->getRead(uint32_t(readId));
    const OrientedReadId orientedReadId0(readId, 0);
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];
    const uint64_t readMarkerCount = orientedReadMarkers0.size();
    SHASTA_ASSERT(kmers0.size() == readMarkerCount);

    // Loop over all markers.
    for(uint64_t ordinal0=0; ordinal0<readMarkerCount; ordinal0++) {
        Kmer kmer0;
        extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
        kmers0[ordinal0] = kmer0;
    }



#if 1
    // Check against the KmerIds stored in the markers.
    // These will soon go away.
    for(uint64_t ordinal0=0; ordinal0<readMarkerCount; ordinal0++) {
        SHASTA_ASSERT(kmers0[ordinal0].id(k) == orientedReadMarkers0[ordinal0].kmerId);
    }
#endif
}



void Assembler::getOrientedReadMarkerKmersStrand1(ReadId readId, const span<Kmer>& kmers1) const
{
    const uint64_t k = assemblerInfo->k;

    const auto read = reads->getRead(uint32_t(readId));
    const OrientedReadId orientedReadId0(readId, 0);
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];
    const uint64_t readMarkerCount = orientedReadMarkers0.size();
    SHASTA_ASSERT(kmers1.size() == readMarkerCount);

    // Loop over all markers.
    for(uint64_t ordinal0=0; ordinal0<readMarkerCount; ordinal0++) {
        Kmer kmer0;
        extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
        const Kmer kmer1 = kmer0.reverseComplement(k);
        const uint64_t ordinal1 = readMarkerCount - 1 - ordinal0;
        kmers1[ordinal1] = kmer1;
    }



#if 1
    // Check against the KmerIds stored in the markers.
    // These will soon go away.
    const OrientedReadId orientedReadId1(readId, 1);
    const auto orientedReadMarkers1 = markers[orientedReadId1.getValue()];
    for(uint64_t ordinal1=0; ordinal1<readMarkerCount; ordinal1++) {
        SHASTA_ASSERT(kmers1[ordinal1].id(k) == orientedReadMarkers1[ordinal1].kmerId);
    }
#endif
}



// Get all marker KmerIds for an oriented read.
void Assembler::getOrientedReadMarkerKmerIds(
    OrientedReadId orientedReadId,
    const span<KmerId>& kmerIds) const
{
    const ReadId readId = orientedReadId.getReadId();
    const Strand strand = orientedReadId.getStrand();

    if(strand == 0) {
        getOrientedReadMarkerKmerIdsStrand0(readId, kmerIds);
    } else {
        getOrientedReadMarkerKmerIdsStrand1(readId, kmerIds);
    }
}



void Assembler::getOrientedReadMarkerKmerIdsStrand0(ReadId readId, const span<KmerId>& kmerIds0) const
{
    const uint64_t k = assemblerInfo->k;

    const auto read = reads->getRead(uint32_t(readId));
    const OrientedReadId orientedReadId0(readId, 0);
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];
    const uint64_t readMarkerCount = orientedReadMarkers0.size();
    SHASTA_ASSERT(kmerIds0.size() == readMarkerCount);

    // Loop over all markers.
    for(uint64_t ordinal0=0; ordinal0<readMarkerCount; ordinal0++) {
        Kmer kmer0;
        extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
        kmerIds0[ordinal0] = KmerId(kmer0.id(k));
    }



#if 1
    // Check against the KmerIds stored in the markers.
    // These will soon go away.
    for(uint64_t ordinal0=0; ordinal0<readMarkerCount; ordinal0++) {
        SHASTA_ASSERT(kmerIds0[ordinal0] == orientedReadMarkers0[ordinal0].kmerId);
    }
#endif
}



void Assembler::getOrientedReadMarkerKmerIdsStrand1(ReadId readId, const span<KmerId>& kmerIds1) const
{
    const uint64_t k = assemblerInfo->k;

    const auto read = reads->getRead(uint32_t(readId));
    const OrientedReadId orientedReadId0(readId, 0);
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];
    const uint64_t readMarkerCount = orientedReadMarkers0.size();
    SHASTA_ASSERT(kmerIds1.size() == readMarkerCount);

    // Loop over all markers.
    for(uint64_t ordinal0=0; ordinal0<readMarkerCount; ordinal0++) {
        Kmer kmer0;
        extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
        const Kmer kmer1 = kmer0.reverseComplement(k);
        const uint64_t ordinal1 = readMarkerCount - 1 - ordinal0;
        kmerIds1[ordinal1] = KmerId(kmer1.id(k));
    }



#if 1
    // Check against the KmerIds stored in the markers.
    // These will soon go away.
    const OrientedReadId orientedReadId1(readId, 1);
    const auto orientedReadMarkers1 = markers[orientedReadId1.getValue()];
    for(uint64_t ordinal1=0; ordinal1<readMarkerCount; ordinal1++) {
        SHASTA_ASSERT(kmerIds1[ordinal1] == orientedReadMarkers1[ordinal1].kmerId);
    }
#endif
}



// Get all marker Kmers for a read in both orientations.
void Assembler::getReadMarkerKmers(
    ReadId readId,
    const span<Kmer>& kmers0,
    const span<Kmer>& kmers1) const
{
    const uint64_t k = assemblerInfo->k;

    // Access the information we need for this read.
    const auto read = reads->getRead(uint32_t(readId));
    const OrientedReadId orientedReadId0(uint32_t(readId), 0);
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];
    const uint64_t readMarkerCount = orientedReadMarkers0.size();
    SHASTA_ASSERT(kmers0.size() == readMarkerCount);
    SHASTA_ASSERT(kmers1.size() == readMarkerCount);

    // Loop over all markers.
    for(uint64_t ordinal0=0; ordinal0<readMarkerCount; ordinal0++) {

        // Strand 0.
        Kmer kmer0;
        extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
        kmers0[ordinal0] = kmer0;

        // Strand 1.
        const Kmer kmer1 = kmer0.reverseComplement(k);
        const uint64_t ordinal1 = readMarkerCount - 1 - ordinal0;
        kmers1[ordinal1] = kmer1;
    }



#if 1
    // Check against the KmerIds stored in the markers.
    // These will soon go away.
    for(uint64_t ordinal0=0; ordinal0<readMarkerCount; ordinal0++) {
        SHASTA_ASSERT(kmers0[ordinal0].id(k) == orientedReadMarkers0[ordinal0].kmerId);
    }
    const OrientedReadId orientedReadId1(uint32_t(readId), 1);
    const auto orientedReadMarkers1 = markers[orientedReadId1.getValue()];
    for(uint64_t ordinal1=0; ordinal1<readMarkerCount; ordinal1++) {
        SHASTA_ASSERT(kmers1[ordinal1].id(k) == orientedReadMarkers1[ordinal1].kmerId);
    }
#endif
}



// Get all marker KmerIds for a read in both orientations.
void Assembler::getReadMarkerKmerIds(
    ReadId readId,
    const span<KmerId>& kmerIds0,
    const span<KmerId>& kmerIds1) const
{
    // Get the marker length.
    const uint64_t k = assemblerInfo->k;

    // Access the information we need for this read.
    const auto read = reads->getRead(uint32_t(readId));
    const OrientedReadId orientedReadId0(uint32_t(readId), 0);
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];
    const uint64_t readMarkerCount = orientedReadMarkers0.size();
    SHASTA_ASSERT(kmerIds0.size() == readMarkerCount);
    SHASTA_ASSERT(kmerIds1.size() == readMarkerCount);

    // Loop over all markers.
    for(uint64_t ordinal0=0; ordinal0<readMarkerCount; ordinal0++) {

        // Strand 0.
        Kmer kmer0;
        extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
        kmerIds0[ordinal0] = KmerId(kmer0.id(k));

        // Strand 1.
        const Kmer kmer1 = kmer0.reverseComplement(k);
        const uint64_t ordinal1 = readMarkerCount - 1 - ordinal0;
        kmerIds1[ordinal1] = KmerId(kmer1.id(k));
    }



#if 1
    // Check against the KmerIds stored in the markers.
    // These will soon go away.
    for(uint64_t ordinal0=0; ordinal0<readMarkerCount; ordinal0++) {
        SHASTA_ASSERT(kmerIds0[ordinal0] == orientedReadMarkers0[ordinal0].kmerId);
    }
    const OrientedReadId orientedReadId1(uint32_t(readId), 1);
    const auto orientedReadMarkers1 = markers[orientedReadId1.getValue()];
    for(uint64_t ordinal1=0; ordinal1<readMarkerCount; ordinal1++) {
        SHASTA_ASSERT(kmerIds1[ordinal1] == orientedReadMarkers1[ordinal1].kmerId);
    }
#endif
}



// Get the Kmer for an oriented read at a given marker ordinal.
Kmer Assembler::getOrientedReadMarkerKmer(OrientedReadId orientedReadId, uint64_t ordinal) const
{
    const uint64_t k = assemblerInfo->k;

    const ReadId readId = orientedReadId.getReadId();
    const Strand strand = orientedReadId.getStrand();
    const auto read = reads->getRead(readId);
    const OrientedReadId orientedReadId0(uint32_t(readId), 0);
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];

    if(strand == 0) {

        const uint64_t ordinal0 = ordinal;
        Kmer kmer0;
        extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
        return kmer0;

    } else {

        const uint64_t ordinal0 = orientedReadMarkers0.size() - 1 - ordinal;
        Kmer kmer0;
        extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
        return kmer0.reverseComplement(k);

    }
}



// Get the KmerId for an oriented read at a given marker ordinal.
KmerId Assembler::getOrientedReadMarkerKmerId(OrientedReadId orientedReadId, uint64_t ordinal) const
{
    const uint64_t k = assemblerInfo->k;

    const ReadId readId = orientedReadId.getReadId();
    const Strand strand = orientedReadId.getStrand();
    const auto read = reads->getRead(readId);
    const OrientedReadId orientedReadId0(uint32_t(readId), 0);
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];

    if(strand == 0) {

        const uint64_t ordinal0 = ordinal;
        Kmer kmer0;
        extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
        return KmerId(kmer0.id(k));

    } else {

        const uint64_t ordinal0 = orientedReadMarkers0.size() - 1 - ordinal;
        Kmer kmer0;
        extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
        return KmerId(kmer0.reverseComplement(k).id(k));

    }
}

