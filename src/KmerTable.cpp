// Shasta.
#include "KmerTable.hpp"
#include "AssemblerOptions.hpp"
#include "deduplicate.hpp"
#include "Kmer.hpp"
#include "Reads.hpp"
using namespace shasta;

// Standard library.
#include "fstream.hpp"
#include <random>

// Explicit template instantiations.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<KmerTable1>;
template class MultithreadedObject<KmerTable2>;
template class MultithreadedObject<KmerTable4>;



// Randomly select the k-mers to be used as markers.
KmerTable0::KmerTable0(
    uint64_t k,
    double probability, // The probability that a k-mer is selected as a marker.
    int seed ,           // For random number generator.
    const MappedMemoryOwner& mappedMemoryOwner
    ) : KmerTable(k, true, mappedMemoryOwner)
{

    // The total number of k-mers of this length.
    // This includes both RLE and non-RLE k-mers.
    const size_t kmerCount = 1ULL << (2ULL*k);

    // Sanity check on the requested fraction.
    // It can be 1 at most. If it is 1, all k-mers
    // are guaranteed to be selected (because we use <=
    // comparison in the main loop.
    if(probability<0. || probability>1.) {
        throw runtime_error("Invalid k-mer probability " +
            to_string(probability) + " requested.");
    }



    // Compute the probability p with which we select
    // each k-mer and its reverse complement
    // in order to achieve the required k-mer fraction.
    // For k-mers that are not self-complementary,
    // the probability of not being selected
    // is (1-p)^2 (there are two chances,
    // with and without reverse complement).
    // So, probability = 1 - (1-p)^2, and therefore p=1-sqrt(1-probability).
    // For simplicity, we use the same p for k-mers that are
    // self-complementary. They are a small minority, and because
    // of this they are chose with lower probability.
    // Probably a good thing anyway.
    const double p = 1. - sqrt(1. - probability);
    if(probability == 1.) {
        SHASTA_ASSERT(p == 1.);
    }



    // Prepare to generate uniformly distributed numbers between 0 and 1.
    std::mt19937 randomSource(seed);
    std::uniform_real_distribution<> uniformDistribution;

    // Pick each k-mer and its reverse complement with probability p.
    // Use <= comparison, so if probability=1, p=1, all k-mers are kept.
    for(uint64_t kmerId=0; kmerId<kmerCount; kmerId++) {
        const double x = uniformDistribution(randomSource);
        if(x <= p) {
            kmerTable[kmerId].isMarker = true;
            const uint64_t reverseComplementedKmerId = kmerTable[kmerId].reverseComplementedKmerId;
            kmerTable[reverseComplementedKmerId].isMarker = true;
        }
    }


}



KmerTable::KmerTable(
    uint64_t k,
    bool createNew,
    const MappedMemoryOwner& mappedMemoryOwner) :
    MappedMemoryOwner(mappedMemoryOwner), k(k)
{
    if(createNew) {
        createKmerTable();
    } else {
        accessKmerTable();
    }
}



void KmerTable::createKmerTable()
{
    SHASTA_ASSERT(k <= Kmer16::capacity);

    // Create the kmer table with the necessary size.
    kmerTable.createNew(largeDataName("Kmers"), largeDataPageSize);
    const size_t kmerCount = 1ULL << (2ULL*k);
    kmerTable.resize(kmerCount);

    // Store the reverse complement of each k-mer.
    for(uint64_t kmerId=0; kmerId<kmerCount; kmerId++) {
        const Kmer kmer(kmerId, k);
        const Kmer reverseComplementedKmer = kmer.reverseComplement(k);
        kmerTable[kmerId].reverseComplementedKmerId = KmerId16(reverseComplementedKmer.id(k));
    }
    for(uint64_t kmerId=0; kmerId<kmerCount; kmerId++) {
        const uint64_t reverseComplementedKmerId = kmerTable[kmerId].reverseComplementedKmerId;
        SHASTA_ASSERT(kmerTable[reverseComplementedKmerId].reverseComplementedKmerId == kmerId);
    }



    // Figure out which k-mers are run-length-encoded sequence.
    // They are the ones without repeated bases.
    for(uint64_t kmerId=0; kmerId<kmerCount; kmerId++) {
        kmerTable[kmerId].isRleKmer = true;
        const Kmer kmer(kmerId, k);
        for(size_t i=1; i<k; i++) {
            if(kmer[i-1] == kmer[i]) {
                kmerTable[kmerId].isRleKmer = false;
                break;
            }
        }
    }

}



void KmerTable::accessKmerTable()
{
    kmerTable.accessExistingReadOnly(largeDataName("Kmers"));
    SHASTA_ASSERT(kmerTable.size() == 1ULL << (2ULL*k));
}



// Select marker k-mers randomly, but excluding
// the ones that have high frequency in the reads.
KmerTable1::KmerTable1(

    uint64_t k,

    // The desired marker density
    double markerDensity,

    // Seed for random number generator.
    int seed,

    // Exclude k-mers enriched by more than this amount.
    // Enrichment is the ratio of k-mer frequency in reads
    // over what a random distribution would give.
    double enrichmentThreshold,

    const Reads& reads,

    size_t threadCount,

    const MappedMemoryOwner& mappedMemoryOwner) :

    KmerTable(k, true, mappedMemoryOwner),
    MultithreadedObject<KmerTable1>(*this),
    reads(reads)
{

    // Sanity check.
    if(markerDensity<0. || markerDensity>1.) {
        throw runtime_error("Invalid marker density " +
            to_string(markerDensity) + " requested.");
    }

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Compute the frequency of all k-mers in oriented reads.
    setupLoadBalancing(reads.readCount(), 1000);
    runThreads(&KmerTable1::computeKmerFrequency, threadCount);

    // Compute the total number of k-mer occurrences in reads
    // and the number of k-mers that can possibly occur.
    // This is done by counting all k-mers
    // when using the raw read representation and
    // only RLE k-mers when using the RLE read representation.
    uint64_t totalKmerOccurrences = 0;
    uint64_t possibleKmerCount = 0;
    for(uint64_t kmerId=0; kmerId!=kmerTable.size(); kmerId++) {
        const KmerInfo& info = kmerTable[kmerId];
        totalKmerOccurrences += info.frequency;
        if(reads.representation == 0) {
            ++possibleKmerCount;
        } else {
            if(info.isRleKmer) {
                ++possibleKmerCount;
            }
        }
    }
    const double averageOccurrenceCount =
        double(totalKmerOccurrences) / double(possibleKmerCount);



    if(reads.representation == 0) {

        // We are using raw read representation.
        cout <<
            "K-mer length k " << k << "\n"
            "Total number of k-mers " << kmerTable.size() << "\n"
            "Total number of k-mer occurrences in all oriented reads " << totalKmerOccurrences << "\n"
            "Average number of occurrences per k-mer " <<
            averageOccurrenceCount << endl;

    } else {

        // We are using RLE read representation.
        cout <<
            "K-mer length k " << k << "\n"
            "Total number of k-mers " << kmerTable.size() << "\n"
            "Total number of RLE k-mers " << possibleKmerCount << "\n"
            "Total number of k-mer occurrences in all oriented reads " << totalKmerOccurrences << "\n"
            "Average number of occurrences per RLE k-mer " <<
            averageOccurrenceCount << endl;
    }



    // Convert the enrichment threshold to a frequency.
    const uint64_t frequencyThreshold =
        uint64_t(enrichmentThreshold * averageOccurrenceCount);



    // Write out what we found.
    ofstream csv("KmerFrequencies.csv");
    csv << "KmerId,Kmer,ReverseComplementedKmerId,ReverseComplementedKmer,Frequency,Enrichment,Overenriched?\n";
    for(uint64_t kmerId=0; kmerId<kmerTable.size(); kmerId++) {
        const KmerInfo& info = kmerTable[kmerId];
        const uint64_t frequency = info.frequency;

        const Kmer kmer(kmerId, k);
        const Kmer reverseComplementedKmer(info.reverseComplementedKmerId, k);
        csv << kmerId << ",";
        kmer.write(csv, k);
        csv << ",";
        reverseComplementedKmer.write(csv, k);
        csv << ",";
        csv << info.reverseComplementedKmerId << ",";
        csv << frequency << ",";
        csv << double(frequency) / averageOccurrenceCount;
        csv << ",";
        if(frequency > frequencyThreshold) {
            csv << "Yes";
        } else {
            csv << "No";
        }
        csv << "\n";
    }


    // Gather k-mers that are not overenriched.
    vector<KmerId> candidateKmers;
    for(uint64_t kmerId=0; kmerId<kmerTable.size(); kmerId++) {
        const KmerInfo& info = kmerTable[kmerId];
        if((reads.representation==1) and  (not info.isRleKmer)) {
            continue;
        }
        const uint64_t frequency = info.frequency;
        if(frequency > frequencyThreshold) {
            continue;
        }
        candidateKmers.push_back(KmerId(kmerId));
    }
    cout << possibleKmerCount - candidateKmers.size() << " k-mers were found to be "
        "enriched by more than a factor of " << enrichmentThreshold <<
        " and will not be used as markers." << endl;
    cout << "Markers will be chosen randomly among the remaining pool of " <<
        candidateKmers.size() << " k-mers." << endl;



    // Prepare to generate a random index into this vector of candidate k-mers.
    std::mt19937 randomSource(seed);
    std::uniform_int_distribution<uint64_t> uniformDistribution(0, candidateKmers.size()-1);

    // Flag all k-mers as not markers.
    for(uint64_t kmerId=0; kmerId<kmerTable.size(); kmerId++) {
        kmerTable[kmerId].isMarker = false;
    }


    // Now randomly generate markers from this pool of k-mers
    // until we have enough.
    uint64_t kmerOccurrencesCount = 0;
    uint64_t kmerCount = 0;
    const uint64_t giveUpCount =  uint64_t(0.9 * double(candidateKmers.size()));
    const uint64_t desiredKmerOccurrencesCount =
        uint64_t(markerDensity * double(totalKmerOccurrences));
    while(kmerOccurrencesCount < desiredKmerOccurrencesCount) {

        // Generate a random index into the candidateKmers vector.
        const uint64_t index = uniformDistribution(randomSource);

        // Check that this k-mer is not already selected as a marker.
        const KmerId kmerId = candidateKmers[index];
        KmerInfo& info = kmerTable[uint64_t(kmerId)];
        if(info.isMarker) {
            continue;
        }

        // This k-mer is not already selected as a marker.
        // Let's add it.
        info.isMarker = true;
        kmerOccurrencesCount += info.frequency;
        ++kmerCount;

        // If this k-mer is palindromic, we are done.
        if(info.reverseComplementedKmerId == kmerId) {
            continue;
        }

        // This k-mer is not palindromic, so we also add its reverse complement.
        KmerInfo& reverseComplementedInfo = kmerTable[info.reverseComplementedKmerId];
        SHASTA_ASSERT(!reverseComplementedInfo.isMarker);
        SHASTA_ASSERT(reverseComplementedInfo.frequency == info.frequency);
        reverseComplementedInfo.isMarker = true;
        kmerOccurrencesCount += reverseComplementedInfo.frequency;
        ++kmerCount;

        if(kmerCount >= giveUpCount) {
            throw runtime_error("Giving up after selecting as markers 90% of the candidate kmers.");
        }
    }
    cout << "Selected " << kmerCount << " k-mers as markers." << endl;


}



void KmerTable1::computeKmerFrequency(size_t threadId)
{
    // Create a frequency vector for this thread.
    MemoryMapped::Vector<uint64_t> frequency;
    frequency.createNew(
        largeDataName("tmp-KmerFrequency-" + to_string(threadId)),
        largeDataPageSize);
    frequency.resize(kmerTable.size());
    fill(frequency.begin(), frequency.end(), 0);



    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all reads assigned to this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); readId++) {

            // Access the sequence of this read.
            const LongBaseSequenceView read = reads.getRead(readId);

            // If the read is pathologically short, it has no k-mers.
            if(read.baseCount < k) {
                continue;
            }

            // Loop over k-mers of this read.
            Kmer kmer;
            for(size_t position=0; position<k; position++) {
                kmer.set(position, read[position]);
            }
            for(uint32_t position=0; /*The check is done later */; position++) {

                // Get the k-mer id.
                const KmerId kmerId = KmerId(kmer.id(k));

                // Increment its frequency.
                ++frequency[uint64_t(kmerId)];

                // Also increment the frequency of the reverse complemented k-mer.
                ++frequency[uint64_t(kmerTable[uint64_t(kmerId)].reverseComplementedKmerId)];

                // Check if we reached the end of the read.
                if(position+k == read.baseCount) {
                    break;
                }

                // Update the k-mer.
                kmer.shiftLeft();
                kmer.set(k-1, read[position+k]);
            }
        }
    }


    // Update the frequency in the k-mer table.
    {
        std::lock_guard<std::mutex> lock(mutex);
        for(uint64_t kmerId=0; kmerId!=frequency.size(); kmerId++) {
            kmerTable[kmerId].frequency += frequency[kmerId];
        }
    }



    // Remove the frequency vector for this thread.
    frequency.remove();
}



// Read the k-mers from file.
KmerTable3::KmerTable3(
    uint64_t k,
    uint64_t readRepresentation,
    const string& fileName,
    const MappedMemoryOwner& mappedMemoryOwner) :
    KmerTable(k, true, mappedMemoryOwner)
{
    if(fileName.empty() or
        fileName[0] != '/') {
        throw runtime_error("Option --Kmers.file must specify an absolute path. "
            "A relative path is not accepted.");
    }

    // Open the file.
    ifstream file(fileName);
    if(not file) {
        throw runtime_error("Error opening " + fileName);
    }



    // Read one k-mer per line.
    uint64_t lineCount = 0;
    while(true) {

        // Read a line.
        string line;
        std::getline(file, line);
        if(not file) {
            break;
        }

        // Check the length.
        if(line.size() != k) {
            throw runtime_error("Unexpected line length in " + fileName + ":\n" + line + "\n" +
                "Expected " + to_string(k) + " characters, found " + to_string(line.size()));
        }

        // Read the k-mer.
        Kmer kmer;
        for(uint64_t i=0; i<k; i++) {
            const char c = line[i];
            const Base base = Base::fromCharacterNoException(c);
            if(not base.isValid()) {
                throw runtime_error("Unexpected base character in " + fileName + ":\n" + line);
            }
            kmer.set(i, base);
        }

        // Sanity checks.
        const KmerId kmerId = KmerId(kmer.id(k));
        SHASTA_ASSERT(kmerId < kmerTable.size());
        KmerInfo& kmerInfo = kmerTable[uint64_t(kmerId)];
        if((readRepresentation==1) and (not kmerInfo.isRleKmer)) {
            throw runtime_error("Non-RLE k-mer (duplicate consecutive bases) in " +
                fileName + ":\n" + line);
        }

        // Flag it as a marker, together with its reverse complement.
        kmerInfo.isMarker = 1;
        kmerTable[kmerInfo.reverseComplementedKmerId].isMarker = 1;
        ++lineCount;

    }
    cout << "Processed " << lineCount << " lines of " << fileName << endl;

    // Count the number of k-mers flagged as markers.
    uint64_t usedKmerCount = 0;
    uint64_t possibleKmerCount = 0;
    for(const KmerInfo& kmerInfo: kmerTable) {
        if(kmerInfo.isMarker) {
            ++usedKmerCount;
        }
        if(readRepresentation == 0) {
            ++possibleKmerCount;
        } else {
            if(kmerInfo.isRleKmer) {
                ++possibleKmerCount;
            }
        }
    }
    cout << "Flagged as markers " << usedKmerCount << " out of " << possibleKmerCount <<
        " possible k-mers of length " << k << endl;
}



// In this version, marker k-mers are selected randomly, but excluding
// any k-mer that is over-enriched even in a single oriented read.
KmerTable2::KmerTable2(
    uint64_t k,
    double markerDensity,
    int seed,
    double enrichmentThreshold,
    const Reads& reads,
    uint64_t threadCount,
    const MappedMemoryOwner& mappedMemoryOwner) :
    KmerTable(k, true, mappedMemoryOwner),
    MultithreadedObject<KmerTable2>(*this),
    reads(reads),
    enrichmentThreshold(enrichmentThreshold)
{

    // Sanity check.
    if(markerDensity<0. || markerDensity>1.) {
        throw runtime_error("Invalid marker density " +
            to_string(markerDensity) + " requested.");
    }

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // For each KmerId that is an RLE k-mer, compute the
    // global frequency (total number of occurrences in all
    // oriented reads) and the number of reads in
    // which the k-mer is over-enriched.
    globalFrequency.createNew(
        largeDataName("tmp-SelectKmers2-GlobalFrequency"),  largeDataPageSize);
    overenrichedReadCount.createNew(
        largeDataName("tmp-SelectKmers2-OverenrichedReadCount"),  largeDataPageSize);
    globalFrequency.resize(kmerTable.size());
    overenrichedReadCount.resize(kmerTable.size());
    fill(
        globalFrequency.begin(),
        globalFrequency.end(), 0);
    fill(
        overenrichedReadCount.begin(),
        overenrichedReadCount.end(), 0);
    setupLoadBalancing(reads.readCount(), 100);
    runThreads(&KmerTable2::threadFunction, threadCount);



    // Compute the total number of k-mer occurrences
    // and the number of possible kmers.
    uint64_t totalKmerOccurrences = 0;
    uint64_t possibleKmerCount = 0;
    for(uint64_t kmerId=0; kmerId!=kmerTable.size(); kmerId++) {
        totalKmerOccurrences += globalFrequency[kmerId];
        if(reads.representation == 0) {
            ++ possibleKmerCount;
        } else {
            if(kmerTable[kmerId].isRleKmer) {
                ++possibleKmerCount;
            }
        }
    }
    const double averageOccurrenceCount =
        double(totalKmerOccurrences) / double(possibleKmerCount);



    // Write out what we found.
    ofstream csv("KmerFrequencies.csv");
    csv << "KmerId,Kmer,ReverseComplementedKmerId,ReverseComplementedKmer,"
        "GlobalFrequency,GlobalEnrichment,NumberOfReadsOverenriched\n";
    for(uint64_t kmerId=0; kmerId<kmerTable.size(); kmerId++) {
        const KmerInfo& info = kmerTable[kmerId];
        const uint64_t frequency = globalFrequency[kmerId];

        const Kmer kmer(kmerId, k);
        const Kmer reverseComplementedKmer(info.reverseComplementedKmerId, k);
        csv << kmerId << ",";
        kmer.write(csv, k);
        csv << ",";
        reverseComplementedKmer.write(csv, k);
        csv << ",";
        csv << info.reverseComplementedKmerId << ",";
        csv << frequency << ",";
        csv << double(frequency) / averageOccurrenceCount;
        csv << ",";
        csv << overenrichedReadCount[kmerId];

        csv << "\n";
    }
    csv.close();



    // Gather k-mers that are not overenriched in any read and therefore
    // can be used as markers.
    vector<KmerId> candidateKmers;
    for(uint64_t kmerId=0; kmerId<kmerTable.size(); kmerId++) {
        const bool readIsUsable = (reads.representation==0) ? true : kmerTable[kmerId].isRleKmer;
        if(readIsUsable and overenrichedReadCount[kmerId] == 0) {
            candidateKmers.push_back(KmerId(kmerId));
        }
    }
    cout << "Out of a total " << possibleKmerCount << " possible k-mers, " <<
        possibleKmerCount - candidateKmers.size() <<
        " were found to be over-enriched by more than a factor of " <<
        enrichmentThreshold <<
        " in at least one read and will not be used as markers." << endl;
    cout << "Markers will be chosen randomly from the remaining pool of " <<
        candidateKmers.size() << " k-mers." << endl;
    cout << "The enrichment threshold of " << enrichmentThreshold <<
        " corresponds to one occurrence every " <<
        double(possibleKmerCount) / enrichmentThreshold <<
        " bases";
    if(reads.representation == 1) {
        cout << " (in RLE representation)";
    }
    cout << "." << endl;


    // Prepare to generate a random index into this vector of candidate k-mers.
    std::mt19937 randomSource(seed);
    std::uniform_int_distribution<uint64_t> uniformDistribution(0, candidateKmers.size()-1);

    // Flag all k-mers as not markers.
    for(uint64_t kmerId=0; kmerId<kmerTable.size(); kmerId++) {
        kmerTable[kmerId].isMarker = false;
    }



    // Now randomly generate markers from this pool of k-mers
    // until we have enough.
    uint64_t kmerOccurrencesCount = 0;
    uint64_t kmerCount = 0;
    const uint64_t desiredKmerOccurrencesCount =
        uint64_t(markerDensity * double(totalKmerOccurrences));
    const uint64_t giveUpCount =  uint64_t(0.9 * double(candidateKmers.size()));
    while(kmerOccurrencesCount < desiredKmerOccurrencesCount) {

        // Generate a random index into the candidateKmers vector.
        const uint64_t index = uniformDistribution(randomSource);

        // Check that this k-mer is not already selected as a marker.
        const KmerId kmerId = candidateKmers[index];
        KmerInfo& info = kmerTable[uint64_t(kmerId)];
        if(info.isMarker) {
            continue;
        }

        // This k-mer is not already selected as a marker.
        // Let's add it.
        info.isMarker = true;
        kmerOccurrencesCount += globalFrequency[uint64_t(kmerId)];
        ++kmerCount;

        // If this k-mer is palindromic, we are done.
        if(info.reverseComplementedKmerId == kmerId) {
            continue;
        }

        // This k-mer is not palindromic, so we also add its reverse complement.
        KmerInfo& reverseComplementedInfo = kmerTable[info.reverseComplementedKmerId];
        SHASTA_ASSERT(!reverseComplementedInfo.isMarker);
        SHASTA_ASSERT(reverseComplementedInfo.frequency == info.frequency);
        reverseComplementedInfo.isMarker = true;
        kmerOccurrencesCount += globalFrequency[info.reverseComplementedKmerId];
        ++kmerCount;

        if(kmerCount >= giveUpCount) {
            throw runtime_error("Giving up after selecting as markers 90% of the candidate kmers.");
        }
    }
    cout << "Selected " << kmerCount << " k-mers as markers." << endl;
    cout << "These k-mers have a total " << kmerOccurrencesCount <<
        " occurrences out of a total " << totalKmerOccurrences <<
        " in all oriented reads." << endl;

}



void KmerTable2::threadFunction(size_t threadId)
{
    // Initialize globalFrequency for this thread.
    MemoryMapped::Vector<uint64_t> threadGlobalFrequency;
    threadGlobalFrequency.createNew(
        largeDataName("tmp-KmerTable2-GlobalFrequency-" + to_string(threadId)),
        largeDataPageSize);
    threadGlobalFrequency.resize(kmerTable.size());
    fill(threadGlobalFrequency.begin(), threadGlobalFrequency.end(), 0);

    // Initialize overenrichedReadCount for this thread.
    MemoryMapped::Vector<ReadId> threadOverenrichedReadCount;
    threadOverenrichedReadCount.createNew(
        largeDataName("tmp-KmerTable2-OverenrichedReadCount-" + to_string(threadId)),
        largeDataPageSize);
    threadOverenrichedReadCount.resize(kmerTable.size());
    fill(threadOverenrichedReadCount.begin(), threadOverenrichedReadCount.end(), 0);

    // Vectors to hold KmerIds and their frequencies for a single read.
    vector<KmerId> readKmerIds;
    vector<uint32_t> readKmerIdFrequencies;

    // Compute the total number of possible k-mers.
    // It is needed below for overenrichment computations.
    uint64_t possibleKmerCount = 0;
    for(const KmerInfo& kmerInfo: kmerTable) {
        if(reads.representation == 0) {
            ++possibleKmerCount;
        } else {
            if(kmerInfo.isRleKmer) {
                ++possibleKmerCount;
            }
        }
    }


    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all reads assigned to this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); readId++) {

            // Access the sequence of this read.
            const LongBaseSequenceView read = reads.getRead(readId);

            // If the read is pathologically short, it has no k-mers.
            if(read.baseCount < k) {
                continue;
            }

            // Loop over k-mers of this read.
            Kmer kmer;
            for(size_t position=0; position<k; position++) {
                kmer.set(position, read[position]);
            }
            readKmerIds.clear();
            for(uint32_t position=0; /*The check is done later */; position++) {

                // Get the k-mer id.
                const KmerId kmerId = KmerId(kmer.id(k));
                readKmerIds.push_back(kmerId);

                // Increment its global frequency.
                ++threadGlobalFrequency[uint64_t(kmerId)];

                // Also increment the frequency of the reverse complemented k-mer.
                ++threadGlobalFrequency[uint64_t(kmerTable[uint64_t(kmerId)].reverseComplementedKmerId)];

                // Check if we reached the end of the read.
                if(position+k == read.baseCount) {
                    break;
                }

                // Update the k-mer.
                kmer.shiftLeft();
                kmer.set(k-1, read[position+k]);
            }

            // Compute k-mer frequencies for this read.
            deduplicateAndCount(readKmerIds, readKmerIdFrequencies);

            // Compute the k-mer frequency threshold for an over-enriched k-mer
            // in this read.
            const uint32_t readKmerCount = uint32_t(read.baseCount + 1 - k);
            const uint32_t frequencyThreshold =
                uint32_t(enrichmentThreshold * double(readKmerCount) / double(possibleKmerCount));

            // See if any k-mers are over-enriched in this read.
            SHASTA_ASSERT(readKmerIds.size() == readKmerIdFrequencies.size());
            for(uint64_t i=0; i<readKmerIds.size(); i++) {
                const KmerId kmerId = readKmerIds[i];
                const uint32_t frequency = readKmerIdFrequencies[i];
                if(frequency > frequencyThreshold) {
                    ++threadOverenrichedReadCount[uint64_t(kmerId)];
                    ++threadOverenrichedReadCount[uint64_t(kmerTable[uint64_t(kmerId)].reverseComplementedKmerId)];
                }
            }
        }
    }



    // Add our globalFrequency and overenrichedReadCount
    // to the values computer by the other threads.
    {
        std::lock_guard<std::mutex> lock(mutex);
        for(uint64_t kmerId=0; kmerId!=globalFrequency.size(); kmerId++) {
            globalFrequency[kmerId] += threadGlobalFrequency[kmerId];
            overenrichedReadCount[kmerId] += threadOverenrichedReadCount[kmerId];
        }
    }
    threadGlobalFrequency.remove();
    threadOverenrichedReadCount.remove();
}



// In this version, marker k-mers are selected randomly, but excluding
// k-mers that appear repeated at short distances in any oriented read.
// More precisely, for each k-mer we compute the minimum distance
// (in RLE bases) at which any two copies of that k-mer appear in any oriented read.
// K-mers for which this minimum distance is less than distanceThreshold
// are not used as markers. Marker k-mers are selected randomly among the
// remaining k-mers, until the desired marker density is achieved.
KmerTable4::KmerTable4(
    uint64_t k,
    double markerDensity,
    int seed,
    uint64_t distanceThreshold,
    const Reads& reads,
    uint64_t threadCount,
    const MappedMemoryOwner& mappedMemoryOwner) :
        KmerTable(k, true, mappedMemoryOwner),
    MultithreadedObject<KmerTable4>(*this),
    reads(reads)
{
    const bool debug = false;

    // Sanity check.
    if(markerDensity<0. || markerDensity>1.) {
        throw runtime_error("Invalid marker density " +
            to_string(markerDensity) + " requested.");
    }

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Initialize the global frequency of all k-mers.
    globalFrequency.createNew(
        largeDataName("tmp-KmerTable44-GlobalFrequency"),  largeDataPageSize);
    globalFrequency.resize(kmerTable.size());
    fill(
        globalFrequency.begin(),
        globalFrequency.end(), 0);


    // Initialize the minimumDistance vector, which stores
    // the minimum RLE distance between any two copies of each k-mer
    // in any oriented read.
    minimumDistance.createNew(
        largeDataName("tmp-KmerTable4-minimumDistance"), largeDataPageSize);
    const uint64_t kmerCount = kmerTable.size();
    minimumDistance.resize(kmerCount);
    for(uint64_t i=0; i<kmerCount; i++) {
        minimumDistance[i].second = std::numeric_limits<uint32_t>::max();
    }

    // Compute the minimumDistance vector.
    setupLoadBalancing(reads.readCount(), 100);
    runThreads(&KmerTable4::threadFunction, threadCount);



    // Write out what we found.
    if(debug) {
        const uint64_t totalFrequency = std::accumulate(
            globalFrequency.begin(),
            globalFrequency.end(), 0ULL);
        cout << "Total number of k-mer occurrences in all oriented reads is " << totalFrequency << endl;
        ofstream csv("KmerInfo.csv");
        csv << "KmerId,Kmer,KmerIdRc,KmerRc,Frequency,FrequencyRc,TotalFrequency,"
            "MinDist,MinDistRc,MinMinDist\n";
        for(uint64_t kmerId=0; kmerId<kmerTable.size(); kmerId++) {
            const KmerInfo& info = kmerTable[kmerId];
            if(!info.isRleKmer) {
                continue;
            }

            const uint64_t frequency = globalFrequency[kmerId];
            const uint64_t frequencyReverseComplement = globalFrequency[info.reverseComplementedKmerId];
            const uint64_t totalFrequency = frequency + frequencyReverseComplement;

            const uint32_t kmerMinimumDistance = minimumDistance[kmerId].second;
            const uint32_t kmerMinimumDistanceReverseComplement =
                minimumDistance[info.reverseComplementedKmerId].second;

            const Kmer kmer(kmerId, k);
            const Kmer reverseComplementedKmer(info.reverseComplementedKmerId, k);
            csv << kmerId << ",";
            kmer.write(csv, k);
            csv << ",";
            csv << info.reverseComplementedKmerId << ",";
            reverseComplementedKmer.write(csv, k);
            csv << ",";
            csv << frequency << ",";
            csv << frequencyReverseComplement << ",";
            csv << totalFrequency << ",";
            csv << kmerMinimumDistance << ",";
            csv << kmerMinimumDistanceReverseComplement << ",";
            csv << min(kmerMinimumDistance, kmerMinimumDistanceReverseComplement) << "\n";
        }
    }



    // Compute the total number of k-mer occurrences
    // and the number of RLE kmers.
    uint64_t totalKmerOccurrences = 0;
    uint64_t rleKmerCount = 0;
    for(uint64_t kmerId=0; kmerId!=kmerTable.size(); kmerId++) {
        const KmerInfo& info = kmerTable[kmerId];
        if((reads.representation==1) and (not info.isRleKmer)) {
            SHASTA_ASSERT(globalFrequency[kmerId] == 0);
            continue;
        }
        totalKmerOccurrences += globalFrequency[kmerId];
        if(kmerTable[kmerId].isRleKmer) {
            ++rleKmerCount;
        }
    }
    cout << "K-mer length k " << k << endl;
    cout << "Distance threshold " << distanceThreshold << " RLE bases." << endl;
    cout << "Total number of distinct RLE k-mers " << rleKmerCount << endl;
    cout << "Total number of RLE k-mers in all oriented reads " << totalKmerOccurrences << endl;
    cout << "Requested marker density " << markerDensity << endl;
    const uint64_t requiredMarkerOccurrences = uint64_t(markerDensity * double(totalKmerOccurrences));
    cout << "Required number of marker occurrences in all oriented reads " << requiredMarkerOccurrences << endl;



    // Gather k-mers for which the minimum distance between two copies
    // equals at least distanceThreshold. Exclude palindromic k-mers.
    vector<KmerId> candidateKmers;
    uint64_t candidateFrequency = 0;
    for(uint64_t kmerId=0; kmerId<kmerTable.size(); kmerId++) {
        const KmerInfo& info = kmerTable[kmerId];
        const KmerId kmerIdRc = info.reverseComplementedKmerId;
        if((reads.representation==1) and (not info.isRleKmer)) {
            continue;
        }
        if(kmerIdRc == kmerId) {
            // Palindromic. Exclude.
            continue;
        }
        // Only store the lower KmerId in the pair.
        if(kmerId > kmerIdRc) {
            continue;
        }
        if(minimumDistance[uint64_t(kmerId)].second < distanceThreshold) {
            // Too close. skip.
            continue;
        }
        if(minimumDistance[uint64_t(kmerIdRc)].second < distanceThreshold) {
            // Too close. Skip.
            continue;
        }

        candidateKmers.push_back(KmerId(kmerId));
        candidateFrequency += globalFrequency[uint64_t(kmerId)];
        candidateFrequency += globalFrequency[uint64_t(kmerIdRc)];
    }
    cout << "Markers will be chosen randomly from the a pool of " <<
        2*candidateKmers.size() << " RLE k-mers." << endl;
    cout << "RLE k-mers in this pool occur " <<
        candidateFrequency << " times in all oriented reads." << endl;
    cout << "This is sufficient to achieve marker density up to " <<
        double(candidateFrequency) / double(totalKmerOccurrences) << endl;

    // If these candidates don't have sufficient frequency, we
    // can't achieve the required marker density.
    if(candidateFrequency < requiredMarkerOccurrences) {
        throw runtime_error("Cannot achieve required marker density. "
            "Increase k, or decrease marker density, or decrease distance threshold.");
    }

    // Flag all k-mers as not markers.
    for(uint64_t kmerId=0; kmerId<kmerTable.size(); kmerId++) {
        kmerTable[kmerId].isMarker = false;
    }



    // Randomly pick markers in this vector of candidate k-mers.
    std::mt19937_64 randomSource(seed);
    std::uniform_real_distribution<> uniformDistribution;
    uint64_t markerOccurrencesCount = 0;
    uint64_t markerCount = 0;
    while(true) {

        // Pick a random index in the candidateKmers vector.
        const double x = uniformDistribution(randomSource);  // In [0,1)
        const uint64_t i = uint64_t(x * double(candidateKmers.size()));
        SHASTA_ASSERT(i < candidateKmers.size());

        // This KmerId and its reverse complement  will be used as markers.
        const KmerId kmerId = candidateKmers[i];
        const KmerId kmerIdRc = kmerTable[uint64_t(kmerId)].reverseComplementedKmerId;

        kmerTable[uint64_t(kmerId)].isMarker = true;
        kmerTable[uint64_t(kmerIdRc)].isMarker = true;

        // Increment counters.
        markerCount += 2;
        markerOccurrencesCount += globalFrequency[uint64_t(kmerId)];
        markerOccurrencesCount += globalFrequency[uint64_t(kmerIdRc)];

        // Remove kmerId from the vector of candidates.
        if(i != candidateKmers.size()-1) {
            candidateKmers[i] = candidateKmers.back();
        }
        candidateKmers.resize(candidateKmers.size() - 1);

        if(markerOccurrencesCount >= requiredMarkerOccurrences) {
            break;
        }
    }
    cout << "Selected " << markerCount << " k-mers as markers." << endl;
    cout << "Actual marker density " << double(markerOccurrencesCount) / double(totalKmerOccurrences) << endl;



    // Clean up.
    minimumDistance.remove();
    globalFrequency.remove();

}



void KmerTable4::threadFunction(size_t threadId)
{
    // Initialize globalFrequency for this thread.
    // Having all threads accumulate atomically on the global frequency vector is too slow.
    MemoryMapped::Vector<uint64_t> threadGlobalFrequency;
    threadGlobalFrequency.createNew(
        largeDataName("tmp-KmerTable4-GlobalFrequency-" + to_string(threadId)),
        largeDataPageSize);
    threadGlobalFrequency.resize(kmerTable.size());
    fill(threadGlobalFrequency.begin(), threadGlobalFrequency.end(), 0);

    // Vector to hold pairs(KmerId, RLE position) for one read.
    vector< pair<KmerId, uint32_t> > readKmers;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all reads assigned to this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); readId++) {

            // Access the sequence of this read.
            const LongBaseSequenceView read = reads.getRead(readId);

            // If the read is pathologically short, it has no k-mers.
            if(read.baseCount < k) {
                continue;
            }

            // Loop over k-mers of this read.
            Kmer kmer;
            for(size_t position=0; position<k; position++) {
                kmer.set(position, read[position]);
            }
            readKmers.clear();
            for(uint32_t position=0; /*The check is done later */; position++) {

                // Get the k-mer id.
                const KmerId kmerId = KmerId(kmer.id(k));
                readKmers.push_back(make_pair(kmerId, position));

                // Update the frequency of this k-mer.
                ++threadGlobalFrequency[uint64_t(kmerId)];
                ++threadGlobalFrequency[uint64_t(kmerTable[uint64_t(kmerId)].reverseComplementedKmerId)];

                // Check if we reached the end of the read.
                if(position+k == read.baseCount) {
                    break;
                }

                // Update the k-mer.
                kmer.shiftLeft();
                kmer.set(k-1, read[position+k]);
            }

            // Sort by k-mer, then by position.
            sort(readKmers.begin(), readKmers.end());

            // Update minDistance for each pair of repeated k-mers.
            for(uint64_t i=1; i<readKmers.size(); i++) {
                const auto& p0 = readKmers[i-1];
                const auto& p1 = readKmers[i];
                const KmerId kmerId0 = p0.first;
                const KmerId kmerId1 = p1.first;
                if(kmerId0 != kmerId1) {
                    continue;
                }
                const uint32_t distance = p1.second - p0.second;

                pair<std::mutex, uint32_t>& p = minimumDistance[uint64_t(kmerId0)];
                std::lock_guard<std::mutex> lock(p.first);;
                p.second = min(p.second, distance);
            }
        }
    }

    // Add our globalFrequency to the values computed by the other threads.
    {
        std::lock_guard<std::mutex> lock(mutex);
        for(uint64_t kmerId=0; kmerId!=globalFrequency.size(); kmerId++) {
            globalFrequency[kmerId] += threadGlobalFrequency[kmerId];
        }
    }
    threadGlobalFrequency.remove();
}



KmerTable0::KmerTable0(
    uint64_t k,
    const MappedMemoryOwner& mappedMemoryOwner) :
    KmerTable(k, false, mappedMemoryOwner)
{
}



KmerTable1::KmerTable1(
    uint64_t k,
    const Reads& reads,
    const MappedMemoryOwner& mappedMemoryOwner) :
    KmerTable(k, false, mappedMemoryOwner),
    MultithreadedObject<KmerTable1>(*this),
    reads(reads)
{
}



KmerTable2::KmerTable2(
    uint64_t k,
    const Reads& reads,
    const MappedMemoryOwner& mappedMemoryOwner) :
    KmerTable(k, false, mappedMemoryOwner),
    MultithreadedObject<KmerTable2>(*this),
    reads(reads)
{
}



KmerTable3::KmerTable3(
    uint64_t k,
    const MappedMemoryOwner& mappedMemoryOwner) :
    KmerTable(k, false, mappedMemoryOwner)
{
}



KmerTable4::KmerTable4(
    uint64_t k,
    const Reads& reads,
    const MappedMemoryOwner& mappedMemoryOwner) :
    KmerTable(k, false, mappedMemoryOwner),
    MultithreadedObject<KmerTable4>(*this),
    reads(reads)
{
}

