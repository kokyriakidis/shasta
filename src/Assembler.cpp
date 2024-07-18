#include "Assembler.hpp"
#include "AssemblerOptions.hpp"
#include "buildId.hpp"
#include "Coverage.hpp"
#include "KmerCheckerFactory.hpp"
#include "MedianConsensusCaller.hpp"
#include "MurmurHash2.hpp"
#include "Reads.hpp"
#include "SimpleConsensusCaller.hpp"
#include "SimpleBayesianConsensusCaller.hpp"
using namespace shasta;

#include "MultithreadedObject.tpp"
template class MultithreadedObject<Assembler>;


// Constructor to be called one to create a new run.
Assembler::Assembler(
    const string& largeDataFileNamePrefixArgument,
    bool createNew,
    uint64_t readRepresentation, // 0 = raw sequence, 1 = RLE sequence. Only used if createNew.
    size_t largeDataPageSizeArgument) :
    MultithreadedObject(*this)
{
    largeDataFileNamePrefix = largeDataFileNamePrefixArgument;

    if(createNew) {

        // Create a new assembly.
        assemblerInfo.createNew(largeDataName("Info"), largeDataPageSizeArgument);
        assemblerInfo->readRepresentation = readRepresentation;
        assemblerInfo->largeDataPageSize = largeDataPageSizeArgument;
        largeDataPageSize = largeDataPageSizeArgument;

        reads = make_unique<Reads>();
        reads->createNew(
            assemblerInfo->readRepresentation,
            largeDataName("Reads"),
            largeDataName("ReadNames"),
            largeDataName("ReadMetaData"),
            largeDataName("ReadRepeatCounts"),
            largeDataName("ReadFlags"),
            largeDataName("ReadIdsSortedByName"),
            largeDataPageSize
        );
        // cout << "Created a new assembly with page size " << largeDataPageSize << endl;

    } else {

        // Access an existing assembly.
        assemblerInfo.accessExistingReadWrite(largeDataName("Info"));
        largeDataPageSize = assemblerInfo->largeDataPageSize;

        reads = make_unique<Reads>();
        reads->access(
            assemblerInfo->readRepresentation,
            largeDataName("Reads"),
            largeDataName("ReadNames"),
            largeDataName("ReadMetaData"),
            largeDataName("ReadRepeatCounts"),
            largeDataName("ReadFlags"),
            largeDataName("ReadIdsSortedByName")
        );
        // cout << "Accessed an existing assembly with page size " << largeDataPageSize << endl;

    }
    SHASTA_ASSERT(largeDataPageSize == assemblerInfo->largeDataPageSize);

    // In both cases, assemblerInfo, reads, readNames, readRepeatCounts are all open for write.

    fillServerFunctionTable();
}



// Set up the ConsensusCaller used to compute the "best"
// base and repeat count at each assembly position.
// The argument to setupConsensusCaller specifies
// the consensus caller to be used.
// It can be one of the following:
// - Modal
//   Selects the SimpleConsensusCaller.
// - Median
//   Selects the MedianConsensusCaller.
// - Bayesian:fileName
//   Selects the SimpleBayesianConsensusCaller,
//   using fileName as the configuration file.
//   Filename must be an absolute path (it must begin with "/").
void Assembler::setupConsensusCaller(const string& s)
{

    // Parse the argument.
    // The portion up to the first colon is the type string (Modal, Median, or Bayesian).
    // The portion following the first colon is the constructor string.
    const size_t colonPosition = s.find_first_of(':');
    string typeString, constructorString;
    if(colonPosition==string::npos) {
        // There is no colon.
        typeString = s;
        constructorString = "";
    } else {
        // The colon is present.
        typeString =  s.substr(0, colonPosition);
        constructorString =  s.substr(colonPosition+1);
    }



    // The Modal consensus caller (class SimpleConsensusCaller)
    // takes no constructor string (nothing after the colon).
    if(typeString == "Modal") {

        //  The constructor string must be empty.
        if(constructorString.size() > 0) {
            throw runtime_error("Invalid consensus caller " + s);
        }

        consensusCaller = std::make_shared<SimpleConsensusCaller>();
        return;
    }



    // The Median consensus caller (class MedianConsensusCaller)
    // takes no  constructor string (nothing after the colon).
    if(typeString == "Median") {

        //  The constructor string must be empty.
        if(constructorString.size() > 0) {
            throw runtime_error("Invalid consensus caller " + s);
        }

        consensusCaller = std::make_shared<MedianConsensusCaller>();
        return;
    }



    // The Bayesian consensus caller (class SimpleBayesianConsensusCaller)
    // requires a constructor string (portion following the first colon).
    // The constructor string is either the name of a built-in Bayesian model
    // or the absolute path to the configuration file.
    if(typeString == "Bayesian") {
        if(constructorString.size() == 0) {
            throw runtime_error("Invalid consensus caller " + s +
                ". Bayesian:builtinName required or Bayesian::absolutePath");
        }

        consensusCaller = std::make_shared<SimpleBayesianConsensusCaller>(constructorString);
        return;
    }



    // If getting here, the argument does not specify a supported
    // consensus caller.
    throw runtime_error("Invalid consensus caller " + s +
        ". Valid choices are: Modal, Median, Bayesian:absolutePathToConfigFile.");
}



// Store assembly time.
void Assembler::storeAssemblyTime(
    double elapsedTimeSeconds,
    double averageCpuUtilization)
{
    assemblerInfo->assemblyElapsedTimeSeconds = elapsedTimeSeconds;
    assemblerInfo->averageCpuUtilization = averageCpuUtilization;
}

void Assembler::storePeakMemoryUsage(uint64_t peakMemoryUsage) {
    assemblerInfo->peakMemoryUsage = peakMemoryUsage;
}



void Assembler::createKmerChecker(
    const KmersOptions& kmersOptions,
    uint64_t threadCount)
{
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    assemblerInfo->k = kmersOptions.k;
    assemblerInfo->kmerGenerationMethod = kmersOptions.generationMethod;

    kmerChecker = KmerCheckerFactory::createNew(
        kmersOptions,
        threadCount,
        getReads(),
        *this);
}



void Assembler::accessKmerChecker()
{
    kmerChecker = KmerCheckerFactory::createFromBinaryData(
        assemblerInfo->k,
        assemblerInfo->kmerGenerationMethod,
        getReads(),
        *this);
}



// Hash a KmerId in such a way that it has the same hash as its reverse
// complement. This is used by alignment method 3 to downsample markers.
uint32_t Assembler::hashKmerId(KmerId kmerId) const
{
    const uint64_t k = assemblerInfo->k;

    // Construct the k-mer and its reverse complement.
    const Kmer kmer(kmerId, k);
    const Kmer kmerRc = kmer.reverseComplement(k);

    // Compute the id of the reverse complement k-mer.
    const KmerId kmerIdRc = KmerId(kmerRc.id(k));

    // Hash the sum of the two KmerIds.
    // This guarantees that we return the same hash
    // for a k-mer and its reverse complement.
    const uint64_t sum = kmerId + kmerIdRc;

    return MurmurHash2(&sum, sizeof(sum), 13477);
}



void Assembler::createSaveBinaryDataDirectory(const string& memoryMode)
{
    if(memoryMode == "anonymous") {
        saveBinaryDataDirectory = "Data/";
    } else if(memoryMode == "filesystem") {
        saveBinaryDataDirectory = "DataOnDisk/";
    } else {
        SHASTA_ASSERT(0);
    }

    std::filesystem::create_directory(saveBinaryDataDirectory);

    // initiateSaveBinaryData(&Assembler::saveMarkers);
    // SHASTA_ASSERT(0);
}



void Assembler::initiateSaveBinaryData(SaveBinaryDataFunction function)
{
    if(not saveBinaryDataDirectory.empty()) {
        saveBinaryDataThreads.push_back(std::thread(function, this));
    }
}

void Assembler::waitForSaveBinaryDataThreads()
{
    for(std::thread& t: saveBinaryDataThreads) {
        t.join();
    }
}

void Assembler::saveMarkers() const
{
    markers.save(saveBinaryDataDirectory + "Markers");
}




