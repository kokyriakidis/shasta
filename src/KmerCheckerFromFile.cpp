#include "KmerCheckerFromFile.hpp"
#include "Kmer.hpp"
using namespace shasta;

#include "algorithm.hpp"
#include "fstream.hpp"



// Initial creation.
KmerCheckerFromFile::KmerCheckerFromFile(
    uint64_t k,
    const string& fileName,
    const MappedMemoryOwner& mappedMemoryOwner) :
    MappedMemoryOwner(mappedMemoryOwner),
    k(k)
{
    // Check that fileName specifies an absolute path.
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

    // Initialize the marker k-mer table.
    markerKmers.createNew(largeDataName("KmerCheckerFromFile"), largeDataPageSize);

    // Read one k-mer per line.
    uint64_t lineCount = 1;
    string line;
    while(true) {

        // Read a line.
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

        // Compute the KmerId of this Kmer and its reverse complement.
        const Kmer kmerRc = kmer.reverseComplement(k);
        const KmerId kmerId = KmerId(kmer.id(k));
        const KmerId kmerIdRc = KmerId(kmerRc.id(k));

        // Store the lowest of the two in the table.
        markerKmers.push_back(min(kmerId, kmerIdRc));

        ++lineCount;

    }

    // Sort and deduplicate the marker k-mers.
    sort(markerKmers.begin(), markerKmers.end());
    auto it = std::unique(markerKmers.begin(), markerKmers.end());
    markerKmers.resize(it - markerKmers.begin());
    markerKmers.unreserve();
}



// Creation from binary data.
KmerCheckerFromFile::KmerCheckerFromFile(
    uint64_t k,
    const MappedMemoryOwner& mappedMemoryOwner) :
    MappedMemoryOwner(mappedMemoryOwner),
    k(k)
{
    markerKmers.accessExistingReadOnly(largeDataName("KmerCheckerFromFile"));
}



bool KmerCheckerFromFile::isMarker(KmerId kmerId) const
{
    // Find the reverse complement KmerId.
    const Kmer kmer(kmerId, k);
    const Kmer kmerRc = kmer.reverseComplement(k);
    const KmerId kmerIdRc = KmerId(kmerRc.id(k));

    // Look for the lowest of the two in our sorted table of marker kmers.
    return std::binary_search(markerKmers.begin(), markerKmers.end(), min(kmerId, kmerIdRc));
}
