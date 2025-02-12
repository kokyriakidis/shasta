#include "ShortBaseSequence.hpp"
#include "ShortBaseSequenceEdit.hpp"
#include "algorithm.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;

#include <iomanip>



void shasta::testShortBaseSequence()
{
    while(true) {

        // Read the k-mer string.
        string kmerString;
        cout << "Enter a k-mer." << endl;
        cin >> kmerString;
        const uint64_t k = kmerString.size();
        if(k > 64) {
            cout << "Can be at most 64 bases long." << endl;
        }

        // Create the k-mer.
        ShortBaseSequence64 kmer;
        for(uint64_t i=0; i<k; i++) {
            kmer.set(i, Base::fromCharacter(kmerString[i]));
        }

        // Write out the maximum homopolymer length.
        cout << "Maximum homopolymer length is " <<
            kmer.maxHomopolymerLength(k) << endl;
    }
}
