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

        // Test maxHomopolymerLength.
        cout << "maxHomopolymerLength returned " <<
            kmer.maxHomopolymerLength(k) << endl;

        // Test countExactRepeatCopies.
        cout << "countExactRepeatCopies<1> returned " <<
            kmer.countExactRepeatCopies<1>(k) << endl;
        cout << "countExactRepeatCopies<2> returned " <<
            kmer.countExactRepeatCopies<2>(k) << endl;
        cout << "countExactRepeatCopies<3> returned " <<
            kmer.countExactRepeatCopies<3>(k) << endl;
        cout << "countExactRepeatCopies<4> returned " <<
            kmer.countExactRepeatCopies<4>(k) << endl;
        cout << "countExactRepeatCopies<5> returned " <<
            kmer.countExactRepeatCopies<5>(k) << endl;
        cout << "countExactRepeatCopies<6> returned " <<
            kmer.countExactRepeatCopies<6>(k) << endl;

    }
}
