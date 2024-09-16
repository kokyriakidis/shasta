#include "ShortBaseSequence.hpp"
#include "algorithm.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;

#include <iomanip>



void shasta::testShortBaseSequence()
{
    ShortBaseSequence16 s;
    s.set(0, Base::fromCharacter('T'));
    s.set(1, Base::fromCharacter('C'));
    s.set(2, Base::fromCharacter('G'));
    s.set(3, Base::fromCharacter('T'));
    cout << s << " " << s.reverseComplement(6) << endl;
    s.shiftLeft();
    cout << s << endl;

    // const auto oldFill = cout.fill('0');
    for(const uint16_t x: s.data) {
        cout << std::setw(2) << std::hex << int(x) << endl;
        // cout << int(x) << endl;
    }
    // cout.fill(oldFill);

    // Check that constructor from id does the inverse of function id().
    const ShortBaseSequence16 t(s.id(4), 4);
    SHASTA_ASSERT(t == s);


    // Verify that the KmerId for a k-mer of given length k is the
    // same regardless of how the k-mer is stored.
    {
        const string sequenceString = "TCGAGCTTAG";
        const uint64_t k = sequenceString.size();

        ShortBaseSequence16 s16;
        ShortBaseSequence32 s32;
        ShortBaseSequence64 s64;
        for(uint64_t i=0; i<k; i++) {
            const Base base = Base::fromCharacter(sequenceString[i]);
            s16.set(i, base);
            s32.set(i, base);
            s64.set(i, base);
        }
        const uint32_t kmerId16 = s16.id(k);
        const uint64_t kmerId32 = s32.id(k);
        const __uint128_t kmerId64 = s64.id(k);

        cout << kmerId16 << " " << kmerId32 << " " << kmerId64 << endl;
        SHASTA_ASSERT(kmerId16 == kmerId32);
        SHASTA_ASSERT(kmerId32 == kmerId64);

    }
}
