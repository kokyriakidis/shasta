#pragma once

#include "deduplicate.hpp"
#include "ShortBaseSequence.hpp"

#include "vector.hpp"

namespace shasta
{
    template<class Int> void applySingleEdit(
        const ShortBaseSequence<Int>&,
        uint64_t k,
        vector< ShortBaseSequence<Int> >&
        );
}



template<class Int> void shasta::applySingleEdit(
    const ShortBaseSequence<Int>& x,
    uint64_t k,
    vector< ShortBaseSequence<Int> >& v
    )
{
    v.clear();

    // Apply one base mismatches.
    for(uint64_t i=0; i<k; i++) {
        const uint64_t ibx = x[i].value;
        for(uint64_t iby=0; iby<4; iby++) {
            if(iby != ibx) {
                ShortBaseSequence<Int> y = x;
                y.set(i, Base::fromInteger(iby));
                v.push_back(y);
            }
        }
    }


    // Apply 1-base deletions.
    // We need to keep the length the same, so this means that
    // we have to also insert a base at the beginning or the end.
    for(uint64_t deletionPosition=0; deletionPosition<k; deletionPosition++) {
        ShortBaseSequence<Int> y;
        for(uint64_t j=0; j<deletionPosition; j++) {
            y.set(j, x[j]);
        }
        for(uint64_t j=deletionPosition+1; j<k; j++) {
            y.set(j-1, x[j]);
        }

        // Add one base at the end.
        ShortBaseSequence<Int> z = y;
        for(uint64_t ibz=0; ibz<4; ibz++) {
            z.set(k-1, Base::fromInteger(ibz));
            if(z != x) {
                v.push_back(z);
            }
        }

        // Add one base at the beginning.
        y.shiftRight();
        z = y;
        for(uint64_t ibz=0; ibz<4; ibz++) {
            z.set(0, Base::fromInteger(ibz));
            if(z != x) {
                v.push_back(z);
            }
        }
    }


    // Apply 1-base insertions.
    // We need to keep the length the same, so this means that
    // we have to also delete a base at the beginning or the end.
    for(uint64_t insertionPosition=0; insertionPosition<k; insertionPosition++) {
        for(uint64_t iby=0; iby<4; iby++) {
            ShortBaseSequence<Int> y;
            for(uint64_t j=0; j<insertionPosition; j++) {
                y.set(j, x[j]);
            }
            y.set(insertionPosition, Base::fromInteger(iby));
            for(uint64_t j=insertionPosition; j<k; j++) {
                y.set(j+1, x[j]);
            }

            // Remove the last base.
            ShortBaseSequence<Int> y1 = y;
            y1.set(k, Base::fromInteger(Int(0)));
            if(y1 != x) {
                v.push_back(y1);
            }

            // Remove the first base.
            ShortBaseSequence<Int> y2 = y;
            y2.shiftLeft();
            if(y2 != x) {
                v.push_back(y2);
            }
        }
    }

    deduplicate(v);
}
