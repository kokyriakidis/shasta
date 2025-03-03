#include "bitReversal.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;

#include <algorithm.hpp>
#include <array.hpp>
#include <iostream.hpp>
#include <random>

void shasta::testBitReversal()
{
    array<int, 128> bits;
    array<int, 128> reversedBits;

    std::uniform_int_distribution distribution(0, 1);
    std::mt19937 generator;

    for(uint64_t iteration=0; iteration<10; iteration++) {

        // Get random bits.
        for(int& bit: bits) {
            bit = distribution(generator);
        }

        // Create a __uint128_t with these bits.
        __uint128_t x = 0;
        for(uint64_t i=0; i<bits.size(); i++) {
            if(bits[i]) {
                x |= __uint128_t(1ULL) << i;
            }
        }

        // Bit reversal.
        const __uint128_t y = bitReversal(x);

        // Extract the reversed bits.
        for(uint64_t i=0; i<bits.size(); i++) {
            reversedBits[i] = (y >> i) & 1ULL;
        }

        // Reverse them again.
        std::reverse(reversedBits.begin(), reversedBits.end());

        for(const int bit: bits) {
            cout << int(bit);
        }
        cout << endl;

        for(const int bit: reversedBits) {
            cout << int(bit);
        }
        cout << endl;

        // They should be the same as what we started with.
        SHASTA_ASSERT(reversedBits == bits);

    }
}
