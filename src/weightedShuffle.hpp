#pragma once

// Weighted shuffle using the algorithm by Efraimidis and Spirakis, 2005
// https://utopia.duth.gr/~pefraimi/research/data/2007EncOfAlg.pdf
// After shuffling, elements with a large weight are more likely
// to be placed at the beginning of the vector.

#include <algorithm.hpp>
#include <cmath>
#include <random>
#include <vector.hpp>

namespace shasta {
    template<class T> class WeightedShuffleItem;
    template<class T, class RandomSource> void weightedShuffle(
        vector<WeightedShuffleItem<T> >&,
        RandomSource&);
}



template<class T> class shasta::WeightedShuffleItem {
public:

    // One of the item to be shuffled.
    T t;

    // Its weight, given on input to weightedShuffle.
    double w;

    // Space for the sorting key in by weightedShuffle.
    // This does not need to be initialized before calling weightedShuffle.
    double k;

    WeightedShuffleItem(const T& t, double w) :
        t(t), w(w) {}

    // Sort with large weights closer to the beginning of the vector.
    bool operator<(const WeightedShuffleItem& that) const
    {
        return k > that.k;
    }
};



template<class T, class RandomSource> void shasta::weightedShuffle(
    vector< WeightedShuffleItem<T> >& v,
    RandomSource& randomSource)
{
    std::uniform_real_distribution<double> distribution;

    // Fill in the keys. See https://utopia.duth.gr/~pefraimi/research/data/2007EncOfAlg.pdf
    for(WeightedShuffleItem<T>& x: v) {
        x.k = std::pow(distribution(randomSource), 1. / x.w);
    }

    // Sort by the keys computed as above.
    // Elemenets with large weight are more likely to end up near the end of the vector.
    sort(v.begin(), v.end());
}
