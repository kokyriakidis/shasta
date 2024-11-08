/********************************************************************************

Creation of Anchors from json input.

********************************************************************************/

#include "mode3-Anchor.hpp"
#include "Base.hpp"
using namespace shasta;
using namespace mode3;



Anchors::Anchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const string& /* jsonFileName */,
    uint64_t /* minPrimaryCoverage */,
    uint64_t /* maxPrimaryCoverage */,
    uint64_t /* threadCount */) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    k(k),
    markers(markers)
{
    throw runtime_error("Anchor creation from json is not yet implemented.");
}
