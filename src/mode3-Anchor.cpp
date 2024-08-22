#include "mode3-Anchor.hpp"
#include "MarkerGraph.hpp"
using namespace shasta;
using namespace mode3;


Anchors::Anchors(const MarkerGraph& markerGraph) :
    anchorMarkerIntervals(markerGraph.edgeMarkerIntervals)
{
}


Anchor Anchors::operator[](AnchorId anchorId) const
{
    return anchorMarkerIntervals[anchorId];
}



uint64_t Anchors::size() const
{
    return anchorMarkerIntervals.size();
}
