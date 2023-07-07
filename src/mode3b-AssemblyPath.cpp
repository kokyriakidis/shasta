#include "mode3b-AssemblyPath.hpp"
#include "mode3b-PathFiller1.hpp"
#include "mode3b-PathFinder.hpp"
#include "Assembler.hpp"
#include "MarkerInterval.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;
using namespace mode3b;

#include <iostream.hpp>



AssemblyPath::AssemblyPath(
    const Assembler& assembler,
    MarkerGraphEdgeId startEdgeId,
    uint64_t direction  // 0 = forward, 1 = backward, 2=bidirectional
    ) :
    assembler(assembler)
{
    create(startEdgeId, direction);
    assemble();
}



// Create the primaryEdges and the steps.
void AssemblyPath::create(
    MarkerGraphEdgeId startEdgeId,
    uint64_t direction  // 0 = forward, 1 = backward, 2=bidirectional
    )
{

    // Forward.
    if(direction == 0) {

        vector< pair<MarkerGraphEdgeId, MarkerGraphEdgePairInfo> > otherPrimaryEdges;
        PathFinder pathFinder(assembler, startEdgeId, direction, otherPrimaryEdges);

        // Create the primaryEdges and the steps.
        primaryEdges.push_back(startEdgeId);
        for(const auto& p: otherPrimaryEdges) {
            primaryEdges.push_back(p.first);
            steps.push_back(Step(p.second));
        }
    }



    // Backward.
    else if(direction == 1) {

        vector< pair<MarkerGraphEdgeId, MarkerGraphEdgePairInfo> > otherPrimaryEdges;
        PathFinder pathFinder(assembler, startEdgeId, direction, otherPrimaryEdges);

        // Reverse the other primary edges.
        reverse(otherPrimaryEdges.begin(), otherPrimaryEdges.end());
        for(auto& p: otherPrimaryEdges) {
            p.second.reverse();
        }

        // Create the primaryEdges and the steps.
        for(const auto& p: otherPrimaryEdges) {
            primaryEdges.push_back(p.first);
            steps.push_back(Step(p.second));
        }
        primaryEdges.push_back(startEdgeId);
    }



    // Bidirectional
    else if(direction == 2) {

        // Move forward.
        vector< pair<MarkerGraphEdgeId, MarkerGraphEdgePairInfo> > otherPrimaryEdgesForward;
        PathFinder forwardPathFinder(assembler, startEdgeId, 0, otherPrimaryEdgesForward);

        // Move backward.
        vector< pair<MarkerGraphEdgeId, MarkerGraphEdgePairInfo> > otherPrimaryEdgesBackward;
        PathFinder backWardPathFinder(assembler, startEdgeId, 1, otherPrimaryEdgesBackward);
        reverse(otherPrimaryEdgesBackward.begin(), otherPrimaryEdgesBackward.end());
        for(auto& p: otherPrimaryEdgesBackward) {
            p.second.reverse();
        }

        // Combine them.
        for(const auto& p: otherPrimaryEdgesBackward) {
            primaryEdges.push_back(p.first);
            steps.push_back(Step(p.second));
        }
        primaryEdges.push_back(startEdgeId);
        for(const auto& p: otherPrimaryEdgesForward) {
            primaryEdges.push_back(p.first);
            steps.push_back(Step(p.second));
        }
    }



    SHASTA_ASSERT(primaryEdges.size() == steps.size() + 1);
}



// Assemble the sequence of each Step.
void AssemblyPath::assemble()
{
    for(uint64_t i=0; i<steps.size(); i++) {
        const MarkerGraphEdgeId edgeIdA = primaryEdges[i];
        const MarkerGraphEdgeId edgeIdB = primaryEdges[i+1];
        cout << "Assembling between primary edges " <<
            edgeIdA << " " << edgeIdB << endl;
        Step& step = steps[i];
        ostream html(0);
        PathFiller1 pathFiller(assembler, edgeIdA, edgeIdB, html, false, false, false, false, false);
        pathFiller.getSecondarySequence(step.sequence);
        cout << "Coverage " << pathFiller.coverage() << " , assembled length " << step.sequence.size() << endl;
    }
}



AssemblyPath::Step::Step(const MarkerGraphEdgePairInfo& info) :
    info(info)
{}



void AssemblyPath::getSequence(vector<Base>& sequence) const
{
    sequence.clear();
    for(uint64_t i=0; /* Check later */ ; i++) {

        // Append the primary edge sequence.
        const MarkerGraphEdgeId edgeId = primaryEdges[i];
        const auto edgeSequence = assembler.markerGraph.edgeSequence[edgeId];
        copy(edgeSequence.begin(), edgeSequence.end(), back_inserter(sequence));

        // If this is the last primary edge, we are done.
        if(i == primaryEdges.size() - 1) {
            break;
        }

        // Append the step sequence.
        const Step& step = steps[i];
        copy(step.sequence.begin(), step.sequence.end(), back_inserter(sequence));
    }
}


void AssemblyPath::writeFasta(ostream& fasta) const
{
    vector<Base> sequence;
    getSequence(sequence);

    fasta << ">Path " << sequence.size() << "\n";
    copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(fasta));
    fasta << "\n";
}



void AssemblyPath::writeCsv(ostream& csv) const
{
    csv << "Step,EdgeId,Begin,End,Length,Sequence\n";

    uint64_t positionBegin = 0;

    for(uint64_t i=0; /* Check later */ ; i++) {

        // Write a line for this primary edge.
        {
            const MarkerGraphEdgeId edgeId = primaryEdges[i];
            const auto edgeSequence = assembler.markerGraph.edgeSequence[edgeId];
            const uint64_t length = edgeSequence.size();
            const uint64_t positionEnd = positionBegin + length;

            csv << i << ",";
            csv << edgeId << ",";
            csv << positionBegin << ",";
            csv << positionEnd << ",";
            csv << length << ",";
            copy(edgeSequence.begin(), edgeSequence.end(), ostream_iterator<Base>(csv));
            csv << "\n";

            positionBegin = positionEnd;
        }



        // If this is the last primary edge, we are done.
        if(i == primaryEdges.size() - 1) {
            break;
        }

        // Write a line with the sequence of the step between this primary edge and the next.
        {
            const Step& step = steps[i];
            const uint64_t length = step.sequence.size();
            const uint64_t positionEnd = positionBegin + length;

            csv << i << ",";
            csv << ",";
            csv << positionBegin << ",";
            csv << positionEnd << ",";
            csv << length << ",";

            // The sequence can be long, so for convenience we write it on multiple lines.
            // To keep all lines in the same cell, we must quote the entire sequence.
            csv << "\"";
            const uint64_t basesPerLine = 100;
            for(uint64_t i=0; i<step.sequence.size(); i++) {
                if(i != 0 and ((i % basesPerLine) == 0)) {
                    csv << "\n";
                }
                csv << step.sequence[i];
            }
            csv << "\"\n";

            positionBegin = positionEnd;
        }
    }

}
