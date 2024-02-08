#include "mode3b-AssemblyPath.hpp"
// #include "mode3b-PathFiller1.hpp"
// #include "mode3b-PathFiller2.hpp"
#include "mode3b-PathFiller3.hpp"
#include "mode3b-PathFinder.hpp"
#include "Assembler.hpp"
#include "MarkerInterval.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;
using namespace mode3b;

#include <iostream.hpp>

#include "MultithreadedObject.tpp"
template class MultithreadedObject<AssemblyPath>;

// Create the assembly path starting from a given primary edge.
AssemblyPath::AssemblyPath(
    const Assembler& assembler,
    MarkerGraphEdgeId startEdgeId,
    uint64_t direction  // 0 = forward, 1 = backward, 2=bidirectional
    ) :
    MultithreadedObject<AssemblyPath>(*this),
    assembler(assembler)
{
    create(startEdgeId, direction);
    assembleSequential();
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



// Create the assembly path given n primary edges and
// the n-1 MarkerGraphEdgePairInfo between consecutive primary edges.
AssemblyPath::AssemblyPath(
    const Assembler& assembler,
    const vector<MarkerGraphEdgeId>& primaryEdges,
    const vector<MarkerGraphEdgePairInfo>& infos,
    uint64_t threadCount) :
    MultithreadedObject<AssemblyPath>(*this),
    assembler(assembler),
    primaryEdges(primaryEdges)
{
    // Sanity check.
    SHASTA_ASSERT(not infos.empty());
    SHASTA_ASSERT(infos.size() == primaryEdges.size() - 1);

    // Fill in the steps. The primaryEdges were already filled in above.
    for(const MarkerGraphEdgePairInfo& info: infos) {
        steps.push_back(Step(info));
    }

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Assemble in parallel.
    assembleParallel(threadCount);
}



// Assemble the sequence of each Step.
void AssemblyPath::assembleSequential()
{
    for(uint64_t i=0; i<steps.size(); i++) {
        assembleStep(i);
    }
}



void AssemblyPath::assembleStep(uint64_t i)
{
    const MarkerGraphEdgeId edgeIdA = primaryEdges[i];
    const MarkerGraphEdgeId edgeIdB = primaryEdges[i+1];
    // cout << "Assembling between primary edges " << edgeIdA << " " << edgeIdB << endl;
    Step& step = steps[i];
    ostream html(0);

#if 0
    PathFiller1 pathFiller(assembler, edgeIdA, edgeIdB, html, false, false, false, false, false);
    pathFiller.getSecondarySequence(step.sequence);
#endif
#if 0
    PathFiller2 pathFiller(assembler, edgeIdA, edgeIdB, html);
    pathFiller.getSecondarySequence(step.sequence);
#endif

    try {
        PathFiller3 pathFiller(assembler, edgeIdA, edgeIdB, 0, html);
        pathFiller.getSecondarySequence(step.sequence);
    } catch (...) {
        cout << "Error occurred when assembling between marker graph edges " <<
            edgeIdA << " and " << edgeIdB << endl;
        throw;
    }

}



void AssemblyPath::assembleParallel(uint64_t threadCount)
{
    setupLoadBalancing(steps.size(), 1);
    runThreads(&AssemblyPath::assembleThreadFunction, threadCount);
}



void AssemblyPath::assembleThreadFunction(uint64_t threadId)
{
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(uint64_t i=begin; i<end; ++i) {
            assembleStep(i);
        }
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



void AssemblyPath::getInternalSequence(vector<Base>& sequence) const  // Excluding the begin/end edges.
{
    sequence.clear();
    for(uint64_t i=0; /* Check later */ ; i++) {

        // Append the primary edge sequence.
        if(i>0 and i<primaryEdges.size()-1) {
            const MarkerGraphEdgeId edgeId = primaryEdges[i];
            const auto edgeSequence = assembler.markerGraph.edgeSequence[edgeId];
            copy(edgeSequence.begin(), edgeSequence.end(), back_inserter(sequence));
        }

        // If this is the last primary edge, we are done.
        if(i == primaryEdges.size() - 1) {
            break;
        }

        // Append the step sequence.
        const Step& step = steps[i];
        copy(step.sequence.begin(), step.sequence.end(), back_inserter(sequence));
    }

}



void AssemblyPath::writeFasta(ostream& fasta, const string& name) const
{
    vector<Base> sequence;
    getSequence(sequence);

    fasta << ">" << name << " " << sequence.size() << "\n";
    copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(fasta));
    fasta << "\n";
}



void AssemblyPath::writeCsv(ostream& csv, const string& name) const
{

    uint64_t positionBegin = 0;

    for(uint64_t i=0; /* Check later */ ; i++) {

        // Write a line for this primary edge.
        {
            const MarkerGraphEdgeId edgeId = primaryEdges[i];
            const auto edgeSequence = assembler.markerGraph.edgeSequence[edgeId];
            const uint64_t length = edgeSequence.size();
            const uint64_t positionEnd = positionBegin + length;

            csv << name << ",";
            csv << i << ",";
            csv << edgeId << ",";
            csv << positionBegin << ",";
            csv << positionEnd << ",";
            csv << length << ",";
            // copy(edgeSequence.begin(), edgeSequence.end(), ostream_iterator<Base>(csv));
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

            csv << name << ",";
            csv << i << ",";
            csv << ",";
            csv << positionBegin << ",";
            csv << positionEnd << ",";
            csv << length << ",";

            // The sequence can be long, so for convenience we write it on multiple lines.
            // To keep all lines in the same cell, we must quote the entire sequence.
#if 0
            csv << "\"";
            const uint64_t basesPerLine = 100;
            for(uint64_t i=0; i<step.sequence.size(); i++) {
                if(i != 0 and ((i % basesPerLine) == 0)) {
                    csv << "\n";
                }
                csv << step.sequence[i];
            }
            csv << "\""
#endif
            csv << "\n";

            positionBegin = positionEnd;
        }
    }

}
