// Shasta.
#include "mode3-LocalAssemblyGraph.hpp"
#include "computeLayout.hpp"
#include "html.hpp"
#include "HttpServer.hpp"
#include "MurmurHash2.hpp"
#include "platformDependent.hpp"
#include "runCommandWithTimeout.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include <queue>



LocalAssemblyGraph::LocalAssemblyGraph(
    const AssemblyGraphPostprocessor& assemblyGraph,
    const vector<ChainIdentifier>& startingChains,
    uint64_t maxDistance,
    const string& assemblyStage) :
    assemblyGraph(assemblyGraph),
    maxDistance(maxDistance),
    assemblyStage(assemblyStage)
{
    addVertices(startingChains, maxDistance);
    addEdges();
}



void LocalAssemblyGraph::addVertices(
    const vector<ChainIdentifier>& startingChains,
    uint64_t maxDistance)
{
    LocalAssemblyGraph& localAssemblyGraph = *this;

    // We will use a BFS to create the vertices.
    std::queue<vertex_descriptor> q;

    // Initialize the queue by adding vertices corresponding to the initial chains.
    /// These vertices have distance 0.
    for(const ChainIdentifier& chainIdentifier: startingChains) {

        /*
        cout << "Creating starting vertices for Chain " <<
            assemblyGraph.chainStringId(chainIdentifier.e, chainIdentifier.positionInBubbleChain, chainIdentifier.indexInBubble)
            << endl;
        */

        // Add a LocalAssemblyGraphVertex for the source of this Chain.
        {
            const LocalAssemblyGraphVertex vertex = createVertexAtChainSource(chainIdentifier);
            if(not vertexMap.contains(vertex)) {
                const vertex_descriptor v = add_vertex(vertex, localAssemblyGraph);
                vertexMap.insert({vertex, v});
                q.push(v);
                /*
                cout << "Added starting vertex: ";
                writeVertex(vertex, cout);
                cout << endl;
                */
            }
        }

        // Add a LocalAssemblyGraphVertex for the target of this Chain.
        {
            const LocalAssemblyGraphVertex vertex = createVertexAtChainTarget(chainIdentifier);
            if(not vertexMap.contains(vertex)) {
                const vertex_descriptor v = add_vertex(vertex, localAssemblyGraph);
                vertexMap.insert({vertex, v});
                q.push(v);
                /*
                cout << "Added starting vertex: ";
                writeVertex(vertex, cout);
                cout << endl;
                */
            }
        }
    }



    // BFS loop.
    vector<LocalAssemblyGraphVertex> neighbors;
    while(not q.empty()) {

        const vertex_descriptor v0 = q.front();
        q.pop();
        const LocalAssemblyGraphVertex& vertex0 = localAssemblyGraph[v0];

        /*
        cout << "Dequeued ";
        writeVertex(vertex0, cout);
        cout << endl;
        */

        neighbors.clear();
        getNeighbors(vertex0, neighbors);

        for(const LocalAssemblyGraphVertex& vertex1: neighbors) {
            if(vertexMap.contains(vertex1)) {
                continue;
            }

            /*
            cout << "Found ";
            writeVertex(vertex1, cout);
            cout << endl;
            */

            const vertex_descriptor v1 = add_vertex(vertex1, localAssemblyGraph);
            LocalAssemblyGraphVertex& addedVertex1 = localAssemblyGraph[v1];
            addedVertex1.distance = vertex0.distance + 1;
            vertexMap.insert({vertex1, v1});
            if(addedVertex1.distance < maxDistance) {
                q.push(v1);
            }
        }
    }


    /*
    cout << "Vertices: " << endl;
    BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph) {
        const LocalAssemblyGraphVertex& vertex = localAssemblyGraph[v];
        writeVertex(vertex, cout);
        cout << " distance " << vertex.distance << endl;
    }
    */
}



void LocalAssemblyGraph::addEdges()
{
    LocalAssemblyGraph& localAssemblyGraph = *this;

    BGL_FORALL_VERTICES(v0, localAssemblyGraph, LocalAssemblyGraph) {
        const LocalAssemblyGraphVertex& vertex0 = localAssemblyGraph[v0];

        if(vertex0.isTypeA()) {

            // Vertex0 corresponds to a vertex of the assembly graph.
            // Loop over all outgoing chains.
            const AssemblyGraph::vertex_descriptor av0 = vertex0.v;
            BGL_FORALL_OUTEDGES(av0, ae, assemblyGraph, AssemblyGraph) {
                const AssemblyGraphEdge& aEdge = assemblyGraph[ae];
                const BubbleChain& bubbleChain = aEdge;
                const Bubble& firstBubble = bubbleChain.front();
                for(uint64_t indexInBubble=0; indexInBubble<firstBubble.size(); indexInBubble++) {
                    const LocalAssemblyGraphEdge edge(ChainIdentifier(ae, 0, indexInBubble));
                    if(bubbleChain.size() == 1) {
                        const LocalAssemblyGraphVertex& vertex1(target(ae, assemblyGraph));
                        auto it = vertexMap.find(vertex1);
                        if(it != vertexMap.end()) {
                            const vertex_descriptor v1 = it->second;
                            add_edge(v0, v1, edge, localAssemblyGraph);
                        }
                    } else {
                        const LocalAssemblyGraphVertex vertex1(ae, 0);
                        auto it = vertexMap.find(vertex1);
                        if(it != vertexMap.end()) {
                            const vertex_descriptor v1 = it->second;
                            add_edge(v0, v1, edge, localAssemblyGraph);
                        }

                    }
                }
            }

        } else {

            /*
            cout << "Adding edges with source vertex ";
            writeVertex(vertex0, cout);
            cout << endl;
            */

            // Vertex0 does not correspond to a vertex of the assembly graph.
            // Loop over all outgoing chains.
            const AssemblyGraph::edge_descriptor ae = vertex0.e;
            const AssemblyGraphEdge& aEdge = assemblyGraph[ae];
            const BubbleChain& bubbleChain = aEdge;
            const Bubble& nextBubble = bubbleChain[vertex0.positionInBubbleChain + 1];
            for(uint64_t indexInBubble=0; indexInBubble<nextBubble.size(); indexInBubble++) {
                const LocalAssemblyGraphEdge edge(ChainIdentifier(ae, vertex0.positionInBubbleChain + 1, indexInBubble));
                if(vertex0.positionInBubbleChain == bubbleChain.size() - 2) {
                    // cout << "Case 1" << endl;
                    const LocalAssemblyGraphVertex& vertex1(target(ae, assemblyGraph));
                    auto it = vertexMap.find(vertex1);
                    if(it != vertexMap.end()) {
                        // cout << "Case 1 added" << endl;
                        const vertex_descriptor v1 = it->second;
                        add_edge(v0, v1, edge, localAssemblyGraph);
                    }
                } else {
                    // cout << "Case 2" << endl;
                    const LocalAssemblyGraphVertex vertex1(ae, vertex0.positionInBubbleChain + 1);
                    auto it = vertexMap.find(vertex1);
                    if(it != vertexMap.end()) {
                        // cout << "Case 2 added." << endl;
                        const vertex_descriptor v1 = it->second;
                        add_edge(v0, v1, edge, localAssemblyGraph);
                    }
                }
            }
        }


    }


    /*
    cout << " Edges:" << endl;
    BGL_FORALL_EDGES(e, localAssemblyGraph, LocalAssemblyGraph) {
        writeEdge(e, cout);
        cout << endl;
    }
    */
}


LocalAssemblyGraphVertex LocalAssemblyGraph::createVertexAtChainSource(
    const ChainIdentifier& chainIdentifier) const
{
    if(chainIdentifier.positionInBubbleChain == 0) {

        // Type A vertex.
        const AssemblyGraph::vertex_descriptor av = source(chainIdentifier.e, assemblyGraph);
        return LocalAssemblyGraphVertex(av);

    } else {

        // Type B vertex.
        return LocalAssemblyGraphVertex(chainIdentifier.e, chainIdentifier.positionInBubbleChain - 1);

    }
}



LocalAssemblyGraphVertex LocalAssemblyGraph::createVertexAtChainTarget(
    const ChainIdentifier& chainIdentifier) const
{
    const BubbleChain& bubbleChain = assemblyGraph[chainIdentifier.e];

    if(chainIdentifier.positionInBubbleChain == bubbleChain.size() - 1) {

        // Type A vertex.
        const AssemblyGraph::vertex_descriptor av = target(chainIdentifier.e, assemblyGraph);
        return LocalAssemblyGraphVertex(av);

    } else {

        // Type B vertex.
        return LocalAssemblyGraphVertex(chainIdentifier.e, chainIdentifier.positionInBubbleChain);
    }
}



void LocalAssemblyGraph::getNeighbors(
    const LocalAssemblyGraphVertex& vertex,
    vector<LocalAssemblyGraphVertex>& neighbors) const
{
    getChildren(vertex, neighbors);
    getParents(vertex, neighbors);
}



void LocalAssemblyGraph::getChildren(
    const LocalAssemblyGraphVertex& vertex0,
    vector<LocalAssemblyGraphVertex>& children) const
{

    if(vertex0.isTypeA()) {
        BGL_FORALL_OUTEDGES(vertex0.v, ae, assemblyGraph, AssemblyGraph) {
            ChainIdentifier chainIdentifier;
            chainIdentifier.e = ae;
            chainIdentifier.positionInBubbleChain = 0;
            chainIdentifier.indexInBubble = 0;
            children.push_back(createVertexAtChainTarget(chainIdentifier));
        }
    } else {
        ChainIdentifier chainIdentifier;
        chainIdentifier.e = vertex0.e;
        chainIdentifier.positionInBubbleChain = vertex0.positionInBubbleChain + 1;
        chainIdentifier.indexInBubble = 0;
        children.push_back(createVertexAtChainTarget(chainIdentifier));
    }
}



void LocalAssemblyGraph::getParents(
    const LocalAssemblyGraphVertex& vertex0,
    vector<LocalAssemblyGraphVertex>& parents) const
{
    /*
    cout << "LocalAssemblyGraph::getParents called for ";
    writeVertex(vertex0, cout);
    cout << endl;
    */

    if(vertex0.isTypeA()) {
        BGL_FORALL_INEDGES(vertex0.v, ae, assemblyGraph, AssemblyGraph) {
            const BubbleChain& bubbleChain = assemblyGraph[ae];
            ChainIdentifier chainIdentifier;
            chainIdentifier.e = ae;
            chainIdentifier.positionInBubbleChain = bubbleChain.size() - 1;
            chainIdentifier.indexInBubble = 0;
            parents.push_back(createVertexAtChainSource(chainIdentifier));

            /*
            cout << "LocalAssemblyGraph::getParents added ";
            writeVertex(parents.back(), cout);
            cout << endl;
            */
        }
    } else {
        ChainIdentifier chainIdentifier;
        chainIdentifier.e = vertex0.e;
        chainIdentifier.positionInBubbleChain = vertex0.positionInBubbleChain;
        chainIdentifier.indexInBubble = 0;
        parents.push_back(createVertexAtChainSource(chainIdentifier));

        /*
        cout << "LocalAssemblyGraph::getParents added ";
        writeVertex(parents.back(), cout);
        cout << endl;
        */
    }
    // cout << "LocalAssemblyGraph::getParents ends" << endl;;
}



void LocalAssemblyGraph::writeVertex(
    const LocalAssemblyGraphVertex& vertex,
    ostream& s) const
{
    if(vertex.isTypeA()) {
        const AssemblyGraphVertex aVertex = assemblyGraph[vertex.v];
        s << "Type A vertex with AnchorId " << anchorIdToString(aVertex.anchorId);
    } else {
        const AssemblyGraphEdge& aEdge = assemblyGraph[vertex.e];
        s << "Type B vertex with BubbleChain " << aEdge.id <<
            " position in BubbleChain " << vertex.positionInBubbleChain;
    }
}


void LocalAssemblyGraph::writeEdge(edge_descriptor e, ostream& s) const
{
    const LocalAssemblyGraph& localAssemblyGraph = *this;
    const LocalAssemblyGraphEdge& edge = localAssemblyGraph[e];
    s << assemblyGraph.getChainStringId(edge);
}



void LocalAssemblyGraph::writeGraphviz(
    const string& fileName,
    const LocalAssemblyGraphDisplayOptions& options) const
{
    ofstream file(fileName);
    writeGraphviz(file, options);
}



void LocalAssemblyGraph::writeGraphviz(
    ostream& s,
    const LocalAssemblyGraphDisplayOptions& options ) const
{
    const LocalAssemblyGraph& localAssemblyGraph = *this;

    s << "digraph LocalAssemblyGraph {\n";

    // Write the vertices.
    BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph) {
        const LocalAssemblyGraphVertex& vertex = localAssemblyGraph[v];
        const AnchorId anchorId = getAnchorId(v);
        const string anchorIdString = anchorIdToString(anchorId);
        const uint64_t coverage = assemblyGraph.anchors[anchorId].coverage();

        // For the Graphviz vertex name we use the vertex descriptor (an integer).
        // It has no particular meaning.
        // We cannot use the AnchorId because there is no guarantee that an AnchorId
        // appears only once in the AssemblyGraph.
        s << v;

        // Begin vertex attributes.
        s << " [";

        // Label.
        if(options.vertexLabels) {
            s << "label=\"" << anchorIdString << "\\n" << coverage << "\"";
        }

        // URL
        s << " URL=\"exploreAnchor?anchorIdString=" << HttpServer::urlEncode(anchorIdString) << "\"";

        // Tooltip.
        s << " tooltip=\"" << anchorIdString << " " << coverage << ", distance " << vertex.distance << "\"";

        // Color.
        if(vertex.distance == 0) {
            s << " color=blue";
        } else if(vertex.distance == maxDistance) {
            s << " color=cyan";
        }

        // Size.
        if(not options.vertexLabels) {
            const double displaySize = 20. * (
                (options.vertexSizeByCoverage ?
                options.vertexSize * sqrt(0.1 * double(coverage)) :
                options.vertexSize) / 72.);
            s << " width=" << displaySize ;
            s << " penwidth=" << 0.5 * displaySize;
        }

        // End vertex attributes.
        s << "]";

        // End the line for this vertex.
        s << ";\n";
    }


    BGL_FORALL_EDGES(e, localAssemblyGraph, LocalAssemblyGraph) {
        const LocalAssemblyGraphEdge& edge = localAssemblyGraph[e];
        const vertex_descriptor v0 = source(e, localAssemblyGraph);
        const vertex_descriptor v1 = target(e, localAssemblyGraph);
        const string chainStringId = assemblyGraph.getChainStringId(edge);

        s << v0 << "->" << v1 << " [";

        // Tooltip.
        s << " tooltip=\"" << chainStringId << "\"";

        // Label.
        if(options.edgeLabels) {
            s << " label=\"" << chainStringId << "\"";
        }

        // URL
        s << "URL=\"exploreSegment?assemblyStage=" << assemblyStage << "&segmentName=" <<
            HttpServer::urlEncode(chainStringId) << "\"";

        // Color.
        if(options.edgeColoring == "random") {
            // To decide the color, hash the chainStringId.
            // This way we always get the same color for the same chain.
            const uint32_t hashValue = MurmurHash2(chainStringId.data(), int(chainStringId.size()), 759);
            const double hue = double(hashValue % 360) / 360.;
            s << " color=\"" << std::fixed << std::setprecision(2) << hue << " 0.66 0.75\"";
        }

        // Thickness.
        s << " penwidth=" << 10.*options.edgeThickness;

        // Arrow size.
        s << " arrowsize=" << options.arrowSize;
        s << "];\n";
    }

    s << "}\n";
}



LocalAssemblyGraphDisplayOptions::LocalAssemblyGraphDisplayOptions(const vector<string>& request)
{
    // Figure out if command "customLayout" is available.
    const int commandStatus = std::system("which customLayout > /dev/null");
    SHASTA_ASSERT(WIFEXITED(commandStatus));
    const int returnCode = WEXITSTATUS(commandStatus);
    const bool customLayoutIsAvailable = (returnCode == 0);


    sizePixels = 600;
    HttpServer::getParameterValue(request, "sizePixels", sizePixels);

    layoutMethod = (customLayoutIsAvailable ? "custom" : "sfdp");
    HttpServer::getParameterValue(request, "layoutMethod", layoutMethod);

    vertexSize =  1.;
    HttpServer::getParameterValue(request, "vertexSize", vertexSize);

    string vertexSizeByCoverageString;
    vertexSizeByCoverage = HttpServer::getParameterValue(request,
        "vertexSizeByCoverage", vertexSizeByCoverageString);

    string vertexLabelsString;
    vertexLabels = HttpServer::getParameterValue(request,
        "vertexLabels", vertexLabelsString);

    edgeColoring = "random";
    HttpServer::getParameterValue(request, "edgeColoring", edgeColoring);

    minimumEdgeLength = 1.;
    HttpServer::getParameterValue(request, "minimumEdgeLength", minimumEdgeLength);

    additionalEdgeLengthPerKb = 0.001;
    HttpServer::getParameterValue(request, "additionalEdgeLengthPerKb", additionalEdgeLengthPerKb);

    edgeThickness = 1.;
    HttpServer::getParameterValue(request, "edgeThickness", edgeThickness);

    arrowSize = 1.;
    HttpServer::getParameterValue(request, "arrowSize", arrowSize);

    string edgeLabelsString;
    edgeLabels = HttpServer::getParameterValue(request,
        "edgeLabels", edgeLabelsString);
}



void LocalAssemblyGraphDisplayOptions::writeForm(ostream& html) const
{
    // Figure out if command "customLayout" is available.
    const int commandStatus = std::system("which customLayout > /dev/null");
    SHASTA_ASSERT(WIFEXITED(commandStatus));
    const int returnCode = WEXITSTATUS(commandStatus);
    const bool customLayoutIsAvailable = (returnCode == 0);

    html <<
        "<tr>"
        "<th title='Graphics size in pixels. "
        "Changing this works better than zooming. Make it larger if the graph is too crowded."
        " Ok to make it much larger than screen size.'>"
        "Graphics size in pixels"
        "<td class=centered><input type=text required name=sizePixels size=8 style='text-align:center'" <<
        " value='" << sizePixels <<
        "'>";

    html <<
        "<tr>"
        "<th>Layout method"
        "<td class=left>"
        "<input type=radio required name=layoutMethod value='sfdp'" <<
        (layoutMethod == "sfdp" ? " checked=on" : "") <<
        ">sfdp"
        "<br><input type=radio required name=layoutMethod value='fdp'" <<
        (layoutMethod == "fdp" ? " checked=on" : "") <<
        ">fdp"
        "<br><input type=radio required name=layoutMethod value='neato'" <<
        (layoutMethod == "neato" ? " checked=on" : "") <<
        ">neato"
        "<br><input type=radio required name=layoutMethod value='dot'" <<
        (layoutMethod == "dot" ? " checked=on" : "") <<
        ">dot";

    // If command "customLayout" is available, add an option for that.
    if(customLayoutIsAvailable) {
        html <<
            "<br><input type=radio required name=layoutMethod value='custom'" <<
            (layoutMethod == "custom" ? " checked=on" : "") <<
            ">custom";
    }

    html <<
        "<tr>"
        "<th>Vertices"
        "<td class=left>"
        "<input type=text name=vertexSize style='text-align:center' required size=6 value=" <<
        vertexSize << "> Vertex size (arbitrary units)"
        "<br><input type=checkbox name=vertexSizeByCoverage" <<
        (vertexSizeByCoverage ? " checked" : "") <<
        "> Size proportional to coverage"

        "<br><input type=checkbox name=vertexLabels" <<
        (vertexLabels ? " checked" : "") <<
        "> Labels (dot layout only)";



    html <<
        "<tr>"
        "<th>Edges"
        "<td class=left>"

        "<b>Edge coloring</b>"
        "<br><input type=radio required name=edgeColoring value='black'" <<
        (edgeColoring == "black" ? " checked=on" : "") << ">Black"
        "<br><input type=radio required name=edgeColoring value='random'" <<
        (edgeColoring == "random" ? " checked=on" : "") << ">Random"
        "<hr>"

        "<b>Edge graphics</b>"

        "<br><input type=text name=edgeThickness style='text-align:center' required size=6 value=" <<
        edgeThickness << "> Thickness (arbitrary units)"

        "<br><input type=text name=minimumEdgeLength style='text-align:center' required size=6 value=" <<
        minimumEdgeLength << "> Minimum edge length (arbitrary units)"

        "<br><input type=text name=additionalEdgeLengthPerKb style='text-align:center' required size=6 value=" <<
        additionalEdgeLengthPerKb << "> Additional edge length per Kb (arbitrary units)"

        "<br><input type=text name=arrowSize style='text-align:center' required size=6 value=" <<
        arrowSize << "> Arrow size (arbitrary units)"

        "<hr>"
        "<input type=checkbox name=edgeLabels" <<
        (edgeLabels ? " checked" : "") <<
        "> Labels (dot layout only)";

}



void LocalAssemblyGraph::writeHtml(
    ostream& html,
    const LocalAssemblyGraphDisplayOptions& options,
    const string& assemblyStage)
{
    if((options.layoutMethod == "dot") and (options.vertexLabels or options.edgeLabels)) {

        // Use svg output from graphviz.
        writeHtml1(html, options);

    } else {

        // Compute graph layout and use it to generate svg.
        writeHtml2(html, options, assemblyStage);

    }
}



// This is the code that uses svg output from graphviz.
void LocalAssemblyGraph::writeHtml1(
    ostream& html,
    const LocalAssemblyGraphDisplayOptions& options) const
{


    // Write it out in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    writeGraphviz(dotFileName, options);

    // Use graphviz to compute the layout.
    const string svgFileName = dotFileName + ".svg";
    const string shape = options.vertexLabels ? "rectangle" : "point";
    string command =
        options.layoutMethod +
        " -T svg " + dotFileName + " -o " + svgFileName +
        " -Nshape=" + shape +
        " -Gsize=" + to_string(options.sizePixels/72) + " -Gratio=expand ";
    if(options.vertexLabels) {
        command += " -Goverlap=false";
    }
    // cout << "Running command: " << command << endl;
    const int timeout = 30;
    bool timeoutTriggered = false;
    bool signalOccurred = false;
    int returnCode = 0;
    runCommandWithTimeout(command, timeout, timeoutTriggered, signalOccurred, returnCode);
    if(signalOccurred) {
        html << "Error during graph layout. Command was<br>" << endl;
        html << command;
        return;
    }
    if(timeoutTriggered) {
        html << "Timeout during graph layout." << endl;
        return;
    }
    if(returnCode!=0 ) {
        html << "Error during graph layout. Command was<br>" << endl;
        html << command;
        return;
    }
    std::filesystem::remove(dotFileName);



    // Write the svg to html.
    html << "<p><div style='border:solid;display:inline-block'>";
    ifstream svgFile(svgFileName);
    html << svgFile.rdbuf();
    svgFile.close();
    html << "</div>";

    // Remove the .svg file.
    std::filesystem::remove(svgFileName);

    // Add drag and zoom.
    addSvgDragAndZoom(html);
}



// This is the code that computes the graph layout,
// then creates the svg.
void LocalAssemblyGraph::writeHtml2(
    ostream& html,
    const LocalAssemblyGraphDisplayOptions& options,
    const string& assemblyStage)
{
    // Use scientific notation because svg does not accept floating points
    // ending with a decimal point.
    html << std::scientific;

    computeLayout(options);
    computeLayoutBoundingBox();

    Box viewportBox = boundingBox;
    viewportBox.extend(0.05);
    viewportBox.makeSquare();

    // Begin the svg.
    const string svgId = "LocalssemblyGraph";
    html <<
        "\n<br><div style='display:inline-block;vertical-align:top;'>"
        "<svg id='" << svgId <<
        "' width='" <<  options.sizePixels <<
        "' height='" << options.sizePixels <<
        "' viewbox='" << viewportBox.xMin << " " << viewportBox.yMin << " " <<
        viewportBox.xSize() << " " <<
        viewportBox.ySize() << "'"
        " style='background-color:#f0f0f0'"
        ">\n";

    // Write the edges first so they don't obscure the vertices.
    writeEdges(html, options, assemblyStage);

    // Write the vertices.
    writeVertices(html, options);

    // Finish the svg.
    html << "</svg></div>";

    // Side panel.
    html << "<div style='display:inline-block;margin-left:20px'>";
    writeSvgControls(html, options);
    html << "</div>";

}



void LocalAssemblyGraph::writeSvgControls(
    ostream& html,
    const LocalAssemblyGraphDisplayOptions& /* options */) const
{
    html <<
        "<p><table>";

    // Add drag and zoom.
    addSvgDragAndZoom(html);



    // Buttons to change vertex size.
    html << R"stringDelimiter(
    <tr><th class=left>Vertex size<td>
    <button type='button' onClick='changeVertexSize(0.1)' style='width:3em'>---</button>
    <button type='button' onClick='changeVertexSize(0.5)' style='width:3em'>--</button>
    <button type='button' onClick='changeVertexSize(0.8)' style='width:3em'>-</button>
    <button type='button' onClick='changeVertexSize(1.25)' style='width:3em'>+</button>
    <button type='button' onClick='changeVertexSize(2.)' style='width:3em'>++</button>
    <button type='button' onClick='changeVertexSize(10.)' style='width:3em'>+++</button>
        <script>
        function changeVertexSize(factor)
        {
            var vertexGroup = document.getElementById('vertices');
            var vertices = vertexGroup.getElementsByTagName('circle');
            for(i=0; i<vertices.length; i++) {
                v = vertices[i];
                v.setAttribute('r', factor * v.getAttribute('r'));
            }
        }
        </script>
        )stringDelimiter";



    // Buttons to change edge thickness.
    html << R"stringDelimiter(
    <tr><th class=left>Edge thickness<td>
    <button type='button' onClick='changeThickness(0.1)' style='width:3em'>---</button>
    <button type='button' onClick='changeThickness(0.5)' style='width:3em'>--</button>
    <button type='button' onClick='changeThickness(0.8)' style='width:3em'>-</button>
    <button type='button' onClick='changeThickness(1.25)' style='width:3em'>+</button>
    <button type='button' onClick='changeThickness(2.)' style='width:3em'>++</button>
    <button type='button' onClick='changeThickness(10.)' style='width:3em'>+++</button>
        <script>
        function changeThickness(factor)
        {
            var edgeGroup = document.getElementById('edges');
            var edges = edgeGroup.getElementsByTagName('path');
            for(i=0; i<edges.length; i++) {
                e = edges[i];
                e.setAttribute('stroke-width', factor * e.getAttribute('stroke-width'));
            }

            var arrowsGroup = document.getElementById('arrows');
            var arrows = arrowsGroup.getElementsByTagName('line');
            for(i=0; i<arrows.length; i++) {
                a = arrows[i];
                a.setAttribute('stroke-width', factor * a.getAttribute('stroke-width'));
            }
        }
        </script>
        )stringDelimiter";



    // Zoom buttons.
    html << R"stringDelimiter(
        <tr title='Or use the mouse wheel.'><th class=left>Zoom<td>
        <button type='button' onClick='zoomSvg(0.1)' style='width:3em'>---</button>
        <button type='button' onClick='zoomSvg(0.5)' style='width:3em'>--</button>
        <button type='button' onClick='zoomSvg(0.8)' style='width:3em'>-</button>
        <button type='button' onClick='zoomSvg(1.25)' style='width:3em'>+</button>
        <button type='button' onClick='zoomSvg(2.)' style='width:3em'>++</button>
        <button type='button' onClick='zoomSvg(10.)' style='width:3em'>+++</button>
    )stringDelimiter";



    // Buttons to highlight an anchor and zoom to an anchor.
    html << R"stringDelimiter(
        <tr><td colspan=2>
        <button onClick='highlightAnchor()'>Highlight</button>
        <button onClick='zoomToAnchor()'>Zoom to</button>anchor
        <input id=selectedAnchorId type=text size=10 style='text-align:center'>
    <script>
    function zoomToAnchor()
    {
        // Get the anchor id from the input field.
        var anchorId = document.getElementById("selectedAnchorId").value;
        zoomToGivenAnchor(anchorId);
    }
    function zoomToGivenAnchor(anchorId)
    {
        var element = document.getElementById(anchorId);
        // Find the bounding box and its center.
        var box = element.getBBox();
        var xCenter = box.x + 0.5 * box.width;
        var yCenter = box.y + 0.5 * box.height;

        // Change the viewbox of the svg to be a bit larger than a square
        // containing the bounding box.
        var enlargeFactor = 5.;
        var size = enlargeFactor * Math.max(box.width, box.height);
        var factor = size / width;
        width = size;
        height = size;
        x = xCenter - 0.5 * size;
        y = yCenter - 0.5 * size;
        var svg = document.querySelector('svg');
        svg.setAttribute('viewBox', `${x} ${y} ${size} ${size}`);
        ratio = size / svg.getBoundingClientRect().width;
        svg.setAttribute('font-size', svg.getAttribute('font-size') * factor);
    }
    function highlightAnchor()
    {
        // Get the anchor id  from the input field.
        var anchorId = document.getElementById("selectedAnchorId").value;
        var element = document.getElementById(anchorId);

        element.style.fill = "Magenta";
    }
    </script>
    )stringDelimiter";


    html << "</table>";

    // Scroll down to the svg.
    const string svgId = "LocalAssemblyGraph";
    html <<
        "<script>"
        "document.getElementById('" << svgId << "').scrollIntoView({block:'center'});"
        "</script>";

    html <<
        "<p>Use Ctrl+Click to pan."
        "<p>Use Ctrl-Wheel or the above buttons to zoom.";
}



void LocalAssemblyGraph::writeVertices(
    ostream& html,
    const LocalAssemblyGraphDisplayOptions& options) const
{
    const LocalAssemblyGraph& localAssemblyGraph = *this;

    html << "\n<g id='vertices' style='stroke:none'>";

    BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph) {
        const LocalAssemblyGraphVertex& vertex = localAssemblyGraph[v];
        const AnchorId anchorId = getAnchorId(v);
        const string anchorIdString = anchorIdToString(anchorId);
        const uint64_t coverage = assemblyGraph.anchors[anchorId].coverage();

        // Get the position of this vertex in the computed layout.
        const auto& xy = vertex.xy;
        const double x = xy[0];
        const double y = xy[1];


        // Choose the color for this vertex.
        string color;
        if(vertex.distance == maxDistance) {
            color = "Cyan";
        } else if(vertex.distance == 0) {
            color = "Blue";
        }

        // Hyperlink.
        html << "\n<a href='exploreAnchor?anchorIdString=" <<
            HttpServer::urlEncode(anchorIdString) << "'>";

        // Write the vertex.
        html << "<circle cx='" << x << "' cy='" << y <<
            "' fill='" << color <<
            "' r='" << vertexRadius(v, options) <<
            "' id='" << anchorIdString << "'>"
            "<title>" << anchorIdString << ", coverage " << coverage <<
            "</title></circle>";

        // End the hyperlink.
        html << "</a>";
    }
    html << "\n</g>";
}



double LocalAssemblyGraph::vertexRadius(
    vertex_descriptor v,
    const LocalAssemblyGraphDisplayOptions& options) const
{
    const double scalingFactor =
        (options.layoutMethod == "sfdp") ? 0.01 : 0.05;

    if(options.vertexSizeByCoverage) {
        const AnchorId anchorId = getAnchorId(v);
        const uint64_t coverage = assemblyGraph.anchors[anchorId].coverage();
        return options.vertexSize * (scalingFactor * sqrt(double(coverage)));
    } else {
        return options.vertexSize * (scalingFactor * 3.);
    }
}


void LocalAssemblyGraph::writeEdges(
    ostream& html,
    const LocalAssemblyGraphDisplayOptions& options,
    const string& assemblyStage) const
{
    const LocalAssemblyGraph& localAssemblyGraph = *this;
    const double scalingFactor =
        (options.layoutMethod == "sfdp") ? 0.001 : 0.005;


    html << "\n<g id=edges>";

    BGL_FORALL_EDGES(e, localAssemblyGraph, LocalAssemblyGraph) {
        const LocalAssemblyGraphEdge& edge = localAssemblyGraph[e];
        const string chainStringId = assemblyGraph.getChainStringId(edge);

        const vertex_descriptor v0 = source(e, localAssemblyGraph);
        const vertex_descriptor v1 = target(e, localAssemblyGraph);

        const LocalAssemblyGraphVertex& vertex0 = localAssemblyGraph[v0];
        const LocalAssemblyGraphVertex& vertex1 = localAssemblyGraph[v1];

        const auto& xy0 = vertex0.xy;
        const auto& xy1 = vertex1.xy;

        const double x0 = xy0[0];
        const double y0 = xy0[1];
        const double x1 = xy1[0];
        const double y1 = xy1[1];

        // The middle point of the quadratic Bezier spline used ot display the
        // edge is the xy stored in the edge.
        const auto& xym = edge.xy;
        const double xm = xym[0];
        const double ym = xym[1];


        string color = "Black";

        if(options.edgeColoring == "random") {
            // To decide the color, hash the chainStringId.
            // This way we always get the same color for the same edge.
            const uint32_t hashValue = MurmurHash2(chainStringId.data(), int(chainStringId.size()), 759);
            const uint32_t hue = hashValue % 360;
            color = "hsl(" + to_string(hue) + ",50%,50%)";
        }

        // Hyperlink.
        html << "\n<a href='"
            "exploreSegment?assemblyStage=" << assemblyStage <<
            "&segmentName=" << chainStringId << "'>";

        html <<
            "\n<path d="
            "'M " << x0 << " " << y0 <<
            " Q " << xm << " " << ym <<
            " " << x1 << " " << y1 <<
            "' stroke='" << color <<
            "' stroke-width='" << scalingFactor * options.edgeThickness <<
            "' fill=none>"
            "<title>" << chainStringId << "</title>"
            "</path>";

        // End the hyperlink.
        // html << "</a>";
    }
    html << "</g>";



    // Write the "arrows" to show edge directions.
    // They are just short lines near the target vertex of each edge.
    html << "\n<g id=arrows";
    if(options.edgeColoring == "black") {
        html << " stroke=white";
    } else {
        html << " stroke=black";
    }
    html << ">";
    BGL_FORALL_EDGES(e, localAssemblyGraph, LocalAssemblyGraph) {
        const LocalAssemblyGraphEdge& edge = localAssemblyGraph[e];
        const string chainStringId = assemblyGraph.getChainStringId(edge);

        const vertex_descriptor v0 = source(e, localAssemblyGraph);
        const vertex_descriptor v1 = target(e, localAssemblyGraph);

        const LocalAssemblyGraphVertex& vertex0 = localAssemblyGraph[v0];
        const LocalAssemblyGraphVertex& vertex1 = localAssemblyGraph[v1];

        const auto& xy0 = vertex0.xy;
        const auto& xy1 = vertex1.xy;

        // The positions of v0 and v1.
        const double x0 = xy0[0];
        const double y0 = xy0[1];
        const double x1 = xy1[0];
        const double y1 = xy1[1];

        // The distance between x0 and x1.
        const double d01x = x1 - x0;
        const double d01y = y1 - y0;
        const double d01 = sqrt(d01x * d01x + d01y * d01y);

        // The middle point of the quadratic Bezier spline used to display the
        // edge is the xy stored in the edge.
        const auto& xym = edge.xy;
        const double xm = xym[0];
        const double ym = xym[1];

        // Compute a unit vector in the direction xym - x1.
        const double d1mx = xm - x1;
        const double d1my = ym - y1;
        const double d1m = sqrt(d1mx * d1mx + d1my * d1my);
        const double e1mx = d1mx / d1m;
        const double e1my = d1my / d1m;

        // The "arrow" begins at xy1 and ends at xy2 computed here.
        const double relativeArrowLength = 0.2;
        const double arrowLength = relativeArrowLength * d01;
        const double x2 = x1 + e1mx * arrowLength;
        const double y2 = y1 + e1my * arrowLength;


        html <<
            "\n<line x1='" << x1 << "' y1='" << y1 <<
            "' x2='" << x2 << "' y2='" << y2 <<
            "' stroke-width='" << 0.2 * scalingFactor * options.edgeThickness <<
            "' />";

    }
    html << "</g>";
}



AnchorId LocalAssemblyGraph::getAnchorId(vertex_descriptor v) const
{
    const LocalAssemblyGraph& localAssemblyGraph = *this;
    const LocalAssemblyGraphVertex& vertex = localAssemblyGraph[v];

    if(vertex.isTypeA()) {
        const AssemblyGraph::vertex_descriptor av = vertex.v;
        return assemblyGraph[av].anchorId;
    } else {

        const BubbleChain& bubbleChain = assemblyGraph[vertex.e];
        const Bubble& bubble = bubbleChain[vertex.positionInBubbleChain];
        const Chain& chain = bubble.front();
        return chain.back();
    }
}



void LocalAssemblyGraph::computeLayout(const LocalAssemblyGraphDisplayOptions& options)
{
    LocalAssemblyGraph& localAssemblyGraph = *this;

    // For readably display of bubbles, we want to add an auxiliary vertex in
    // the middle of each vertex for the purpose of layout computation.
    // For this purpose we need to create an auxiliary graph as follows:
    // - The first num_vertices(localAssemblyGraph) correspond to the LocalAssemblyGraph vertices.
    // - The next num_edges(localAssemblyGraph) correspond to the LocalAssemblyGraph edges.

    using AuxiliaryGraph = boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::bidirectionalS>;
    AuxiliaryGraph auxiliaryGraph(num_vertices(localAssemblyGraph) + num_edges(localAssemblyGraph));

    // Vertex descriptors of the LocalAssemblyGraph are identical to vertex descriptors
    // of the AuxiliaryGraph.
    // We need to map edge descriptors of the localAssemblyGraph to vertex descriptors of
    // the auxiliary graph.
    std::map<edge_descriptor, AuxiliaryGraph::vertex_descriptor> edgeMap;
    AuxiliaryGraph::vertex_descriptor i = num_vertices(localAssemblyGraph);
    BGL_FORALL_EDGES(e, localAssemblyGraph, LocalAssemblyGraph) {
        edgeMap.insert({e, i++});
    }
    SHASTA_ASSERT(i == num_vertices(auxiliaryGraph));

    // Add edges of the AuxiliaryGraph.
    // Each edge of the LocalAssemblyGraph generates two edges of the AuxiliaryGraph.
    std::map<AuxiliaryGraph::edge_descriptor, double> edgeLengthMap;
    BGL_FORALL_EDGES(e, localAssemblyGraph, LocalAssemblyGraph) {
        const AuxiliaryGraph::vertex_descriptor v0 = source(e, localAssemblyGraph);
        const AuxiliaryGraph::vertex_descriptor v1 = target(e, localAssemblyGraph);
        const AuxiliaryGraph::vertex_descriptor vMiddle = edgeMap[e];
        const double displayLength =
            options.minimumEdgeLength +
            options.additionalEdgeLengthPerKb * 0.001 * double(getChainLength(e));

        AuxiliaryGraph::edge_descriptor e0;
        tie(e0, ignore) = add_edge(v0, vMiddle, auxiliaryGraph);
        edgeLengthMap.insert({e0, displayLength / 2.});

        AuxiliaryGraph::edge_descriptor e1;
        tie(e1, ignore) = add_edge(vMiddle, v1, auxiliaryGraph);
        edgeLengthMap.insert({e1, displayLength / 2.});
    }


    // Compute the layout of the AuxiliaryGraph.
    std::map<AuxiliaryGraph::vertex_descriptor, array<double, 2> > layout;
    const double timeout = 30.;
    if(options.layoutMethod == "custom") {
        const int quality = 2;
        computeLayoutCustom(
            auxiliaryGraph,
            edgeLengthMap,
            layout,
            quality,
            timeout);
    } else {
        const string additionalOptions = "";
        computeLayoutGraphviz(
            auxiliaryGraph,
            options.layoutMethod,
            timeout,
            layout,
            additionalOptions,
            &edgeLengthMap);
        // If the layout is dot, reverse the y coordinates so the arrows point down.
        if(options.layoutMethod == "dot") {
            for(auto& p: layout) {
                auto& y = p.second[1];
                y = -y;
            }
        }
    }

    // Store the layout in the vertices.
    BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph) {
        localAssemblyGraph[v].xy = layout[v];
    }


    // Store the layout in the edges.
    BGL_FORALL_EDGES(e, localAssemblyGraph, LocalAssemblyGraph) {
        const auto it = edgeMap.find(e);
        SHASTA_ASSERT(it != edgeMap.end());
        const AuxiliaryGraph::vertex_descriptor vAux = it->second;

        const auto jt = layout.find(vAux);
        SHASTA_ASSERT(jt != layout.end());
        localAssemblyGraph[e].xy = jt->second;
    }
}



uint64_t LocalAssemblyGraph::getChainLength(edge_descriptor e) const
{
    const LocalAssemblyGraph& localAssemblyGraph = *this;
    const LocalAssemblyGraphEdge& edge = localAssemblyGraph[e];
    const Chain& chain = assemblyGraph.getChain(edge);

    if(assemblyGraph.sequenceWasAssembled) {
        return chain.sequence.size();
    } else {
        return assemblyGraph.chainOffset(chain);
    }
}



void LocalAssemblyGraph::computeLayoutBoundingBox()
{
    const LocalAssemblyGraph& localAssemblyGraph = *this;

    boundingBox.xMin = std::numeric_limits<double>::max();
    boundingBox.xMax = std::numeric_limits<double>::min();
    boundingBox.yMin = boundingBox.xMin;
    boundingBox.yMax = boundingBox.xMax;

    BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph) {
        const array<double, 2>& xy = localAssemblyGraph[v].xy;
        const double x = xy[0];
        const double y = xy[1];
        boundingBox.xMin = min(boundingBox.xMin, x);
        boundingBox.xMax = max(boundingBox.xMax, x);
        boundingBox.yMin = min(boundingBox.yMin, y);
        boundingBox.yMax = max(boundingBox.yMax, y);
    }

    BGL_FORALL_EDGES(e, localAssemblyGraph, LocalAssemblyGraph) {
        const array<double, 2>& xy = localAssemblyGraph[e].xy;
        const double x = xy[0];
        const double y = xy[1];
        boundingBox.xMin = min(boundingBox.xMin, x);
        boundingBox.xMax = max(boundingBox.xMax, x);
        boundingBox.yMin = min(boundingBox.yMin, y);
        boundingBox.yMax = max(boundingBox.yMax, y);
    }

}



void LocalAssemblyGraph::Box::makeSquare()
{
    if(xSize() > ySize()) {
        const double delta = (xSize() - ySize()) / 2.;
        yMin -= delta;
        yMax += delta;
    } else {
        const double delta = (ySize() - xSize()) / 2.;
        xMin -= delta;
        xMax += delta;
    }
}



void LocalAssemblyGraph::Box::extend(double factor)
{
    const double extend = factor * max(xSize(), ySize());
    xMin -= extend;
    xMax += extend;
    yMin -= extend;
    yMax += extend;
}
