#include "html.hpp"
using namespace shasta;

#include "iostream.hpp"



void shasta::writeHtmlBegin(ostream& html, const string& title)
{
    html <<
        "<!DOCTYPE html>"
        "<html>"
        "<head>"
        "<meta charset='UTF-8'>"
        "<title>" << title << "</title>";
    writeStyle(html);
    html << "</head>";
}



void shasta::writeHtmlEnd(ostream& html)
{
    html << "</html>";
}



void shasta::writeStyle(ostream& html)
{
    html << R"%(
<style>
    body {
        font-family: Arial;
    }
    pre {
        font-family: courier;
    }
    p, input {
        font-size: 16px;
    }
    h1, h2, h3 {
        color: DarkSlateBlue;
    }
    table {
        border-collapse: collapse;
    }
    th, td {
        border: 1px solid #b8b5c7d9;
        padding: 2px;
    }
    th {
        font-weight: bold;
        text-align: center;
    }
    th.left {
        text-align: left;
    }
    td.centered {
        text-align: center;
    }
    td.right {
        text-align: right;
    }
    a {
        color: DarkSlateBlue;
    }

    /* This can be used to get vertical text in table cells. */
    span.rotated 
    {
      writing-mode: vertical-rl;
      transform: rotate(180deg);
    }
</style>
    )%";

}



void shasta::addSvgDragAndZoom(ostream& html)
{
    html << R"zzz(  
<script>

var svg = document.querySelector('svg');
svg.scrollIntoView();
svg.addEventListener('pointerdown', onPointerDown); 
svg.addEventListener('pointerup', onPointerUp); 
svg.addEventListener('pointerleave', onPointerLeave); 
svg.addEventListener('pointermove', onPointerMove); 
svg.addEventListener('wheel', onMouseWheel); 

var pointerIsDown = false;

var xOrigin = 0;
var yOrigin = 0;

// The current viewbox.
var x = svg.viewBox.baseVal.x;
var y = svg.viewBox.baseVal.y;
var width = svg.viewBox.baseVal.width;
var height = svg.viewBox.baseVal.height;

var xNew = 0;
var yNew = 0;

var ratio = width / svg.getBoundingClientRect().width;

function onPointerDown(event) {
    // If the CTRL key is not pressed, don't do anything.
    if(!event.ctrlKey) {
        return;
    }
    event.stopPropagation();
    event.preventDefault();
    pointerIsDown = true;
    xOrigin = event.clientX;
    yOrigin = event.clientY;

    return false;
}

function onPointerMove (event) {
    // If the CTRL key is not pressed, don't do anything.
    if(!event.ctrlKey) {
        return;
    }
    event.stopPropagation();
    event.preventDefault();
    if (!pointerIsDown) {
        return;
    }

    xNew = x - (event.clientX - xOrigin) * ratio;
    yNew = y - (event.clientY - yOrigin) * ratio;

    svg.setAttribute('viewBox', `${xNew} ${yNew} ${width} ${height}`);

    return false;
}

function onPointerUp(event) {
    // If the CTRL key is not pressed, don't do anything.
    if(!event.ctrlKey) {
        return;
    }
    event.stopPropagation();
    event.preventDefault();
    pointerIsDown = false;
    x = xNew;
    y = yNew;

    return false;
}

function onPointerLeave(event) {
    // If the CTRL key is not pressed, don't do anything.
    if(!event.ctrlKey) {
        return;
    }
    event.stopPropagation();
    event.preventDefault();
    if(pointerIsDown) {
        pointerIsDown = false;
        x = xNew;
        y = yNew;
    }

    return false;
}

function onMouseWheel(event) {
    // If the CTRL key is not pressed, don't do anything.
    if(!event.ctrlKey) {
        return;
    }
    event.stopPropagation();
    event.preventDefault();  
    var value = event.wheelDelta / 120.;
    var factor = Math.pow(1.1, value);
    zoomSvg(factor);
}

function zoomSvg(factor)
{
    ratio /= factor;

    // Adjust the viewbox so the center does not move.
    var xCenter = x + 0.5 * width;
    var yCenter = y + 0.5 * height;  
    width /= factor;
    height /= factor;
    x = xCenter - 0.5 * width;
    y = yCenter - 0.5 * height;

    svg.setAttribute('viewBox', `${x} ${y} ${width} ${height}`);
    svg.setAttribute('font-size', svg.getAttribute('font-size') / factor);

    return false;
}

</script>
    )zzz";
}




void shasta::writeInformationIcon(ostream& html, const string& message)
{
    html << "<span style='color:Blue;font-weight:bold' title=\"" <<
            message << "\">&#9432;</span>";
}
