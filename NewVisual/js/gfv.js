// GLOBAL
var padding = 100;

// exon graphic properties
var exonStartSize = 20;
var exonGrowPerEdge = 5;

// make SVG Container for the visualization
var svg = d3.select("body").append("svg").attr("width", screen.width).attr("height", screen.height);


// create different goups for exons and edges so that they can render in correct order
var edgesLayer = svg.append('g');
var exonLayer = svg.append('g');

// array of exon objects
var exons = [];

exons.push(new createExon(1, 100));
exons.push(new createExon(3, 100));
exons.push(new createExon(3, 300));
exons.push(new createExon(5, 100));
exons.push(new createExon(5, 300));

createEdge(exons[0], exons[1]);
createEdge(exons[0], exons[1]);
createEdge(exons[0], exons[2]);
createEdge(exons[2], exons[3]);
createEdge(exons[0], exons[1]);
createEdge(exons[0], exons[1]);
createEdge(exons[3], exons[4]);
// createEdge(exons[0], exons[1]);
// var test = createEdge(exons[1], exons[2]);




// function to create an exon
function createExon(hPos, ypos) {
    // positions and dimensions of exon
    var xPos = hPos * padding;
    var yPos = ypos;

    // radius of exonGraphic
    var radius = 20;

    // arrays of exons connected by incoming and outgoing edges
    this.inExons = [];
    this.outExons = [];

    // arrays of exons connected by incoming and outgoing edges
    this.inEdges = [];
    this.outEdges = [];

    // create the exonGraphic
    this.exonGraphic = exonLayer.append("ellipse")
                                .attr("cx", xPos)
                                .attr("cy", yPos)
                                .attr("rx", radius)
                                .attr("ry", radius)
                                .style("stroke-width", 2)
                                .style("stroke", "black")
                                .style("fill", "lightblue");

    // this.exonGraphic.style("stroke-dasharray", ("10, 3"))
    return this;
}


// function to make an edge between two exons
function createEdge (exon1, exon2) {
    // get positions
    var x1 = exon1.exonGraphic.attr("cx");
    var y1 = exon1.exonGraphic.attr("cy");
    var x2 = exon2.exonGraphic.attr("cx");
    var y2 = exon2.exonGraphic.attr("cy");

    // add the to exons to in/out lists
    exon1.outExons.push(exon2);
    exon2.outExons.push(exon1);


    // create the line (edge graphic)
    var edge = edgesLayer.append("line")
                  .attr("x1", x1)
                  .attr("y1", y1)
                  .attr("x2", x2)
                  .attr("y2", y2)
                  .attr("stroke-width", 2)
                  .attr("stroke", "black");

    // add the edge to the exons list
    exon1.outEdges.push(edge);
    exon2.inEdges.push(edge);

    // increase the size of the exons based on the number of edges
    exon1.exonGraphic.attr("rx", exonStartSize + Math.max(exon1.inEdges.length, exon1.outEdges.length) * exonGrowPerEdge);
    exon1.exonGraphic.attr("ry", exonStartSize + Math.max(exon1.inEdges.length, exon1.outEdges.length) * exonGrowPerEdge);
    exon2.exonGraphic.attr("rx", exonStartSize + Math.max(exon2.inEdges.length, exon2.outEdges.length) * exonGrowPerEdge);
    exon2.exonGraphic.attr("ry", exonStartSize + Math.max(exon2.inEdges.length, exon2.outEdges.length) * exonGrowPerEdge);

    // the margin between the group of edges and the top and bottom of exons
    var edgeMargin = 20;

    // if there is more than one outgoing edge, determine where the edges should be positioned
    if (exon1.outEdges.length > 1) {
        // get the coordinates for the top of the exon
        var topEdgeHeight = exon1.exonGraphic.attr("cy") - (exon1.exonGraphic.attr("ry"));
        // calculae the proper space between each edge
        var edgePadding = (exon1.exonGraphic.attr("ry") * 2 - (edgeMargin*2) ) / (exon1.outEdges.length-1);

        // loop through all outgoing edges and set position
        for (var i = 0; i < exon1.outEdges.length; i++) {
            exon1.outEdges[i].attr("y1", topEdgeHeight + i * edgePadding + edgeMargin);
        }
    }

    if (exon2.inEdges.length > 1) {
        // do the same for the exon2's incoming edges
        topEdgeHeight = exon2.exonGraphic.attr("cy") - (exon2.exonGraphic.attr("ry"));
        edgePadding = (exon2.exonGraphic.attr("ry") * 2 - (edgeMargin*2) ) / (exon2.inEdges.length-1);

        // loop through all incoming edges and set position
        for (var i = 0; i < exon2.inEdges.length; i++) {
            exon2.inEdges[i].attr("y2", topEdgeHeight + i * edgePadding + edgeMargin);
        }
    }

    return edge;
}
