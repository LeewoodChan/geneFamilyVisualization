// get the of the screen
var center = screen.height / 2.5;

// space between exons
var padding = 150;
var yPadding = 100;

// exon graphic properties
var exonStartSize = 20;
var exonGrowPerEdge = 3;
var hoverGrow = 20;

// make SVG Container for the visualization
var svg = d3.select("body").append("svg").attr("width", screen.width).attr("height", screen.height);

// create different goups for exons and edges so that they can render in correct order
var edgesLayer = svg.append('g');
var exonLayer = svg.append('g');

// array of exon objects
var exons = [];


// 2D array of exons at each x position
var exonsAtX = new Array(50);
for (var i = 0; i < 20; i++) {
  exonsAtX[i] = new Array();
}


// create exons
addExon(1);
addExon(1);
addExon(1);
addExon(1);
addExon(2);
addExon(2);
addExon(2);
addExon(2);
addExon(2);
addExon(2);
addExon(2);
addExon(3);
addExon(3);
addExon(3);
addExon(4);
addExon(4);
addExon(5);
addExon(6);
addExon(7);
addExon(8);
addExon(8);
addExon(8);
addExon(8);
addExon(8);
addExon(8);
addExon(8);
addExon(8);

addGeneFamily([0, 4, 11, 14, 16, 17, 18, 19, 20, 22, 24], 'blue');
addGeneFamily([5, 11, 15, 16, 17, 18, 19, 20, 22, 24], 'red');
addGeneFamily([6, 12, 15, 16, 17, 18, 19, 20, 22, 24], 'orange');
addGeneFamily([1, 7, 12, 15, 16, 17, 18, 19, 21, 24], 'lightblue');
addGeneFamily([2, 8, 12, 15, 16, 17, 18, 19, 21, 23, 24], 'green');
addGeneFamily([3, 9, 12, 15, 16, 17, 18, 19, 21, 23, 24], 'brown');


// console.log(exons);


// createEdge(exons[0], exons[1]);
// createEdge(exons[0], exons[1]);
// createEdge(exons[0], exons[2]);
// createEdge(exons[2], exons[3]);
// createEdge(exons[1], exons[3]);
// createEdge(exons[0], exons[1]);
// var test = createEdge(exons[1], exons[2]);



svg.append("text")
        .attr("x", (screen.width / 2))
        .attr("y", 0 + (30))
        .attr("text-anchor", "middle")
        .style("font-size", "19px")
        .style("font-family", "Futura")
        .style("text-decoration", "underline")
        .text("Gene Family Visualization");



// function to add an exon to the graphi
function addExon(x) {
    exons.push(new Exon(x));
}



// funciton to create an exon family
function addGeneFamily (exonNumbers, color) {
  for (var i = 0; i < exonNumbers.length-2; i++) {
    createEdge(exons[exonNumbers[i]], exons[exonNumbers[i+1]]);
    // createEdge(exons[i], exons[i + 1]);
  }
}

// function to create an exon
function Exon(x, ypos) {
    // positions and dimensions of exon
    var xPos = x * padding;
    var yPos = center;

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
                                .style("fill", "white")
                                .on('mouseenter', function(d){
                                    d3.select(this).style("fill", "lightblue")
                                    d3.select(this).transition().attr("ry", parseFloat(d3.select(this).attr("ry")) + hoverGrow)
                                                                .attr("rx", parseFloat(d3.select(this).attr("rx")) + hoverGrow)
                                                                .duration(1.4);
                                })
                                .on('mouseleave', function(d){
                                    d3.select(this).style("fill", "white")
                                    d3.select(this).transition().attr("ry", parseFloat(d3.select(this).attr("ry")) - hoverGrow)
                                                                .attr("rx", parseFloat(d3.select(this).attr("rx")) - hoverGrow)
                                });



    // this.exonGraphic.style("stroke-dasharray", ("10, 3"))


    // add exon to ExonAtX array
    exonsAtX[x].push(this);

    // TODO:

    // if there are more than one exons at an x position, adjust positions of exons
    if (exonsAtX.length > 1) {
        var totalHeight = yPadding * (exonsAtX[x].length-1);
        var topExonPos = center - (totalHeight * .5);
        console.log('totalHeight', totalHeight);
        console.log('center', center);
        console.log('topExon', topExonPos);

        // reposition all of the exons at proper positions
        for (var i = 0; i < exonsAtX[x].length; i++) {
            exonsAtX[x][i].exonGraphic.attr("cy", topExonPos + (yPadding * i));
            console.log(exonsAtX[x][i].exonGraphic.attr("cy"));

        }
    }

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
