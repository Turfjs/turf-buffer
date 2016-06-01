var featurecollection = require('turf-featurecollection');
var destination = require('turf-destination');
var bearing = require('turf-bearing');
var helpers = require('turf-helpers');

module.exports = function(feature, radius, units, resolution){
  if(!resolution) resolution = 32;
  var geom = feature.geometry
  if(geom.type === 'Point') {
    return pointBuffer(feature, radius, units, resolution);
  } else if(geom.type === 'MultiPoint') {
    var buffers = [];
    geom.coordinates.forEach(function(coords) {
      buffers.push(pointBuffer(helpers.point(coords[0], coords[1]), radius, units, resolution));
    });
    return featurecollection(buffers)
  } else if(geom.type === 'LineString') {
    return lineBuffer(feature, radius, units, resolution);
  } else if(geom.type === 'MultiLineString') {
    var buffers = [];
    geom.coordinates.forEach(function(line){
      buffers.push(lineBuffer(feature, radius, units, resolution));
    });
  } else if(geom.type === 'Polygon') {

  } else if(geom.type === 'MultiPolygon') {

  }
}

function pointBuffer (pt, radius, units, resolution) {
  var ring = []
  var resMultiple = 360/resolution;
  for(var i  = 0; i < resolution; i++) {
    var spoke = destination(pt, radius, i*resMultiple, units);
    ring.push(spoke.geometry.coordinates);
  }
  if((ring[0][0] !== ring[ring.length-1][0]) && (ring[0][1] != ring[ring.length-1][1])) {
    ring.push([ring[0][0], ring[0][1]]);
  }
  return helpers.polygon([ring])
}

function lineBuffer (line, radius, units, resolution) {

  var lineBuffer = [];
  var firstLinePoint = helpers.point(line.geometry.coordinates[0]);
  var firstLineBearing = bearing(helpers.point(line.geometry.coordinates[0]), helpers.point(line.geometry.coordinates[1]));
  var firstBufferPoint = destination(firstLinePoint, radius, firstLineBearing + 90, units);
  var lastLinePoint = helpers.point(line.geometry.coordinates[line.geometry.coordinates.length-1]);
  var lastLineBearing = bearing(helpers.point(line.geometry.coordinates[line.geometry.coordinates.length-2]), helpers.point(line.geometry.coordinates[line.geometry.coordinates.length-1]));
  lineBuffer.push.apply(lineBuffer,[firstBufferPoint.geometry.coordinates]);
  lineBuffer.push.apply(lineBuffer,lineBufferOneSide(line, radius, units, resolution, true, false).geometry.coordinates);
  lineBuffer.push.apply(lineBuffer,arc(lastLinePoint, radius, lastLineBearing + 90, lastLineBearing - 90, units, resolution, true).geometry.coordinates);
  lineBuffer.push.apply(lineBuffer,lineBufferOneSide(line, radius, units, resolution, true, true).geometry.coordinates);
  lineBuffer.push.apply(lineBuffer,arc(firstLinePoint, radius, firstLineBearing - 90, firstLineBearing + 90, units, resolution, true).geometry.coordinates);

  return helpers.polygon([lineBuffer]);
  /*
  var lineBuffers = featurecollection([])
  //break line into segments
  var segments = [];
  for(var i = 0; i < line.geometry.coordinates.length-1; i++) {
    segments.push([line.geometry.coordinates[i], line.geometry.coordinates[i+1]]);
  }
  */
  /*create a set of boxes parallel to the segments

    ---------

 ((|¯¯¯¯¯¯¯¯¯|))
(((|---------|)))
 ((|_________|))

  */
  /*
  for(var i = 0; i < segments.length; i++) {
    var bottom = point(segments[i][0][0], segments[i][0][1])
    var top = point(segments[i][1][0], segments[i][1][1])

    var direction = bearing(bottom, top);

    var bottomLeft = destination(bottom, radius, direction - 90, units);
    var bottomRight = destination(bottom, radius, direction + 90, units);
    var topLeft = destination(top, radius, direction - 90, units);
    var topRight = destination(top, radius, direction + 90, units);

    var poly = polygon([[bottomLeft.geometry.coordinates, topLeft.geometry.coordinates]]);

    // add top curve
    var spokeNum = Math.floor(resolution/2);
    var topStart = bearing(top, topLeft);
    for(var k = 1; k < spokeNum; k++) {
      var spokeDirection = topStart + (180 * (k/spokeNum))
      var spoke = destination(top, radius, spokeDirection, units);
      poly.geometry.coordinates[0].push(spoke.geometry.coordinates);
    }
    // add right edge
    poly.geometry.coordinates[0].push(topRight.geometry.coordinates)
    poly.geometry.coordinates[0].push(bottomRight.geometry.coordinates)
    //add bottom curve
    var bottomStart = bearing(bottom, bottomRight);
    for(var k = 1; k < spokeNum; k++) {
      var spokeDirection = (bottomStart + (180 * (k/spokeNum)))
      var spoke = destination(bottom, radius, spokeDirection, units);
      poly.geometry.coordinates[0].push(spoke.geometry.coordinates);
    }

    poly.geometry.coordinates[0].push(bottomLeft.geometry.coordinates)
    lineBuffers.features.push(poly);
  }
  console.log(JSON.stringify(lineBuffers))
  return lineBuffers;
  */
}

function lineBufferOneSide (line, radius, units, resolution, right, reverse) {
  if(right === undefined) var right = true;
  if(reverse === undefined) var reverse = false;
  var coords = line.geometry.coordinates;
  if(reverse) coords = coords.reverse();
  var lineBuffer = [];
  if (coords.length == 2) return lineBuffer;
  var currentLinePoint = helpers.point(coords[1]);
  var currentLineBearing = bearing(helpers.point(coords[0]), helpers.point(coords[1]));
  for (var i = 1; i < coords.length-1; i++) {
    var nextLinePoint = helpers.point(coords[i+1]);
    var nextLineBearing = bearing(currentLinePoint, nextLinePoint);
    lineBuffer.push.apply(lineBuffer,arc(currentLinePoint, radius, currentLineBearing + Math.pow(-1, right + 1) * 90, nextLineBearing + Math.pow(-1, right + 1) * 90, units, resolution, right).geometry.coordinates);
    var currentLinePoint = nextLinePoint;
    var currentLineBearing = nextLineBearing;
  }
  return helpers.lineString(lineBuffer)
}

function arc (pt, radius, bearing1, bearing2, units, resolution, right) {
  if(right === undefined) var right = true;
  var arc = [];
  var resMultiple = 360/resolution;
  var steps = Math.floor((Math.pow(-1, right + 1) * (bearing1 - bearing2)).mod(360)/resMultiple);
  for(var i  = 0; i < steps; i++) {
    var spoke = destination(pt, radius, bearing1 + Math.pow(-1, !right + 1) * i * resMultiple, units);
    arc.push(spoke.geometry.coordinates);

  }
  if(bearing1 - i * resMultiple != bearing2) {
    var spoke = destination(pt, radius, bearing2, units);
    arc.push(spoke.geometry.coordinates);
  }
  return helpers.lineString(arc)
}




// Function to compare Arrays of numbers. From http://stackoverflow.com/questions/7837456/how-to-compare-arrays-in-javascript
// Warn if overriding existing method
if(Array.prototype.equals)
    console.warn("Overriding existing Array.prototype.equals. Possible causes: New API defines the method, there's a framework conflict or you've got double inclusions in your code.");
// attach the .equals method to Array's prototype to call it on any array
Array.prototype.equals = function (array) {
    // if the other array is a falsy value, return
    if (!array)
        return false;

    // compare lengths - can save a lot of time
    if (this.length != array.length)
        return false;

    for (var i = 0, l=this.length; i < l; i++) {
        // Check if we have nested arrays
        if (this[i] instanceof Array && array[i] instanceof Array) {
            // recurse into the nested arrays
            if (!this[i].equals(array[i]))
                return false;
        }
        else if (this[i] != array[i]) {
            // Warning - two different object instances will never be equal: {x:20} != {x:20}
            return false;
        }
    }
    return true;
}
// Hide method from for-in loops
Object.defineProperty(Array.prototype, "equals", {enumerable: false});

// Fix Javascript modulo for negative number. From http://stackoverflow.com/questions/4467539/javascript-modulo-not-behaving
Number.prototype.mod = function(n) {
    return ((this%n)+n)%n;
}
