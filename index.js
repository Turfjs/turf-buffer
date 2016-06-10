var simplepolygon = require('../simplepolygon');
var destination = require('turf-destination');
var bearing = require('turf-bearing');
var helpers = require('turf-helpers');
var union = require('turf-union');
var difference = require('turf-difference');

module.exports = function(feature, radius, units, resolution){
  if (!resolution) resolution = 32;
  if (radius <= 0) throw new Error("The buffer radius must be positive");
  var geom = feature.geometry;
  if (geom === null) return feature;
  if(geom.type === 'Point') {
    return pointBuffer(feature, radius, units, resolution);
  } else if(geom.type === 'MultiPoint') {
    var buffers = [];
    geom.coordinates.forEach(function(coords) {
      buffers.push(pointBuffer(helpers.point(coords), radius, units, resolution));
    });
    return helpers.featureCollection(buffers)
  } else if(geom.type === 'LineString') {
    return lineBuffer(feature, radius, units, resolution);
  } else if(geom.type === 'MultiLineString') {
    var buffers = [];
    geom.coordinates.forEach(function(coords) {
      buffers.push(lineBuffer(helpers.lineString(coords), radius, units, resolution));
    });
    return helpers.featureCollection(buffers)
  } else if(geom.type === 'Polygon') {
    return polygonBuffer(feature, radius, units, resolution);
  } else if(geom.type === 'MultiPolygon') {
    var buffers = [];
    geom.coordinates.forEach(function(coords) {
      buffers.push(polygonBuffer(helpers.polygon(coords), radius, units, resolution));
    });
    return helpers.featureCollection(buffers)
  }
}

function pointBuffer(pt, radius, units, resolution) {
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

function lineBuffer(line, radius, units, resolution) {
  var lineOffset = [];

  if (!(line.geometry.coordinates[0].equals(line.geometry.coordinates[line.geometry.coordinates.length-1]))) {

    // situation at current point = point 0
    var currentLinePoint = helpers.point(line.geometry.coordinates[0]);
    var nextLineBearing = bearing(helpers.point(line.geometry.coordinates[0]), helpers.point(line.geometry.coordinates[1]));
    var currentBufferPoint = destination(currentLinePoint, radius, nextLineBearing + 90, units);
    var previousLinePoint = helpers.point(line.geometry.coordinates[line.geometry.coordinates.length-1]);
    var previousLineBearing = bearing(helpers.point(line.geometry.coordinates[line.geometry.coordinates.length-2]), helpers.point(line.geometry.coordinates[line.geometry.coordinates.length-1]));

    lineOffset.push.apply(lineOffset,[currentBufferPoint.geometry.coordinates]); // Add first buffer point in order to close ring
    lineOffset.push.apply(lineOffset,lineOffsetOneSide(line, radius, units, resolution, false, true).geometry.coordinates);
    lineOffset.push.apply(lineOffset,arc(previousLinePoint, radius, previousLineBearing + 90, previousLineBearing - 90, units, resolution, true).geometry.coordinates);
    lineOffset.push.apply(lineOffset,lineOffsetOneSide(line, radius, units, resolution, true, true).geometry.coordinates);
    lineOffset.push.apply(lineOffset,arc(currentLinePoint, radius, nextLineBearing - 90, nextLineBearing + 90, units, resolution, true).geometry.coordinates);

    return offsetToBuffer(helpers.polygon([lineOffset]));

  } else {

    lineOffset.push(ringOffsetOneSide(line, radius, units, resolution, false, true).geometry.coordinates);
    lineOffset.push(ringOffsetOneSide(line, radius, units, resolution, true, true).geometry.coordinates);

    return offsetToBuffer(helpers.polygon(lineOffset));
  }
}

function polygonBuffer(poly, radius, units, resolution) {
  var polygonOffset = [];

  polygonOffset.push(ringOffsetOneSide(helpers.lineString(poly.geometry.coordinates[0]), radius, units, resolution, false, true).geometry.coordinates);
  for (var i = 1; i < poly.geometry.coordinates.length; i++) {
    polygonOffset.push(ringOffsetOneSide(helpers.lineString(poly.geometry.coordinates[i]), radius, units, resolution, false, true).geometry.coordinates);
  }

  return offsetToBuffer(helpers.polygon(polygonOffset));
}

function lineOffsetOneSide(line, radius, units, resolution, reverse, right) {
  if (reverse === undefined) var reverse = false;
  if (right === undefined) var right = true;
  if (reverse) line.geometry.coordinates = line.geometry.coordinates.reverse();
  var coords = line.geometry.coordinates;
  var lineOffset = [];
  if (coords.length == 2) return helpers.lineString(lineOffset)

  var currentLinePoint = helpers.point(coords[1]);
  var previousLineBearing = bearing(helpers.point(coords[0]), helpers.point(coords[1]));
  for (var i = 1; i < coords.length-1; i++) {
    var nextLinePoint = helpers.point(coords[i+1]);
    var nextLineBearing = bearing(currentLinePoint, nextLinePoint);
    lineOffset.push.apply(lineOffset, arc(currentLinePoint, radius, previousLineBearing + Math.pow(-1, right + 1) * 90, nextLineBearing + Math.pow(-1, right + 1) * 90, units, resolution, right, true).geometry.coordinates);
    var currentLinePoint = nextLinePoint;
    var previousLineBearing = nextLineBearing;
  }

  return helpers.lineString(lineOffset)
}

function ringOffsetOneSide(ring, radius, units, resolution, reverse, right) {
  if (reverse === undefined) var reverse = false;
  if (right === undefined) var right = true;
  if (reverse) ring.geometry.coordinates = ring.geometry.coordinates.reverse();
  var coords = ring.geometry.coordinates; // ring is a linestring
  var ringOffset = [];

  // situation at current point = point 0
  var currentRingPoint = helpers.point(coords[0]);
  var nextRingBearing = bearing(helpers.point(coords[0]), helpers.point(coords[1]));
  var currentBufferPoint = destination(currentRingPoint, radius, nextRingBearing + 90, units);
  var previousRingPoint = helpers.point(coords[coords.length-1]);
  var previousRingBearing = bearing(helpers.point(coords[coords.length-2]), helpers.point(coords[coords.length-1]));

  ringOffset.push.apply(ringOffset, [currentBufferPoint.geometry.coordinates]); // Add first buffer point in order to close ring
  ringOffset.push.apply(ringOffset, lineOffsetOneSide(ring, radius, units, resolution, false, right).geometry.coordinates);
  ringOffset.push.apply(ringOffset, arc(currentRingPoint, radius, previousRingBearing + Math.pow(-1, right + 1) * 90, nextRingBearing + Math.pow(-1, right + 1) * 90, units, resolution, right, true).geometry.coordinates);

  return helpers.lineString(ringOffset)
}

function arc(pt, radius, bearing1, bearing2, units, resolution, right, shortcut) {
  if (right === undefined) var right = true;
  if (shortcut === undefined) var shortcut = false;
  var arc = [];
  var resMultiple = 360/resolution;
  if (right) {
      var bearing = Math.floor(bearing1/resMultiple)*resMultiple;
  } else {
    var bearing = Math.ceil(bearing1/resMultiple)*resMultiple;
  }
  var angle = (Math.pow(-1, right + 1) * (bearing1 - bearing2)).mod(360);
  var numSteps = Math.ceil((Math.pow(-1, right + 1) * (bearing - bearing2)).mod(360)/resMultiple);
  var step = numSteps;
  // Add spoke for bearing1
  var spoke = destination(pt, radius, bearing1, units);
  arc.push(spoke.geometry.coordinates);
  if (!(angle > 180 && shortcut)) {
    while (step) {
      // Only add spokes for bearings different than the first and last
      if (!(bearing == bearing1) || (bearing == bearing2)) {
        spoke = destination(pt, radius, bearing, units);
        arc.push(spoke.geometry.coordinates);
      }
      // The value of bearing is independent of bearing1 and bearing2, such that multiple arcs with the same centerpoint coincide, instead of zigzag overlapping
      bearing = bearing + Math.pow(-1, !right + 1) * resMultiple;
      step--;
    }
  }
  // Only add spoke for bearing2 if it's not a straight line (equivalently: angle == 0)
  if (!((bearing == bearing1) && (bearing1 == bearing2))) {
    spoke = destination(pt, radius, bearing2, units);
    arc.push(spoke.geometry.coordinates);
  }
  return helpers.lineString(arc)
}

function filterNetWinding(fc, filterFn) {
  var i = fc.features.length;
  while (i--) {
    if (!filterFn(fc.features[i].properties.netWinding)) {
        fc.features.splice(i, 1);
    }
  }
  return fc;
}

function unionFeatureCollection(fc) {
  if (fc.features.length == 0) return {type: "Feature", geometry: null};
  var incrementalUnion = fc.features[0];
  if (fc.features.length == 1) return incrementalUnion;
  for (var i = 1; i < fc.features.length; i++) {
    incrementalUnion = union(incrementalUnion, fc.features[i]);
  }
  return incrementalUnion
}

function offsetToBuffer(polygonOffset) {
  var unionWithWindingOne = unionFeatureCollection(filterNetWinding(simplepolygon(polygonOffset), function (netWinding){return netWinding == 1}));
  var unionWithWindingZero = unionFeatureCollection(filterNetWinding(simplepolygon(polygonOffset), function (netWinding){return netWinding == 0}));
  // This last one might have winding -1, so we might have to rewind it if the difference algorithm requires so

  if (unionWithWindingOne.geometry == null) return {type: "Feature", geometry: null};
  if (unionWithWindingZero.geometry == null) return unionWithWindingOne;
  return difference(unionWithWindingOne, unionWithWindingZero);
}

function winding(poly){
  // compute winding of first ring
  var coords = poly.geometry.coordinates[0];
  var leftVtx = 0;
  for (var i = 0; i < coords.length-1; i++) { if (coords[i][0] < coords[leftVtx][0]) leftVtx = i; }
  return (coords[(leftVtx-1).mod(coords.length-1)][1] > coords[(leftVtx+1).mod(coords.length-1)][1]) ? 1 : -1;
}


// Function to compare Arrays of numbers. From http://stackoverflow.com/questions/7837456/how-to-compare-arrays-in-javascript
// Warn if overriding existing method
// if(Array.prototype.equals) console.warn("Overriding existing Array.prototype.equals. Possible causes: New API defines the method, there's a framework conflict or you've got double inclusions in your code.");
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
