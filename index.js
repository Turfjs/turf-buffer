var simplepolygon = require('simplepolygon');
var destination = require('turf-destination');
var bearing = require('turf-bearing');
var helpers = require('turf-helpers');
var union = require('turf-union');
var difference = require('turf-difference');

module.exports = function(feature, radius, units, resolution){
  if (!resolution) resolution = 32; // Same value as JSTS
  if (radius < 0) throw new Error("The buffer radius must be positive");
  if (radius == 0) return feature;
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
  var pointOffset = [[]];
  var resMultiple = 360/resolution;
  for(var i  = 0; i < resolution; i++) {
    var spoke = destination(pt, radius, i*resMultiple, units);
    pointOffset[0].push(spoke.geometry.coordinates);
  }
  if(!(equalArrays(pointOffset[0][0],pointOffset[0][pointOffset[0].length-1]))) {
    pointOffset[0].push(pointOffset[0][0]);
  }
  return helpers.polygon(pointOffset)
}

function lineBuffer(line, radius, units, resolution) {
  var lineOffset = [];

  line.geometry.coordinates = removeDuplicates(line.geometry.coordinates);

  if (!(equalArrays(line.geometry.coordinates[0],line.geometry.coordinates[line.geometry.coordinates.length-1]))) {

    // situation at current point = point 0
    var currentLinePoint = helpers.point(line.geometry.coordinates[0]);
    var nextLineBearing = bearing(helpers.point(line.geometry.coordinates[0]), helpers.point(line.geometry.coordinates[1]));
    var currentBufferPoint = destination(currentLinePoint, radius, nextLineBearing + 90, units);
    var previousLinePoint = helpers.point(line.geometry.coordinates[line.geometry.coordinates.length-1]);
    var previousLineBearing = bearing(helpers.point(line.geometry.coordinates[line.geometry.coordinates.length-2]), helpers.point(line.geometry.coordinates[line.geometry.coordinates.length-1]));

    lineOffset.push([]);
    lineOffset[0].push.apply(lineOffset[0],[currentBufferPoint.geometry.coordinates]); // Add first buffer point in order to close ring
    lineOffset[0].push.apply(lineOffset[0],lineOffsetOneSide(line, radius, units, resolution, false, true).geometry.coordinates);
    lineOffset[0].push.apply(lineOffset[0],arc(previousLinePoint, radius, previousLineBearing + 90, previousLineBearing - 90, units, resolution, true).geometry.coordinates);
    lineOffset[0].push.apply(lineOffset[0],lineOffsetOneSide(line, radius, units, resolution, true, true).geometry.coordinates);
    lineOffset[0].push.apply(lineOffset[0],arc(currentLinePoint, radius, nextLineBearing - 90, nextLineBearing + 90, units, resolution, true).geometry.coordinates);

    return offsetToBuffer(helpers.polygon(lineOffset));

  } else {

    lineOffset.push(ringOffsetOneSide(line, radius, units, resolution, false, true).geometry.coordinates);
    lineOffset.push(ringOffsetOneSide(line, radius, units, resolution, true, true).geometry.coordinates);

    return offsetToBuffer(helpers.polygon(lineOffset));
  }
}

function polygonBuffer(poly, radius, units, resolution) {
  var polygonOffset = [];

  poly = rewind(poly);

  poly.geometry.coordinates[0] = removeDuplicates(poly.geometry.coordinates[0]);
  for (var i = 1; i < poly.geometry.coordinates.length; i++) {
    poly.geometry.coordinates[i] = removeDuplicates(poly.geometry.coordinates[i]);
  }

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
  var coords = ring.geometry.coordinates; // ring is a lineString
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
  var angle = (Math.pow(-1, right + 1) * (bearing1 - bearing2)).modulo(360);
  var numSteps = Math.floor(angle/resMultiple);
  var step = numSteps; // Counting steps first is easier than checking angle (angle involves checking 'right', 'modulo(360)', lefthandedness of bearings
  var bearing = bearing1;
  // Add spoke for bearing1
  var spoke = destination(pt, radius, bearing1, units);
  arc.push(spoke.geometry.coordinates);
  // Add spokes for all bearings between bearing1 to bearing2
  // But don't add spokes if the angle is reflex and the shortcut preference is set. In that case, just add bearing1 and bearing2. This prevents double, zigzag-overlapping arcs, and potentially non-unique vertices, when a lineOffsetOneSide is run on both sides.
  if (!(angle > 180 && shortcut)) {
    while (step) {
      bearing = bearing + Math.pow(-1, !right + 1) * resMultiple;
      spoke = destination(pt, radius, bearing, units);
      arc.push(spoke.geometry.coordinates);
      step--;
    }
  }
  // Add spoke for bearing 2, but only if this spoke has not been added yet. Do this by checking the destination point, since slightly different bearings can create equal destination points.
  var spokeBearing2 = destination(pt, radius, bearing2, units);
  if (!equalArrays(spokeBearing2.geometry.coordinates,spoke.geometry.coordinates)) {
    arc.push(spokeBearing2.geometry.coordinates);
  }
  return helpers.lineString(arc)
}

function filterNetWinding(fc, filterFn) {
  var i = fc.features.length;
  while (i--) {
    if (!filterFn(fc.features[i].properties.netWinding)) {
        fc.features.splice(i, 1);
    } else {
      fc.features[i].properties = {};
    }
  }
  return fc;
}

function unionFeatureCollection(fc) {
  // Note: union takes a polygon, but return a polygon or multipolygon (which it can not take in). In case of buffes, however, it will always return a polygon
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

// This function awaits possible future use
function winding(poly){
  // compute winding of first ring
  var coords = poly.geometry.coordinates[0];
  var leftVtx = 0;
  for (var i = 0; i < coords.length-1; i++) { if (coords[i][0] < coords[leftVtx][0]) leftVtx = i; }
  return (coords[(leftVtx-1).modulo(coords.length-1)][1] > coords[(leftVtx+1).modulo(coords.length-1)][1]) ? 1 : -1;
}

// This function awaits possible future use
function rewind(poly){
  // outer ring to winding +1, inner rings to winding -1
  if (winding(helpers.polygon([poly.geometry.coordinates[0]])) == -1) poly.geometry.coordinates[0] = poly.geometry.coordinates[0].reverse();
  for (var i = 1; i < poly.geometry.coordinates.length; i++) {
    if (winding(helpers.polygon([poly.geometry.coordinates[i]])) == 1) poly.geometry.coordinates[i] = poly.geometry.coordinates[i].reverse();
  }
  return poly
}

function removeDuplicates(arr) {
  for (var i = arr.length-1; i > 0; i--) {
    if (equalArrays(arr[i],arr[i-1])) {
      arr.splice(i,1);
    }
  }
  return arr;
}


// Function to compare Arrays of numbers. From http://stackoverflow.com/questions/7837456/how-to-compare-arrays-in-javascript
function equalArrays(array1, array2) {
    // if the other array is a falsy value, return
    if (!array1 || !array2)
        return false;

    // compare lengths - can save a lot of time
    if (array1.length != array2.length)
        return false;

    for (var i = 0, l=array1.length; i < l; i++) {
        // Check if we have nested arrays
        if (array1[i] instanceof Array && array2[i] instanceof Array) {
            // recurse into the nested arrays
            if (!equalArrays(array1[i],array2[i]))
                return false;
        }
        else if (array1[i] != array2[i]) {
            // Warning - two different object instances will never be equal: {x:20} != {x:20}
            return false;
        }
    }
    return true;
}

// Fix Javascript modulo for negative number. From http://stackoverflow.com/questions/4467539/javascript-modulo-not-behaving
Number.prototype.modulo = function(n) {
    return ((this%n)+n)%n;
}
