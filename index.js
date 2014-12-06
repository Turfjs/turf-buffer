var featurecollection = require('turf-featurecollection');
var destination = require('turf-destination');
var bearing = require('turf-bearing');
var point = require('turf-point');
var polygon = require('turf-polygon');

module.exports = function(feature, radius, units, resolution){
  if(!resolution) resolution = 36;
  var geom = feature.geometry
  if(geom.type === 'Point') {
    return pointBuffer(feature, radius, units, resolution);
  } else if(geom.type === 'MultiPoint') {
    var buffers = []
    geom.coordinates.forEach(function(coords) {
      buffers.push(pointBuffer(point(coords[0], coords[1]), radius, units, resolution));      
    });
    return featurecollection(buffers)
  } else if(geom.type === 'LineString') {
    return lineBuffer(feature, radius, units, resolution);
  } else if(geom.type === 'MultiLineString') {

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
  return polygon([ring])
}

function lineBuffer (line, radius, units, resolution) {
  //break line into segments
  var segments = [];
  for(var i = 0; i < line.geometry.coordinates.length-1; i++) {
    segments.push([line.geometry.coordinates[i], line.geometry.coordinates[i+1]]);
  }
  /*create a set of boxes parallel to the segments
  
    ---------

 ((|¯¯¯¯¯¯¯¯¯|))
(((|---------|)))
 ((|_________|))

  */
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
    console.log(JSON.stringify(poly))
  }
}