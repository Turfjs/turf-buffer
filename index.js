var featurecollection = require('turf-featurecollection');
var destination = require('turf-destination');
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
    return lineBuffer();
  } else if(geom.type === 'MultiLineString') {

  } else if(geom.type === 'Polygon') { 

  } else if(geom.type === 'MultiPolygon') {

  }
}

function pointBuffer (pt, radius, units, resolution) {
  var ring = []
  var resMultiple = 360/resolution;
  for(var i  = 0; i < resolution; i++) {
    var spoke = destination(pt, radius, i*resMultiple, units)
    ring.push(spoke.geometry.coordinates)
  }
  if((ring[0][0] !== ring[ring.length-1][0]) && (ring[0][1] != ring[ring.length-1][1])) {
    ring.push([ring[0][0], ring[0][1]]);
  }
  
  return polygon([ring])
}

function lineBuffer () {

}