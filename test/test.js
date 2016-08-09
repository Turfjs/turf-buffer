var test = require('tape');
var buffer = require('../');
var fs = require('fs');
var helpers = require('turf-helpers')

test('buffer', function(t){
  var point = JSON.parse(fs.readFileSync(__dirname+'/fixtures/Point.geojson'));
  var multiPoint = JSON.parse(fs.readFileSync(__dirname+'/fixtures/MultiPoint.geojson'));
  var lineString = JSON.parse(fs.readFileSync(__dirname+'/fixtures/LineString.geojson'));
  var multiLineString = JSON.parse(fs.readFileSync(__dirname+'/fixtures/MultiLineString.geojson'));
  var polygon = JSON.parse(fs.readFileSync(__dirname+'/fixtures/Polygon.geojson'));
  var multiPolygon = JSON.parse(fs.readFileSync(__dirname+'/fixtures/MultiPolygon.geojson'));

  // Tests for basic features
  var bufferedPoint = buffer(point, 1, 'kilometers', 100);
  fs.writeFileSync(__dirname+'/fixtures/out/Point.geojson', JSON.stringify(helpers.featureCollection([bufferedPoint, point])));

  var bufferedMultiPoint = buffer(multiPoint, 1, 'kilometers');
  bufferedMultiPoint.features.push(multiPoint)
  fs.writeFileSync(__dirname+'/fixtures/out/MultiPoint.geojson', JSON.stringify(bufferedMultiPoint));

  var bufferedLineString = buffer(lineString, 1, 'kilometers');
  fs.writeFileSync(__dirname+'/fixtures/out/LineString.geojson', JSON.stringify(helpers.featureCollection([bufferedLineString, lineString])));

  var bufferedMultiLineString = buffer(multiLineString, 1, 'kilometers');
  bufferedMultiLineString.features.push(multiLineString)
  fs.writeFileSync(__dirname+'/fixtures/out/MultiLineString.geojson', JSON.stringify(bufferedMultiLineString));

  var bufferedPolygon = buffer(polygon, 1, 'kilometers');
  fs.writeFileSync(__dirname+'/fixtures/out/Polygon.geojson', JSON.stringify(helpers.featureCollection([bufferedPolygon, polygon])));

  var bufferedMultiPolygon = buffer(multiPolygon, 0.2, 'kilometers');
  bufferedMultiPolygon.features.push(multiPolygon)
  fs.writeFileSync(__dirname+'/fixtures/out/MultiPolygon.geojson', JSON.stringify(bufferedMultiPolygon));

  // Test to buffer a circle LineString with buffer radius = circle radius
  var disk = buffer(point, 1, 'kilometers');
  var circle = disk;
  circle.geometry.type = "LineString";
  circle.geometry.coordinates = circle.geometry.coordinates[0];
  var bufferedCircle = buffer(circle, 1, 'kilometers');
  fs.writeFileSync(__dirname+'/fixtures/out/Circle.geojson', JSON.stringify(helpers.featureCollection([bufferedCircle, circle])));
  // the simplepolygon output has 32 rings, each with winding > 0

  t.ok(bufferedPoint, 'should buffer a Point');
  t.ok(bufferedMultiPoint, 'should buffer a MultiPoint');
  t.ok(bufferedLineString, 'should buffer a LineString');
  t.ok(bufferedMultiLineString, 'should buffer a MultiLineString');
  t.ok(bufferedPolygon, 'should buffer a Polygon');
  t.ok(bufferedMultiPolygon, 'should buffer a multiPolygon');
  t.ok(bufferedMultiPolygon, 'should buffer a circle with buffer radius = circle radius, without creating a hole or error in the middle');

  t.end();
});
