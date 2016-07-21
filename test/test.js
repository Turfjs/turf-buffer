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

  var BufferedPoint = buffer(point, 1, 'kilometers', 100);
  fs.writeFileSync(__dirname+'/fixtures/out/Point.geojson', JSON.stringify(helpers.featureCollection([BufferedPoint, point])));

  var BufferedMultiPoint = buffer(multiPoint, 1, 'kilometers');
  BufferedMultiPoint.features.push(multiPoint)
  fs.writeFileSync(__dirname+'/fixtures/out/MultiPoint.geojson', JSON.stringify(BufferedMultiPoint));

  var BufferedLineString = buffer(lineString, 1, 'kilometers');
  fs.writeFileSync(__dirname+'/fixtures/out/LineString.geojson', JSON.stringify(helpers.featureCollection([BufferedLineString, lineString])));

  var BufferedMultiLineString = buffer(multiLineString, 1, 'kilometers');
  BufferedMultiLineString.features.push(multiLineString)
  fs.writeFileSync(__dirname+'/fixtures/out/MultiLineString.geojson', JSON.stringify(BufferedMultiLineString));

  var bufferedPolygon = buffer(polygon, 1, 'kilometers');
  fs.writeFileSync(__dirname+'/fixtures/out/Polygon.geojson', JSON.stringify(helpers.featureCollection([bufferedPolygon, polygon])));

  var bufferedMultiPolygon = buffer(multiPolygon, 0.2, 'kilometers');
  bufferedMultiPolygon.features.push(multiPolygon)
  fs.writeFileSync(__dirname+'/fixtures/out/MultiPolygon.geojson', JSON.stringify(bufferedMultiPolygon));


  t.ok(BufferedPoint, 'should buffer a Point');
  t.ok(BufferedMultiPoint, 'should buffer a MultiPoint');
  t.ok(BufferedLineString, 'should buffer a LineString');
  t.ok(BufferedMultiLineString, 'should buffer a MultiLineString');
  t.ok(bufferedPolygon, 'should buffer a Polygon');
  t.ok(bufferedMultiPolygon, 'should buffer a multiPolygon');

  t.end();
});
