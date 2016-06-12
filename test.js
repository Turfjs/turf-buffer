var test = require('tape');
var buffer = require('./');
var fs = require('fs');
var helpers = require('turf-helpers')

test('buffer', function(t){
  var pt = JSON.parse(fs.readFileSync(__dirname+'/fixtures/point.geojson'));
  var multipt = JSON.parse(fs.readFileSync(__dirname+'/fixtures/multiPoint.geojson'));
  var line = JSON.parse(fs.readFileSync(__dirname+'/fixtures/lineString.geojson'));
  var multiline = JSON.parse(fs.readFileSync(__dirname+'/fixtures/multiLineString.geojson'));
  var polygon = JSON.parse(fs.readFileSync(__dirname+'/fixtures/polygon.geojson'));
  var multipolygon = JSON.parse(fs.readFileSync(__dirname+'/fixtures/multiPolygon.geojson'));

  var buffPt = buffer(pt, 1, 'kilometers', 100);
  fs.writeFileSync(__dirname+'/fixtures/out/point.geojson', JSON.stringify(helpers.featureCollection([buffPt, pt])));

  var buffMultiPt = buffer(multipt, 1, 'kilometers');
  buffMultiPt.features.push(multipt)
  fs.writeFileSync(__dirname+'/fixtures/out/multipoint.geojson', JSON.stringify(buffMultiPt));

  var buffLine = buffer(line, 1, 'kilometers');
  fs.writeFileSync(__dirname+'/fixtures/out/line.geojson', JSON.stringify(helpers.featureCollection([buffLine, line])));

  var buffMultiLine = buffer(multiline, 1, 'kilometers');
  buffMultiLine.features.push(multiline)
  fs.writeFileSync(__dirname+'/fixtures/out/multiline.geojson', JSON.stringify(buffMultiLine));

  var buffPolygon = buffer(polygon, 1, 'kilometers');
  fs.writeFileSync(__dirname+'/fixtures/out/polygon.geojson', JSON.stringify(helpers.featureCollection([buffPolygon, polygon])));

  var buffMultiPolygon = buffer(multipolygon, 0.2, 'kilometers');
  buffMultiPolygon.features.push(multipolygon)
  fs.writeFileSync(__dirname+'/fixtures/out/multipolygon.geojson', JSON.stringify(buffMultiPolygon));


  t.ok(buffPt, 'should buffer a point');
  t.ok(buffMultiPt, 'should buffer a multipoint');
  t.ok(buffLine, 'should buffer a line');
  t.ok(buffMultiLine, 'should buffer a multiline');
  t.ok(buffPolygon, 'should buffer a polygon');
  t.ok(buffMultiPolygon, 'should buffer a multipolygon');

  t.end();
});
