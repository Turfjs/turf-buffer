var test = require('tape');
var buffer = require('./');
var fs = require('fs');
var fc = require('turf-featurecollection')

test('buffer', function(t){
  var pt = JSON.parse(fs.readFileSync(__dirname+'/fixtures/Point.geojson'));
  var multipt = JSON.parse(fs.readFileSync(__dirname+'/fixtures/MultiPoint.geojson'));
  var line = JSON.parse(fs.readFileSync(__dirname+'/fixtures/LineString.geojson'));
  var polygon = JSON.parse(fs.readFileSync(__dirname+'/fixtures/Polygon.geojson'));

  var buffPt = buffer(pt, 10, 'miles', 1000);
  fs.writeFileSync(__dirname+'/fixtures/out/point.geojson', JSON.stringify(fc([buffPt, pt])));
  var buffMultiPt = buffer(multipt, 10, 'miles', 100);
  buffMultiPt.features.push(multipt)
  fs.writeFileSync(__dirname+'/fixtures/out/multipoint.geojson', JSON.stringify(buffMultiPt));
  var buffLine = buffer(line, 0.4, 'miles');
  buffLine.features.push(line)
  fs.writeFileSync(__dirname+'/fixtures/out/linestring.geojson', JSON.stringify(buffLine));
  //var buffPoly = buffer(pt, 10, 'miles');
  //var buffFC = buffer(fc, 10, 'miles'); 

  t.ok(buffPt, 'should buffer a point');
  t.ok(buffMultiPt, 'should buffer a multipoint');
  t.ok(buffLine, 'should buffer a line');
  //t.ok(buffPoly, 'should buffer a polygon');
  //t.ok(buffFC, 'should buffer featurecollection');

  t.end();
});