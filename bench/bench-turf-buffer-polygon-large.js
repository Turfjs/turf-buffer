'use strict';

var runBench = require('./bench-run.js');

var buffer_geodesicconvolution = require('../');
var buffer_master = require('turf-buffer');

var fs = require('fs');

var polygon = JSON.parse(fs.readFileSync(__dirname+'/fixtures/Polygon-brazil.geojson'));

runBench({
    'turf-buffer-polygon-master': function () {
        buffer_master(polygon, 1, 'kilometers');
    },
    'turf-buffer-polygon-geodesicconvolution': function () {
        buffer_geodesicconvolution(polygon, 1, 'kilometers');
    }
});

runBench({
    'turf-buffer-polygon-master': function () {
        buffer_master(polygon, 100, 'kilometers');
    },
    'turf-buffer-polygon-geodesicconvolution': function () {
        buffer_geodesicconvolution(polygon, 100, 'kilometers');
    }
});
