'use strict';

var runBench = require('./bench-run.js');

var buffer_geodesicconvolution = require('../');
var buffer_master = require('turf-buffer');

var fs = require('fs');

var lineString = JSON.parse(fs.readFileSync(__dirname+'/fixtures/LineString-brazil.geojson'));

runBench({
    'turf-buffer-lineString-master': function () {
        buffer_master(lineString, 1, 'kilometers');
    },
    'turf-buffer-lineString-geodesicconvolution': function () {
        buffer_geodesicconvolution(lineString, 1, 'kilometers');
    }
});

runBench({
    'turf-buffer-lineString-master': function () {
        buffer_master(lineString, 100, 'kilometers');
    },
    'turf-buffer-lineString-geodesicconvolution': function () {
        buffer_geodesicconvolution(lineString, 100, 'kilometers');
    }
});
