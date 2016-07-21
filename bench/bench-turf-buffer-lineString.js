'use strict';

var runBench = require('./bench-run.js');

var buffer_geodesicconvolution = require('../');
var buffer_master = require('turf-buffer');

var fs = require('fs');

var lineString = JSON.parse(fs.readFileSync(__dirname+'/fixtures/LineString.geojson'));

runBench({
    'turf-buffer-lineString-master': function () {
        buffer_master(lineString, 1, 'kilometers');
    },
    'turf-buffer-lineString-geodesicconvolution': function () {
        buffer_geodesicconvolution(lineString, 1, 'kilometers');
    }
});
