'use strict';

var runBench = require('./bench-run.js');

var buffer_geodesicconvolution = require('../');
var buffer_master = require('turf-buffer');

var fs = require('fs');

var point = JSON.parse(fs.readFileSync(__dirname+'/fixtures/Point.geojson'));

runBench({
    'turf-buffer-point-master': function () {
        buffer_master(point, 1, 'kilometers');
    },
    'turf-buffer-point-geodesicconvolution': function () {
        buffer_geodesicconvolution(point, 1, 'kilometers');
    }
});
