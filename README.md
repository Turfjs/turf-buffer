turf-buffer
===========
[![Build Status](https://travis-ci.org/Turfjs/turf-buffer.svg)](https://travis-ci.org/Turfjs/turf-buffer)

Buffers a point, linestring, or polygon Feature/FeatureCollection to a given radius. Units supported are miles, kilometers, and degrees.

###Install

```sh
npm install turf-buffer
```

###Parameters

|name|description|
|---|---|
|feature|A Feature or FeatureCollection of any type|
|distance|Distance to draw the buffer|
|unit|'miles' or 'kilometers'|

###Usage

```js
buffer(feature, distance, unit)
```

###Example

```js
var buffer = require('turf-buffer')
var point = require('turf-point')

var pt = point(14.616599, -90.548630)
var unit = 'miles'

var buffered = buffer(pt, 10, unit)

console.log(buffered)
``` 