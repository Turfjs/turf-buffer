// http://stackoverflow.com/questions/839899/how-do-i-calculate-a-point-on-a-circles-circumference
// radians = degrees * (pi/180)
// https://github.com/bjornharrtell/jsts/blob/master/examples/buffer.html

var featurecollection = require('turf-featurecollection')
var polygon = require('turf-polygon')
var combine = require('turf-combine')

module.exports = function(feature, radius, units, done){
  var buffered;

  done = done || function () {};

  switch(units){
    case 'miles':
      radius = radius / 69.047
      break
    case 'kilometers':
      radius = radius / 111.12
      break
    case 'degrees':
      break
  }

  if(feature.type === 'FeatureCollection'){
    var multi = combine(feature);
    multi.properties = {}

    buffered = bufferOp(multi, radius);

    done(null, buffered);
    return buffered;
  }
  else{
    buffered = bufferOp(feature, radius);
    
    done(null, buffered);
    return buffered;
  }
}

var bufferOp = function(feature, radius){
  var reader = new jsts.io.GeoJSONReader()
  var geom = reader.read(JSON.stringify(feature.geometry))
  var buffered = geom.buffer(radius);
  var parser = new jsts.io.GeoJSONParser()
  buffered = parser.write(buffered)

  if(buffered.type === 'MultiPolygon'){
    buffered = {
      type: 'Feature',
      geometry: buffered,
      properties: {}
    }
    buffered = t.featurecollection([buffered])
  }
  else{
    buffered = t.featurecollection([t.polygon(buffered.coordinates)])
  }

  return buffered;
}


// Extracted from JSTS:
// todo: break this into a module

// setup jsts objects
var jsts = {}
jsts.io = {}
jsts.geom = {}
jsts.operation = {}
jsts.operation.buffer = {}
jsts.noding = {}
jsts.noding.snapround = {}
jsts.algorithm = {}
jsts.geomgraph = {}
var javascript = {}
javascript.util = {}



//    JSTS GeoJSONParser

/**
     * Create a new parser for GeoJSON
     *
     * @param {GeometryFactory}
     *          geometryFactory
     * @return An instance of GeoJsonParser.
     */
    jsts.io.GeoJSONParser = function(geometryFactory) {
        this.geometryFactory = geometryFactory || new jsts.geom.GeometryFactory();
        this.geometryTypes = ['Point', 'MultiPoint', 'LineString', 'MultiLineString', 'Polygon', 'MultiPolygon'];
    };

    /**
     * Deserialize a GeoJSON object and return the Geometry or Feature(Collection) with JSTS Geometries
     *
     * @param {}
     *          A GeoJSON object.
     * @return {} A Geometry instance or object representing a Feature(Collection) with Geometry instances.
     */
    jsts.io.GeoJSONParser.prototype.read = function(json) {
        var obj;
        if (typeof json === 'string') {
            obj = JSON.parse(json);
        } else {
            obj = json;
        }

        var type = obj.type;

        if (!this.parse[type]) {
            throw new Error('Unknown GeoJSON type: ' + obj.type);
        }

        if (this.geometryTypes.indexOf(type) != -1) {
            return this.parse[type].apply(this, [obj.coordinates]);
        } else if (type === 'GeometryCollection') {
            return this.parse[type].apply(this, [obj.geometries]);
        }

        // feature or feature collection
        return this.parse[type].apply(this, [obj]);
    };

    jsts.io.GeoJSONParser.prototype.parse = {
        /**
         * Parse a GeoJSON Feature object
         *
         * @param {Object}
         *          obj Object to parse.
         *
         * @return {Object} Feature with geometry/bbox converted to JSTS Geometries.
         */
        'Feature': function(obj) {
            var feature = {};

            // copy features
            for (var key in obj) {
                feature[key] = obj[key];
            }

            // parse geometry
            if (obj.geometry) {
                var type = obj.geometry.type;
                if (!this.parse[type]) {
                    throw new Error('Unknown GeoJSON type: ' + obj.type);
                }
                feature.geometry = this.read(obj.geometry);
            }

            // bbox
            if (obj.bbox) {
                feature.bbox = this.parse.bbox.apply(this, [obj.bbox]);
            }

            return feature;
        },

        /**
         * Parse a GeoJSON FeatureCollection object
         *
         * @param {Object}
         *          obj Object to parse.
         *
         * @return {Object} FeatureCollection with geometry/bbox converted to JSTS Geometries.
         */
        'FeatureCollection': function(obj) {
            var featureCollection = {};

            if (obj.features) {
                featureCollection.features = [];

                for (var i = 0; i < obj.features.length; ++i) {
                    featureCollection.features.push(this.read(obj.features[i]));
                }
            }

            if (obj.bbox) {
                featureCollection.bbox = this.parse.bbox.apply(this, [obj.bbox]);
            }

            return featureCollection;
        },


        /**
         * Convert the ordinates in an array to an array of jsts.geom.Coordinates
         *
         * @param {Array}
         *          array Array with {Number}s.
         *
         * @return {Array} Array with jsts.geom.Coordinates.
         */
        'coordinates': function(array) {
            var coordinates = [];

            for (var i = 0; i < array.length; ++i) {
                var sub = array[i];
                coordinates.push(new jsts.geom.Coordinate(sub[0], sub[1]));
            }

            return coordinates;
        },

        /**
         * Convert the bbox to a jsts.geom.LinearRing
         *
         * @param {Array}
         *          array Array with [xMin, yMin, xMax, yMax].
         *
         * @return {Array} Array with jsts.geom.Coordinates.
         */
        'bbox': function(array) {
            return this.geometryFactory.createLinearRing([
                new jsts.geom.Coordinate(array[0], array[1]),
                new jsts.geom.Coordinate(array[2], array[1]),
                new jsts.geom.Coordinate(array[2], array[3]),
                new jsts.geom.Coordinate(array[0], array[3]),
                new jsts.geom.Coordinate(array[0], array[1])
            ]);
        },


        /**
         * Convert an Array with ordinates to a jsts.geom.Point
         *
         * @param {Array}
         *          array Array with ordinates.
         *
         * @return {jsts.geom.Point} Point.
         */
        'Point': function(array) {
            var coordinate = new jsts.geom.Coordinate(array[0], array[1]);
            return this.geometryFactory.createPoint(coordinate);
        },

        /**
         * Convert an Array with coordinates to a jsts.geom.MultiPoint
         *
         * @param {Array}
         *          array Array with coordinates.
         *
         * @return {jsts.geom.MultiPoint} MultiPoint.
         */
        'MultiPoint': function(array) {
            var points = [];

            for (var i = 0; i < array.length; ++i) {
                points.push(this.parse.Point.apply(this, [array[i]]));
            }

            return this.geometryFactory.createMultiPoint(points);
        },

        /**
         * Convert an Array with coordinates to a jsts.geom.LineString
         *
         * @param {Array}
         *          array Array with coordinates.
         *
         * @return {jsts.geom.LineString} LineString.
         */
        'LineString': function(array) {
            var coordinates = this.parse.coordinates.apply(this, [array]);
            return this.geometryFactory.createLineString(coordinates);
        },

        /**
         * Convert an Array with coordinates to a jsts.geom.MultiLineString
         *
         * @param {Array}
         *          array Array with coordinates.
         *
         * @return {jsts.geom.MultiLineString} MultiLineString.
         */
        'MultiLineString': function(array) {
            var lineStrings = [];

            for (var i = 0; i < array.length; ++i) {
                lineStrings.push(this.parse.LineString.apply(this, [array[i]]));
            }

            return this.geometryFactory.createMultiLineString(lineStrings);
        },

        /**
         * Convert an Array to a jsts.geom.Polygon
         *
         * @param {Array}
         *          array Array with shell and holes.
         *
         * @return {jsts.geom.Polygon} Polygon.
         */
        'Polygon': function(array) {
            // shell
            var shellCoordinates = this.parse.coordinates.apply(this, [array[0]]);
            var shell = this.geometryFactory.createLinearRing(shellCoordinates);

            // holes
            var holes = [];
            for (var i = 1; i < array.length; ++i) {
                var hole = array[i];
                var coordinates = this.parse.coordinates.apply(this, [hole]);
                var linearRing = this.geometryFactory.createLinearRing(coordinates);
                holes.push(linearRing);
            }

            return this.geometryFactory.createPolygon(shell, holes);
        },

        /**
         * Convert an Array to a jsts.geom.MultiPolygon
         *
         * @param {Array}
         *          array Array of arrays with shell and rings.
         *
         * @return {jsts.geom.MultiPolygon} MultiPolygon.
         */
        'MultiPolygon': function(array) {
            var polygons = [];

            for (var i = 0; i < array.length; ++i) {
                var polygon = array[i];
                polygons.push(this.parse.Polygon.apply(this, [polygon]));
            }

            return this.geometryFactory.createMultiPolygon(polygons);
        },

        /**
         * Convert an Array to a jsts.geom.GeometryCollection
         *
         * @param {Array}
         *          array Array of GeoJSON geometries.
         *
         * @return {jsts.geom.GeometryCollection} GeometryCollection.
         */
        'GeometryCollection': function(array) {
            var geometries = [];

            for (var i = 0; i < array.length; ++i) {
                var geometry = array[i];
                geometries.push(this.read(geometry));
            }

            return this.geometryFactory.createGeometryCollection(geometries);
        }
    };

    /**
     * Serialize a Geometry object into GeoJSON
     *
     * @param {jsts.geom.geometry}
     *          geometry A Geometry or array of Geometries.
     * @return {Object} A GeoJSON object represting the input Geometry/Geometries.
     */
    jsts.io.GeoJSONParser.prototype.write = function(geometry) {
        var type = geometry.CLASS_NAME.slice(10);

        if (!this.extract[type]) {
            throw new Error('Geometry is not supported');
        }

        return this.extract[type].apply(this, [geometry]);
    };

    jsts.io.GeoJSONParser.prototype.extract = {
        /**
         * Convert a jsts.geom.Coordinate to an Array
         *
         * @param {jsts.geom.Coordinate}
         *          coordinate Coordinate to convert.
         *
         * @return {Array} Array of ordinates.
         */
        'coordinate': function(coordinate) {
            return [coordinate.x, coordinate.y];
        },

        /**
         * Convert a jsts.geom.Point to a GeoJSON object
         *
         * @param {jsts.geom.Point}
         *          point Point to convert.
         *
         * @return {Array} Array of 2 ordinates (paired to a coordinate).
         */
        'Point': function(point) {
            var array = this.extract.coordinate.apply(this, [point.coordinate]);

            return {
                type: 'Point',
                coordinates: array
            };
        },

        /**
         * Convert a jsts.geom.MultiPoint to a GeoJSON object
         *
         * @param {jsts.geom.MultiPoint}
         *          multipoint MultiPoint to convert.
         *
         * @return {Array} Array of coordinates.
         */
        'MultiPoint': function(multipoint) {
            var array = [];

            for (var i = 0; i < multipoint.geometries.length; ++i) {
                var point = multipoint.geometries[i];
                var geoJson = this.extract.Point.apply(this, [point]);
                array.push(geoJson.coordinates);
            }

            return {
                type: 'MultiPoint',
                coordinates: array
            };
        },

        /**
         * Convert a jsts.geom.LineString to a GeoJSON object
         *
         * @param {jsts.geom.LineString}
         *          linestring LineString to convert.
         *
         * @return {Array} Array of coordinates.
         */
        'LineString': function(linestring) {
            var array = [];

            for (var i = 0; i < linestring.points.length; ++i) {
                var coordinate = linestring.points[i];
                array.push(this.extract.coordinate.apply(this, [coordinate]));
            }

            return {
                type: 'LineString',
                coordinates: array
            };
        },

        /**
         * Convert a jsts.geom.MultiLineString to a GeoJSON object
         *
         * @param {jsts.geom.MultiLineString}
         *          multilinestring MultiLineString to convert.
         *
         * @return {Array} Array of Array of coordinates.
         */
        'MultiLineString': function(multilinestring) {
            var array = [];

            for (var i = 0; i < multilinestring.geometries.length; ++i) {
                var linestring = multilinestring.geometries[i];
                var geoJson = this.extract.LineString.apply(this, [linestring]);
                array.push(geoJson.coordinates);
            }

            return {
                type: 'MultiLineString',
                coordinates: array
            };
        },

        /**
         * Convert a jsts.geom.Polygon to a GeoJSON object
         *
         * @param {jsts.geom.Polygon}
         *          polygon Polygon to convert.
         *
         * @return {Array} Array with shell, holes.
         */
        'Polygon': function(polygon) {
            var array = [];

            // shell
            var shellGeoJson = this.extract.LineString.apply(this, [polygon.shell]);
            array.push(shellGeoJson.coordinates);

            // holes
            for (var i = 0; i < polygon.holes.length; ++i) {
                var hole = polygon.holes[i];
                var holeGeoJson = this.extract.LineString.apply(this, [hole]);
                array.push(holeGeoJson.coordinates);
            }

            return {
                type: 'Polygon',
                coordinates: array
            };
        },

        /**
         * Convert a jsts.geom.MultiPolygon to a GeoJSON object
         *
         * @param {jsts.geom.MultiPolygon}
         *          multipolygon MultiPolygon to convert.
         *
         * @return {Array} Array of polygons.
         */
        'MultiPolygon': function(multipolygon) {
            var array = [];

            for (var i = 0; i < multipolygon.geometries.length; ++i) {
                var polygon = multipolygon.geometries[i];
                var geoJson = this.extract.Polygon.apply(this, [polygon]);
                array.push(geoJson.coordinates);
            }

            return {
                type: 'MultiPolygon',
                coordinates: array
            };
        },

        /**
         * Convert a jsts.geom.GeometryCollection to a GeoJSON object
         *
         * @param {jsts.geom.GeometryCollection}
         *          collection GeometryCollection to convert.
         *
         * @return {Array} Array of geometries.
         */
        'GeometryCollection': function(collection) {
            var array = [];

            for (var i = 0; i < collection.geometries.length; ++i) {
                var geometry = collection.geometries[i];
                var type = geometry.CLASS_NAME.slice(10);
                array.push(this.extract[type].apply(this, [geometry]));
            }

            return {
                type: 'GeometryCollection',
                geometries: array
            };
        }
    };





//    JSTS GeoJSONReader

jsts.io.GeoJSONReader = function(geometryFactory) {
      this.geometryFactory = geometryFactory || new jsts.geom.GeometryFactory();
      this.precisionModel = this.geometryFactory.getPrecisionModel();
      this.parser = new jsts.io.GeoJSONParser(this.geometryFactory);
    };

    /**
     * Reads a GeoJSON representation of a {@link Geometry}
     *
     * @param {object}
     *          geoJson a GeoJSON Object or String.
     * @return {jsts.geom.Geometry} a <code>Geometry.</code>
     */
    jsts.io.GeoJSONReader.prototype.read = function(geoJson) {
      var geometry = this.parser.read(geoJson);

      if (this.precisionModel.getType() === jsts.geom.PrecisionModel.FIXED) {
        this.reducePrecision(geometry);
      }

      return geometry;
    };

    // NOTE: this is a hack
    jsts.io.GeoJSONReader.prototype.reducePrecision = function(geometry) {
      var i, len;

      if (geometry.coordinate) {
        this.precisionModel.makePrecise(geometry.coordinate);
      } else if (geometry.points) {
        for (i = 0, len = geometry.points.length; i < len; i++) {
          this.precisionModel.makePrecise(geometry.points[i]);
        }
      } else if (geometry.geometries) {
        for (i = 0, len = geometry.geometries.length; i < len; i++) {
          this.reducePrecision(geometry.geometries[i]);
        }
      }
    };




//    JSTS GeometryFactory

jsts.geom.GeometryFactory = function(precisionModel) {
  this.precisionModel = precisionModel || new jsts.geom.PrecisionModel();
};

jsts.geom.GeometryFactory.prototype.precisionModel = null;

jsts.geom.GeometryFactory.prototype.getPrecisionModel = function() {
  return this.precisionModel;
};


/**
 * Creates a Point using the given Coordinate; a null Coordinate will create an
 * empty Geometry.
 *
 * @param {Coordinate}
 *          coordinate Coordinate to base this Point on.
 * @return {Point} A new Point.
 */
jsts.geom.GeometryFactory.prototype.createPoint = function(coordinate) {
  var point = new jsts.geom.Point(coordinate, this);

  return point;
};


/**
 * Creates a LineString using the given Coordinates; a null or empty array will
 * create an empty LineString. Consecutive points must not be equal.
 *
 * @param {Coordinate[]}
 *          coordinates an array without null elements, or an empty array, or
 *          null.
 * @return {LineString} A new LineString.
 */
jsts.geom.GeometryFactory.prototype.createLineString = function(coordinates) {
  var lineString = new jsts.geom.LineString(coordinates, this);

  return lineString;
};


/**
 * Creates a LinearRing using the given Coordinates; a null or empty array will
 * create an empty LinearRing. The points must form a closed and simple
 * linestring. Consecutive points must not be equal.
 *
 * @param {Coordinate[]}
 *          coordinates an array without null elements, or an empty array, or
 *          null.
 * @return {LinearRing} A new LinearRing.
 */
jsts.geom.GeometryFactory.prototype.createLinearRing = function(coordinates) {
  var linearRing = new jsts.geom.LinearRing(coordinates, this);

  return linearRing;
};


/**
 * Constructs a <code>Polygon</code> with the given exterior boundary and
 * interior boundaries.
 *
 * @param {LinearRing}
 *          shell the outer boundary of the new <code>Polygon</code>, or
 *          <code>null</code> or an empty <code>LinearRing</code> if the
 *          empty geometry is to be created.
 * @param {LinearRing[]}
 *          holes the inner boundaries of the new <code>Polygon</code>, or
 *          <code>null</code> or empty <code>LinearRing</code> s if the
 *          empty geometry is to be created.
 * @return {Polygon} A new Polygon.
 */
jsts.geom.GeometryFactory.prototype.createPolygon = function(shell, holes) {
  var polygon = new jsts.geom.Polygon(shell, holes, this);

  return polygon;
};


jsts.geom.GeometryFactory.prototype.createMultiPoint = function(points) {
  if (points && points[0] instanceof jsts.geom.Coordinate) {
    var converted = [];
    var i;
    for (i = 0; i < points.length; i++) {
      converted.push(this.createPoint(points[i]));
    }
    points = converted;
  }

  return new jsts.geom.MultiPoint(points, this);
};

jsts.geom.GeometryFactory.prototype.createMultiLineString = function(
    lineStrings) {
  return new jsts.geom.MultiLineString(lineStrings, this);
};

jsts.geom.GeometryFactory.prototype.createMultiPolygon = function(polygons) {
  return new jsts.geom.MultiPolygon(polygons, this);
};


/**
 * Build an appropriate <code>Geometry</code>, <code>MultiGeometry</code>,
 * or <code>GeometryCollection</code> to contain the <code>Geometry</code>s
 * in it. For example:<br>
 *
 * <ul>
 * <li> If <code>geomList</code> contains a single <code>Polygon</code>,
 * the <code>Polygon</code> is returned.
 * <li> If <code>geomList</code> contains several <code>Polygon</code>s, a
 * <code>MultiPolygon</code> is returned.
 * <li> If <code>geomList</code> contains some <code>Polygon</code>s and
 * some <code>LineString</code>s, a <code>GeometryCollection</code> is
 * returned.
 * <li> If <code>geomList</code> is empty, an empty
 * <code>GeometryCollection</code> is returned
 * </ul>
 *
 * Note that this method does not "flatten" Geometries in the input, and hence
 * if any MultiGeometries are contained in the input a GeometryCollection
 * containing them will be returned.
 *
 * @param geomList
 *          the <code>Geometry</code>s to combine.
 * @return {Geometry} a <code>Geometry</code> of the "smallest", "most
 *         type-specific" class that can contain the elements of
 *         <code>geomList</code> .
 */
jsts.geom.GeometryFactory.prototype.buildGeometry = function(geomList) {

  /**
   * Determine some facts about the geometries in the list
   */
  var geomClass = null;
  var isHeterogeneous = false;
  var hasGeometryCollection = false;
  for (var i = geomList.iterator(); i.hasNext();) {
    var geom = i.next();

    var partClass = geom.CLASS_NAME;

    if (geomClass === null) {
      geomClass = partClass;
    }
    if (!(partClass === geomClass)) {
      isHeterogeneous = true;
    }
    if (geom.isGeometryCollectionBase())
      hasGeometryCollection = true;
  }

  /**
   * Now construct an appropriate geometry to return
   */
  // for the empty geometry, return an empty GeometryCollection
  if (geomClass === null) {
    return this.createGeometryCollection(null);
  }
  if (isHeterogeneous || hasGeometryCollection) {
    return this.createGeometryCollection(geomList.toArray());
  }
  // at this point we know the collection is hetereogenous.
  // Determine the type of the result from the first Geometry in the list
  // this should always return a geometry, since otherwise an empty collection
  // would have already been returned
  var geom0 = geomList.get(0);
  var isCollection = geomList.size() > 1;
  if (isCollection) {
    if (geom0 instanceof jsts.geom.Polygon) {
      return this.createMultiPolygon(geomList.toArray());
    } else if (geom0 instanceof jsts.geom.LineString) {
      return this.createMultiLineString(geomList.toArray());
    } else if (geom0 instanceof jsts.geom.Point) {
      return this.createMultiPoint(geomList.toArray());
    }
    jsts.util.Assert.shouldNeverReachHere('Unhandled class: ' + geom0);
  }
  return geom0;
};

jsts.geom.GeometryFactory.prototype.createGeometryCollection = function(
    geometries) {
  return new jsts.geom.GeometryCollection(geometries, this);
};

/**
 * Creates a {@link Geometry} with the same extent as the given envelope. The
 * Geometry returned is guaranteed to be valid. To provide this behaviour, the
 * following cases occur:
 * <p>
 * If the <code>Envelope</code> is:
 * <ul>
 * <li>null : returns an empty {@link Point}
 * <li>a point : returns a non-empty {@link Point}
 * <li>a line : returns a two-point {@link LineString}
 * <li>a rectangle : returns a {@link Polygon}> whose points are (minx, miny),
 * (minx, maxy), (maxx, maxy), (maxx, miny), (minx, miny).
 * </ul>
 *
 * @param {jsts.geom.Envelope}
 *          envelope the <code>Envelope</code> to convert.
 * @return {jsts.geom.Geometry} an empty <code>Point</code> (for null
 *         <code>Envelope</code>s), a <code>Point</code> (when min x = max
 *         x and min y = max y) or a <code>Polygon</code> (in all other cases).
 */
jsts.geom.GeometryFactory.prototype.toGeometry = function(envelope) {
  // null envelope - return empty point geometry
  if (envelope.isNull()) {
    return this.createPoint(null);
  }

  // point?
  if (envelope.getMinX() === envelope.getMaxX() &&
      envelope.getMinY() === envelope.getMaxY()) {
    return this.createPoint(new jsts.geom.Coordinate(envelope.getMinX(),
        envelope.getMinY()));
  }

  // vertical or horizontal line?
  if (envelope.getMinX() === envelope.getMaxX() ||
      envelope.getMinY() === envelope.getMaxY()) {
    return this.createLineString([
        new jsts.geom.Coordinate(envelope.getMinX(), envelope.getMinY()),
        new jsts.geom.Coordinate(envelope.getMaxX(), envelope.getMaxY())]);
  }

  // create a CW ring for the polygon
  return this.createPolygon(this.createLinearRing([
      new jsts.geom.Coordinate(envelope.getMinX(), envelope.getMinY()),
      new jsts.geom.Coordinate(envelope.getMinX(), envelope.getMaxY()),
      new jsts.geom.Coordinate(envelope.getMaxX(), envelope.getMaxY()),
      new jsts.geom.Coordinate(envelope.getMaxX(), envelope.getMinY()),
      new jsts.geom.Coordinate(envelope.getMinX(), envelope.getMinY())]), null);
};



//    JSTS PrecisionModel

jsts.geom.PrecisionModel = function(modelType) {
  if (typeof modelType === 'number') {
    this.modelType = jsts.geom.PrecisionModel.FIXED;
    this.scale = modelType;
    return;
  }

  this.modelType = modelType || jsts.geom.PrecisionModel.FLOATING;

  if (this.modelType === jsts.geom.PrecisionModel.FIXED) {
    this.scale = 1.0;
  }
};


/**
 * @type {string}
 */
jsts.geom.PrecisionModel.FLOATING = 'FLOATING';


/**
 * @type {string}
 */
jsts.geom.PrecisionModel.FIXED = 'FIXED';


/**
 * @type {string}
 */
jsts.geom.PrecisionModel.FLOATING_SINGLE = 'FLOATING_SINGLE';

jsts.geom.PrecisionModel.prototype.scale = null;
jsts.geom.PrecisionModel.prototype.modelType = null;


/**
 * Tests whether the precision model supports floating point
 *
 * @return {boolean} if the precision model supports floating point.
 */
jsts.geom.PrecisionModel.prototype.isFloating = function() {
  return this.modelType === jsts.geom.PrecisionModel.FLOATING ||
      this.modelType === jsts.geom.PrecisionModel.FLOATING_SINLGE;
};

/**
 * Returns the scale factor used to specify a fixed precision model. The number
 * of decimal places of precision is equal to the base-10 logarithm of the scale
 * factor. Non-integral and negative scale factors are supported. Negative scale
 * factors indicate that the places of precision is to the left of the decimal
 * point.
 *
 * @return the scale factor for the fixed precision model.
 */
jsts.geom.PrecisionModel.prototype.getScale = function() {
  return this.scale;
};

/**
 * @return {string} always jsts.geom.PrecisionModel.FLOATING.
 */
jsts.geom.PrecisionModel.prototype.getType = function() {
  return this.modelType;
};

jsts.geom.PrecisionModel.prototype.equals = function(other) {
  return true;

  if (!(other instanceof jsts.geom.PrecisionModel)) {
    return false;
  }
  var otherPrecisionModel = other;
  return this.modelType === otherPrecisionModel.modelType &&
      this.scale === otherPrecisionModel.scale;
};


/**
 * Rounds a numeric value to the PrecisionModel grid. Asymmetric Arithmetic
 * Rounding is used, to provide uniform rounding behaviour no matter where the
 * number is on the number line.
 * <p>
 * This method has no effect on NaN values.
 * <p>
 * <b>Note:</b> Java's <code>Math#rint</code> uses the "Banker's Rounding"
 * algorithm, which is not suitable for precision operations elsewhere in JTS.
 */
jsts.geom.PrecisionModel.prototype.makePrecise = function(val) {
  if (val instanceof jsts.geom.Coordinate) {
    this.makePrecise2(val);
    return;
  }

  // don't change NaN values
  if (isNaN(val))
    return val;

  // TODO: support single precision?
  /*if (this.modelType == FLOATING_SINGLE) {
    float floatSingleVal = (float) val;
    return (double) floatSingleVal;
  }*/
  if (this.modelType === jsts.geom.PrecisionModel.FIXED) {
    return Math.round(val * this.scale) / this.scale;
  }
  // modelType == FLOATING - no rounding necessary
  return val;
};


/**
 * Rounds a Coordinate to the PrecisionModel grid.
 */
jsts.geom.PrecisionModel.prototype.makePrecise2 = function(coord) {
  // optimization for full precision
  if (this.modelType === jsts.geom.PrecisionModel.FLOATING)
    return;

  coord.x = this.makePrecise(coord.x);
  coord.y = this.makePrecise(coord.y);
  // MD says it's OK that we're not makePrecise'ing the z [Jon Aquino]
};


/**
 * Compares this {@link PrecisionModel} object with the specified object for
 * order. A PrecisionModel is greater than another if it provides greater
 * precision. The comparison is based on the value returned by the
 * {@link #getMaximumSignificantDigits} method. This comparison is not strictly
 * accurate when comparing floating precision models to fixed models; however,
 * it is correct when both models are either floating or fixed.
 *
 * @param o
 *          the <code>PrecisionModel</code> with which this
 *          <code>PrecisionModel</code> is being compared.
 * @return a negative integer, zero, or a positive integer as this
 *         <code>PrecisionModel</code> is less than, equal to, or greater than
 *         the specified <code>PrecisionModel.</code>
 */
jsts.geom.PrecisionModel.prototype.compareTo = function(o) {
  var other = o;

  // TODO: needs to be ported for fixed precision

  // var sigDigits = this.getMaximumSignificantDigits();
  // var otherSigDigits = other.getMaximumSignificantDigits();
  // return (new Integer(sigDigits)).compareTo(new Integer(otherSigDigits));

  return 0;
};



//    JSTS GEOMETRY

jsts.geom.Geometry = function(factory) {
  this.factory = factory;
};


/**
 * The bounding box of this <code>Geometry</code>.
 */
jsts.geom.Geometry.prototype.envelope = null;

/**
 * The {@link GeometryFactory} used to create this Geometry
 *
 * @protected
 */
jsts.geom.Geometry.prototype.factory = null;


/**
 * Returns the name of this object's <code>com.vivid.jts.geom</code>
 * interface.
 *
 * @return {string} the name of this <code>Geometry</code>s most specific
 *         <code>jsts.geom</code> interface.
 */
jsts.geom.Geometry.prototype.getGeometryType = function() {
  return 'Geometry';
};


/**
 * Returns true if the array contains any non-empty <code>Geometry</code>s.
 *
 * @param {Geometry[]}
 *          geometries an array of <code>Geometry</code>s; no elements may be
 *          <code>null.</code>
 * @return {boolean} <code>true</code> if any of the <code>Geometry</code>s
 *         <code>isEmpty</code> methods return <code>false.</code>
 */
jsts.geom.Geometry.hasNonEmptyElements = function(geometries) {
  var i;
  for (i = 0; i < geometries.length; i++) {
    if (!geometries[i].isEmpty()) {
      return true;
    }
  }
  return false;
};


/**
 * Returns true if the array contains any <code>null</code> elements.
 *
 * @param {Object[]}
 *          array an array to validate.
 * @return {boolean} <code>true</code> if any of <code>array</code>s
 *         elements are <code>null.</code>
 */
jsts.geom.Geometry.hasNullElements = function(array) {
  var i;
  for (i = 0; i < array.length; i++) {
    if (array[i] === null) {
      return true;
    }
  }
  return false;
};


/**
 * Gets the factory which contains the context in which this geometry was
 * created.
 *
 * @return {jsts.geom.GeometryFactory} the factory for this geometry.
 */
jsts.geom.Geometry.prototype.getFactory = function() {
  // NOTE: Geometry could be created without JSTS constructor so need to check
  // for member data
  // TODO: above should not happen
  if (this.factory === null || this.factory === undefined) {
    this.factory = new jsts.geom.GeometryFactory();
  }

  return this.factory;
};


/**
 * Returns the number of {@link Geometry}s in a {@link GeometryCollection} (or
 * 1, if the geometry is not a collection).
 *
 * @return {number} the number of geometries contained in this geometry.
 */
jsts.geom.Geometry.prototype.getNumGeometries = function() {
  return 1;
};


/**
 * Returns an element {@link Geometry} from a {@link GeometryCollection} (or
 * <code>this</code>, if the geometry is not a collection).
 *
 * @param {number}
 *          n the index of the geometry element.
 * @return {Geometry} the n'th geometry contained in this geometry.
 */
jsts.geom.Geometry.prototype.getGeometryN = function(n) {
  return this;
};


/**
 * Returns the <code>PrecisionModel</code> used by the <code>Geometry</code>.
 *
 * @return {PrecisionModel} the specification of the grid of allowable points,
 *         for this <code>Geometry</code> and all other <code>Geometry</code>s.
 */
jsts.geom.Geometry.prototype.getPrecisionModel = function() {
  return this.getFactory().getPrecisionModel();
};



/**
 * Returns a vertex of this <code>Geometry</code> (usually, but not
 * necessarily, the first one). The returned coordinate should not be assumed to
 * be an actual Coordinate object used in the internal representation.
 *
 * @return {Coordinate} a {@link Coordinate} which is a vertex of this
 *         <code>Geometry</code>. null if this Geometry is empty.
 */
jsts.geom.Geometry.prototype.getCoordinate = function() {
  throw new jsts.error.AbstractMethodInvocationError();
};


/**
 * Returns an array containing the values of all the vertices for this geometry.
 * If the geometry is a composite, the array will contain all the vertices for
 * the components, in the order in which the components occur in the geometry.
 * <p>
 * In general, the array cannot be assumed to be the actual internal storage for
 * the vertices. Thus modifying the array may not modify the geometry itself.
 * Use the {@link CoordinateSequence#setOrdinate} method (possibly on the
 * components) to modify the underlying data. If the coordinates are modified,
 * {@link #geometryChanged} must be called afterwards.
 *
 * @return {Coordinate[]} the vertices of this <code>Geometry.</code>
 * @see geometryChanged
 * @see CoordinateSequence#setOrdinate
 */
jsts.geom.Geometry.prototype.getCoordinates = function() {
  throw new jsts.error.AbstractMethodInvocationError();
};


/**
 * Returns the count of this <code>Geometry</code>s vertices. The
 * <code>Geometry</code> s contained by composite <code>Geometry</code>s
 * must be Geometry's; that is, they must implement <code>getNumPoints</code>
 *
 * @return {number} the number of vertices in this <code>Geometry.</code>
 */
jsts.geom.Geometry.prototype.getNumPoints = function() {
  throw new jsts.error.AbstractMethodInvocationError();
};


/**
 * Tests whether this {@link Geometry} is simple. In general, the SFS
 * specification of simplicity follows the rule:
 * <UL>
 * <LI> A Geometry is simple iff the only self-intersections are at boundary
 * points.
 * </UL>
 * Simplicity is defined for each {@link Geometry} subclass as follows:
 * <ul>
 * <li>Valid polygonal geometries are simple by definition, so
 * <code>isSimple</code> trivially returns true.
 * <li>Linear geometries are simple iff they do not self-intersect at points
 * other than boundary points.
 * <li>Zero-dimensional geometries (points) are simple iff they have no
 * repeated points.
 * <li>Empty <code>Geometry</code>s are always simple
 * <ul>
 *
 * @return {boolean} <code>true</code> if this <code>Geometry</code> has any
 *         points of self-tangency, self-intersection or other anomalous points.
 * @see #isValid
 */
jsts.geom.Geometry.prototype.isSimple = function() {
  this.checkNotGeometryCollection(this);
  var op = new jsts.operation.IsSimpleOp(this);
  return op.isSimple();
};


/**
 * Tests the validity of this <code>Geometry</code>. Subclasses provide their
 * own definition of "valid".
 *
 * @return {boolean} <code>true</code> if this <code>Geometry</code> is
 *         valid.
 *
 * @see IsValidOp
 */
jsts.geom.Geometry.prototype.isValid = function() {
  var isValidOp = new jsts.operation.valid.IsValidOp(this);
  return isValidOp.isValid();
};


/**
 * Returns whether or not the set of points in this <code>Geometry</code> is
 * empty.
 *
 * @return {boolean} <code>true</code> if this <code>Geometry</code> equals
 *         the empty geometry.
 */
jsts.geom.Geometry.prototype.isEmpty = function() {
  throw new jsts.error.AbstractMethodInvocationError();
};


/**
 * Returns the minimum distance between this <code>Geometry</code> and the
 * <code>Geometry</code> g
 *
 * @param {Geometry}
 *          g the <code>Geometry</code> from which to compute the distance.
 * @return {number} the distance between the geometries. 0 if either input
 *         geometry is empty.
 * @throws IllegalArgumentException
 *           if g is null
 */
jsts.geom.Geometry.prototype.distance = function(g) {
  return jsts.operation.distance.DistanceOp.distance(this, g);
};


/**
 * Tests whether the distance from this <code>Geometry</code> to another is
 * less than or equal to a specified value.
 *
 * @param {Geometry}
 *          geom the Geometry to check the distance to.
 * @param {number}
 *          distance the distance value to compare.
 * @return {boolean} <code>true</code> if the geometries are less than
 *         <code>distance</code> apart.
 */
jsts.geom.Geometry.prototype.isWithinDistance = function(geom, distance) {
  var envDist = this.getEnvelopeInternal().distance(geom.getEnvelopeInternal());
  if (envDist > distance) {
    return false;
  }
  return DistanceOp.isWithinDistance(this, geom, distance);
};

jsts.geom.Geometry.prototype.isRectangle = function() {
  // Polygon overrides to check for actual rectangle
  return false;
};

/**
 * Returns the area of this <code>Geometry</code>. Areal Geometries have a
 * non-zero area. They override this function to compute the area. Others return
 * 0.0
 *
 * @return the area of the Geometry.
 */
jsts.geom.Geometry.prototype.getArea = function() {
  return 0.0;
};

/**
 * Returns the length of this <code>Geometry</code>. Linear geometries return
 * their length. Areal geometries return their perimeter. They override this
 * function to compute the area. Others return 0.0
 *
 * @return the length of the Geometry.
 */
jsts.geom.Geometry.prototype.getLength = function() {
  return 0.0;
};

/**
 * Computes the centroid of this <code>Geometry</code>. The centroid is equal
 * to the centroid of the set of component Geometries of highest dimension
 * (since the lower-dimension geometries contribute zero "weight" to the
 * centroid)
 *
 * @return a {@link Point} which is the centroid of this Geometry.
 */
jsts.geom.Geometry.prototype.getCentroid = function() {
  if (this.isEmpty()) {
    return null;
  }
  var cent;
  var centPt = null;
  var dim = this.getDimension();
  if (dim === 0) {
    cent = new jsts.algorithm.CentroidPoint();
    cent.add(this);
    centPt = cent.getCentroid();
  } else if (dim === 1) {
    cent = new jsts.algorithm.CentroidLine();
    cent.add(this);
    centPt = cent.getCentroid();
  } else {
    cent = new jsts.algorithm.CentroidArea();
    cent.add(this);
    centPt = cent.getCentroid();
  }
  return this.createPointFromInternalCoord(centPt, this);

};


/**
 * Computes an interior point of this <code>Geometry</code>. An interior
 * point is guaranteed to lie in the interior of the Geometry, if it possible to
 * calculate such a point exactly. Otherwise, the point may lie on the boundary
 * of the geometry.
 *
 * @return {Point} a {@link Point} which is in the interior of this Geometry.
 */
jsts.geom.Geometry.prototype.getInteriorPoint = function() {
  var intPt;
  var interiorPt = null;
  var dim = this.getDimension();
  if (dim === 0) {
    intPt = new jsts.algorithm.InteriorPointPoint(this);
    interiorPt = intPt.getInteriorPoint();
  } else if (dim === 1) {
    intPt = new jsts.algorithm.InteriorPointLine(this);
    interiorPt = intPt.getInteriorPoint();
  } else {
    intPt = new jsts.algorithm.InteriorPointArea(this);
    interiorPt = intPt.getInteriorPoint();
  }
  return this.createPointFromInternalCoord(interiorPt, this);
};


/**
 * Returns the dimension of this geometry. The dimension of a geometry is is the
 * topological dimension of its embedding in the 2-D Euclidean plane. In the JTS
 * spatial model, dimension values are in the set {0,1,2}.
 * <p>
 * Note that this is a different concept to the dimension of the vertex
 * {@link Coordinate}s. The geometry dimension can never be greater than the
 * coordinate dimension. For example, a 0-dimensional geometry (e.g. a Point)
 * may have a coordinate dimension of 3 (X,Y,Z).
 *
 * @return {number} the topological dimension of this geometry.
 */
jsts.geom.Geometry.prototype.getDimension = function() {
  throw new jsts.error.AbstractMethodInvocationError();
};


/**
 * Returns the boundary, or an empty geometry of appropriate dimension if this
 * <code>Geometry</code> is empty. (In the case of zero-dimensional
 * geometries, ' an empty GeometryCollection is returned.) For a discussion of
 * this function, see the OpenGIS Simple Features Specification. As stated in
 * SFS Section 2.1.13.1, "the boundary of a Geometry is a set of Geometries of
 * the next lower dimension."
 *
 * @return {Geometry} the closure of the combinatorial boundary of this
 *         <code>Geometry.</code>
 */
jsts.geom.Geometry.prototype.getBoundary = function() {
  throw new jsts.error.AbstractMethodInvocationError();
};


/**
 * Returns the dimension of this <code>Geometry</code>s inherent boundary.
 *
 * @return {number} the dimension of the boundary of the class implementing this
 *         interface, whether or not this object is the empty geometry. Returns
 *         <code>Dimension.FALSE</code> if the boundary is the empty geometry.
 */
jsts.geom.Geometry.prototype.getBoundaryDimension = function() {
  throw new jsts.error.AbstractMethodInvocationError();
};


/**
 * Returns this <code>Geometry</code>s bounding box. If this
 * <code>Geometry</code> is the empty geometry, returns an empty
 * <code>Point</code>. If the <code>Geometry</code> is a point, returns a
 * non-empty <code>Point</code>. Otherwise, returns a <code>Polygon</code>
 * whose points are (minx, miny), (maxx, miny), (maxx, maxy), (minx, maxy),
 * (minx, miny).
 *
 * @return {Geometry} an empty <code>Point</code> (for empty
 *         <code>Geometry</code>s), a <code>Point</code> (for
 *         <code>Point</code>s) or a <code>Polygon</code> (in all other
 *         cases).
 */
jsts.geom.Geometry.prototype.getEnvelope = function() {
  return this.getFactory().toGeometry(this.getEnvelopeInternal());
};


/**
 * Returns the minimum and maximum x and y values in this <code>Geometry</code>,
 * or a null <code>Envelope</code> if this <code>Geometry</code> is empty.
 *
 * @return {Envelope} this <code>Geometry</code>s bounding box; if the
 *         <code>Geometry</code> is empty, <code>Envelope#isNull</code> will
 *         return <code>true.</code>
 */
jsts.geom.Geometry.prototype.getEnvelopeInternal = function() {
  if (this.envelope === null) {
    this.envelope = this.computeEnvelopeInternal();
  }
  return this.envelope;
};


/**
 * Tests whether this geometry is disjoint from the specified geometry.
 * <p>
 * The <code>disjoint</code> predicate has the following equivalent
 * definitions:
 * <ul>
 * <li>The two geometries have no point in common
 * <li>The DE-9IM Intersection Matrix for the two geometries matches
 * <code>[FF*FF****]</code>
 * <li><code>! g.intersects(this)</code> (<code>disjoint</code> is the
 * inverse of <code>intersects</code>)
 * </ul>
 *
 * @param {Geometry}
 *          g the <code>Geometry</code> with which to compare this
 *          <code>Geometry.</code>
 * @return {boolean} <code>true</code> if the two <code>Geometry</code>s
 *         are disjoint.
 *
 * @see Geometry#intersects
 */
jsts.geom.Geometry.prototype.disjoint = function(g) {
  return !this.intersects(g);
};


/**
 * Tests whether this geometry touches the specified geometry.
 * <p>
 * The <code>touches</code> predicate has the following equivalent
 * definitions:
 * <ul>
 * <li>The geometries have at least one point in common, but their interiors do
 * not intersect.
 * <li>The DE-9IM Intersection Matrix for the two geometries matches
 * <code>[FT*******]</code> or <code>[F**T*****]</code> or
 * <code>[F***T****]</code>
 * </ul>
 * If both geometries have dimension 0, this predicate returns
 * <code>false</code>
 *
 * @param {Geometry}
 *          g the <code>Geometry</code> with which to compare this
 *          <code>Geometry.</code>
 * @return {boolean} <code>true</code> if the two <code>Geometry</code>s
 *         touch; Returns <code>false</code> if both <code>Geometry</code>s
 *         are points.
 */
jsts.geom.Geometry.prototype.touches = function(g) {
  // short-circuit test
  if (!this.getEnvelopeInternal().intersects(g.getEnvelopeInternal())) {
    return false;
  }
  return this.relate(g).isTouches(this.getDimension(), g.getDimension());
};


/**
 * Tests whether this geometry intersects the specified geometry.
 * <p>
 * The <code>intersects</code> predicate has the following equivalent
 * definitions:
 * <ul>
 * <li>The two geometries have at least one point in common
 * <li>The DE-9IM Intersection Matrix for the two geometries matches
 * <code>[T********]</code> or <code>[*T*******]</code> or
 * <code>[***T*****]</code> or <code>[****T****]</code>
 * <li><code>! g.disjoint(this)</code> (<code>intersects</code> is the
 * inverse of <code>disjoint</code>)
 * </ul>
 *
 * @param {Geometry}
 *          g the <code>Geometry</code> with which to compare this
 *          <code>Geometry.</code>
 * @return {boolean} <code>true</code> if the two <code>Geometry</code>s
 *         intersect.
 *
 * @see Geometry#disjoint
 */
jsts.geom.Geometry.prototype.intersects = function(g) {

  // short-circuit envelope test
  if (!this.getEnvelopeInternal().intersects(g.getEnvelopeInternal())) {
    return false;
  }

  // optimization for rectangle arguments
  if (this.isRectangle()) {
    return RectangleIntersects.intersects(this, g);
  }
  if (g.isRectangle()) {
    return RectangleIntersects.intersects(g, this);
  }
  // general case
  return this.relate(g).isIntersects();
};


/**
 * Tests whether this geometry crosses the specified geometry.
 * <p>
 * The <code>crosses</code> predicate has the following equivalent
 * definitions:
 * <ul>
 * <li>The geometries have some but not all interior points in common.
 * <li>The DE-9IM Intersection Matrix for the two geometries matches
 * <ul>
 * <li><code>[T*T******]</code> (for P/L, P/A, and L/A situations)
 * <li><code>[T*****T**]</code> (for L/P, A/P, and A/L situations)
 * <li><code>[0********]</code> (for L/L situations)
 * </ul>
 * </ul>
 * For any other combination of dimensions this predicate returns
 * <code>false</code>.
 * <p>
 * The SFS defined this predicate only for P/L, P/A, L/L, and L/A situations.
 * JTS extends the definition to apply to L/P, A/P and A/L situations as well,
 * in order to make the relation symmetric.
 *
 * @param {Geometry}
 *          g the <code>Geometry</code> with which to compare this
 *          <code>Geometry.</code>
 * @return {boolean} <code>true</code> if the two <code>Geometry</code>s
 *         cross.
 */
jsts.geom.Geometry.prototype.crosses = function(g) {
  // short-circuit test
  if (!this.getEnvelopeInternal().intersects(g.getEnvelopeInternal())) {
    return false;
  }
  return this.relate(g).isCrosses(this.getDimension(), g.getDimension());
};


/**
 * Tests whether this geometry is within the specified geometry.
 * <p>
 * The <code>within</code> predicate has the following equivalent definitions:
 * <ul>
 * <li>Every point of this geometry is a point of the other geometry, and the
 * interiors of the two geometries have at least one point in common.
 * <li>The DE-9IM Intersection Matrix for the two geometries matches
 * <code>[T*F**F***]</code>
 * <li><code>g.contains(this)</code> (<code>within</code> is the converse
 * of <code>contains</code>)
 * </ul>
 * An implication of the definition is that "The boundary of a Geometry is not
 * within the Geometry". In other words, if a geometry A is a subset of the
 * points in the boundary of a geomtry B, <code>A.within(B) = false</code>
 *
 * @param {Geometry}
 *          g the <code>Geometry</code> with which to compare this
 *          <code>Geometry.</code>
 * @return {boolean} <code>true</code> if this <code>Geometry</code> is
 *         within <code>other.</code>
 *
 * @see Geometry#contains
 */
jsts.geom.Geometry.prototype.within = function(g) {
  return g.contains(this);
};


/**
 * Tests whether this geometry contains the specified geometry.
 * <p>
 * The <code>contains</code> predicate has the following equivalent
 * definitions:
 * <ul>
 * <li>Every point of the other geometry is a point of this geometry, and the
 * interiors of the two geometries have at least one point in common.
 * <li>The DE-9IM Intersection Matrix for the two geometries matches
 * <code>[T*****FF*]</code>
 * <li><code>g.within(this)</code> (<code>contains</code> is the converse
 * of <code>within</code>)
 * </ul>
 * An implication of the definition is that "Geometries do not contain their
 * boundary". In other words, if a geometry A is a subset of the points in the
 * boundary of a geometry B, <code>B.contains(A) = false</code>
 *
 * @param {Geometry}
 *          g the <code>Geometry</code> with which to compare this
 *          <code>Geometry.</code>
 * @return {boolean} <code>true</code> if this <code>Geometry</code>
 *         contains <code>g.</code>
 *
 * @see Geometry#within
 */
jsts.geom.Geometry.prototype.contains = function(g) {
  // short-circuit test
  if (!this.getEnvelopeInternal().contains(g.getEnvelopeInternal())) {
    return false;
  }
  // optimization for rectangle arguments
  if (this.isRectangle()) {
    return RectangleContains.contains(this, g);
  }
  // general case
  return this.relate(g).isContains();
};


/**
 * Tests whether this geometry overlaps the specified geometry.
 * <p>
 * The <code>overlaps</code> predicate has the following equivalent
 * definitions:
 * <ul>
 * <li>The geometries have at least one point each not shared by the other (or
 * equivalently neither covers the other), they have the same dimension, and the
 * intersection of the interiors of the two geometries has the same dimension as
 * the geometries themselves.
 * <li>The DE-9IM Intersection Matrix for the two geometries matches
 * <code>[T*T***T**]</code> (for two points or two surfaces) or
 * <code>[1*T***T**]</code> (for two curves)
 * </ul>
 * If the geometries are of different dimension this predicate returns
 * <code>false</code>.
 *
 * @param {Geometry}
 *          g the <code>Geometry</code> with which to compare this
 *          <code>Geometry.</code>
 * @return {boolean} <code>true</code> if the two <code>Geometry</code>s
 *         overlap.
 */
jsts.geom.Geometry.prototype.overlaps = function(g) {
  // short-circuit test
  if (!this.getEnvelopeInternal().intersects(g.getEnvelopeInternal())) {
    return false;
  }
  return this.relate(g).isOverlaps(this.getDimension(), g.getDimension());
};


/**
 * Tests whether this geometry covers the specified geometry.
 * <p>
 * The <code>covers</code> predicate has the following equivalent definitions:
 * <ul>
 * <li>Every point of the other geometry is a point of this geometry.
 * <li>The DE-9IM Intersection Matrix for the two geometries matches
 * <code>[T*****FF*]</code> or <code>[*T****FF*]</code> or
 * <code>[***T**FF*]</code> or <code>[****T*FF*]</code>
 * <li><code>g.coveredBy(this)</code> (<code>covers</code> is the converse
 * of <code>coveredBy</code>)
 * </ul>
 * If either geometry is empty, the value of this predicate is <tt>false</tt>.
 * <p>
 * This predicate is similar to {@link #contains}, but is more inclusive (i.e.
 * returns <tt>true</tt> for more cases). In particular, unlike
 * <code>contains</code> it does not distinguish between points in the
 * boundary and in the interior of geometries. For most situations,
 * <code>covers</code> should be used in preference to <code>contains</code>.
 * As an added benefit, <code>covers</code> is more amenable to optimization,
 * and hence should be more performant.
 *
 * @param {Geometry}
 *          g the <code>Geometry</code> with which to compare this
 *          <code>Geometry.</code>
 * @return {boolean} <code>true</code> if this <code>Geometry</code> covers
 *         <code>g.</code>
 *
 * @see Geometry#contains
 * @see Geometry#coveredBy
 */
jsts.geom.Geometry.prototype.covers = function(g) {
  // short-circuit test
  if (!this.getEnvelopeInternal().covers(g.getEnvelopeInternal())) {
    return false;
  }
  // optimization for rectangle arguments
  if (this.isRectangle()) {
    // since we have already tested that the test envelope is covered
    return true;
  }
  return this.relate(g).isCovers();
};


/**
 * Tests whether this geometry is covered by the specified geometry.
 * <p>
 * The <code>coveredBy</code> predicate has the following equivalent
 * definitions:
 * <ul>
 * <li>Every point of this geometry is a point of the other geometry.
 * <li>The DE-9IM Intersection Matrix for the two geometries matches
 * <code>[T*F**F***]</code> or <code>[*TF**F***]</code> or
 * <code>[**FT*F***]</code> or <code>[**F*TF***]</code>
 * <li><code>g.covers(this)</code> (<code>coveredBy</code> is the converse
 * of <code>covers</code>)
 * </ul>
 * If either geometry is empty, the value of this predicate is <tt>false</tt>.
 * <p>
 * This predicate is similar to {@link #within}, but is more inclusive (i.e.
 * returns <tt>true</tt> for more cases).
 *
 * @param {Geometry}
 *          g the <code>Geometry</code> with which to compare this
 *          <code>Geometry.</code>
 * @return {boolean} <code>true</code> if this <code>Geometry</code> is
 *         covered by <code>g.</code>
 *
 * @see Geometry#within
 * @see Geometry#covers
 */
jsts.geom.Geometry.prototype.coveredBy = function(g) {
  return g.covers(this);
};


/**
 * Tests whether the elements in the DE-9IM {@link IntersectionMatrix} for the
 * two <code>Geometry</code>s match the elements in
 * <code>intersectionPattern</code>. The pattern is a 9-character string,
 * with symbols drawn from the following set:
 * <UL>
 * <LI> 0 (dimension 0)
 * <LI> 1 (dimension 1)
 * <LI> 2 (dimension 2)
 * <LI> T ( matches 0, 1 or 2)
 * <LI> F ( matches FALSE)
 * <LI> * ( matches any value)
 * </UL>
 * For more information on the DE-9IM, see the <i>OpenGIS Simple Features
 * Specification</i>.
 *
 * @param {Geometry}
 *          other the <code>Geometry</code> with which to compare this
 *          <code>Geometry.</code>
 * @param {string}
 *          intersectionPattern the pattern against which to check the
 *          intersection matrix for the two <code>Geometry</code>s.
 * @return {boolean} <code>true</code> if the DE-9IM intersection matrix for
 *         the two <code>Geometry</code>s match
 *         <code>intersectionPattern.</code>
 * @see IntersectionMatrix
 */
jsts.geom.Geometry.prototype.relate = function(g, intersectionPattern) {
  if (arguments.length === 1) {
    return this.relate2.apply(this, arguments);
  }

  return this.relate2(g).matches(intersectionPattern);
};


/**
 * Returns the DE-9IM {@link IntersectionMatrix} for the two
 * <code>Geometry</code>s.
 *
 * @param {Geometry}
 *          g the <code>Geometry</code> with which to compare this
 *          <code>Geometry.</code>
 * @return {IntersectionMatrix} an {@link IntersectionMatrix} describing the
 *         intersections of the interiors, boundaries and exteriors of the two
 *         <code>Geometry</code>s.
 */
jsts.geom.Geometry.prototype.relate2 = function(g) {
  this.checkNotGeometryCollection(this);
  this.checkNotGeometryCollection(g);
  return jsts.operation.relate.RelateOp.relate(this, g);
};


/**
 * Tests whether this geometry is topologically equal to the argument geometry
 * as defined by the SFS <tt>equals</tt> predicate.
 * <p>
 * The SFS <code>equals</code> predicate has the following equivalent
 * definitions:
 * <ul>
 * <li>The two geometries have at least one point in common, and no point of
 * either geometry lies in the exterior of the other geometry.
 * <li>The DE-9IM Intersection Matrix for the two geometries matches the
 * pattern <tt>T*F**FFF*</tt>
 * <pre>
 * T*F
 * **F
 * FF*
 * </pre>
 *
 * </ul>
 * <b>Note</b> that this method computes <b>topologically equality</b>. For
 * structural equality, see {@link #equalsExact(Geometry)}.
 *
 * @param {Geometry}
 *          g the <code>Geometry</code> with which to compare this
 *          <code>Geometry.</code>
 * @return {boolean} <code>true</code> if the two <code>Geometry</code>s
 *         are topologically equal.
 *
 * @see #equalsExact(Geometry)
 */
jsts.geom.Geometry.prototype.equalsTopo = function(g) {
  // short-circuit test
  if (!this.getEnvelopeInternal().equals(g.getEnvelopeInternal())) {
    return false;
  }
  return this.relate(g).isEquals(this.getDimension(), g.getDimension());
};


/**
 * Tests whether this geometry is structurally and numerically equal to a given
 * <tt>Object</tt>. If the argument <tt>Object</tt> is not a
 * <tt>Geometry</tt>, the result is <tt>false</tt>. Otherwise, the result
 * is computed using {@link #equalsExact(Geometry)}.
 * <p>
 * This method is provided to fulfill the Java contract for value-based object
 * equality. In conjunction with {@link #hashCode()} it provides semantics which
 * are most useful for using <tt>Geometry</tt>s as keys and values in Java
 * collections.
 * <p>
 * Note that to produce the expected result the input geometries should be in
 * normal form. It is the caller's responsibility to perform this where required
 * (using {@link Geometry#norm() or {@link #normalize()} as appropriate).
 *
 * @param {Object}
 *          o the Object to compare.
 * @return {boolean} true if this geometry is exactly equal to the argument.
 *
 * @see #equalsExact(Geometry)
 * @see #hashCode()
 * @see #norm()
 * @see #normalize()
 */
jsts.geom.Geometry.prototype.equals = function(o) {
  if (o instanceof jsts.geom.Geometry || o instanceof jsts.geom.LinearRing ||
      o instanceof jsts.geom.Polygon ||
      o instanceof jsts.geom.GeometryCollection ||
      o instanceof jsts.geom.MultiPoint ||
      o instanceof jsts.geom.MultiLineString ||
      o instanceof jsts.geom.MultiPolygon) {
    return this.equalsExact(o);
  }
  return false;
};

/**
 * Computes a buffer area around this geometry having the given width and with a
 * specified accuracy of approximation for circular arcs, and using a specified
 * end cap style.
 * <p>
 * Mathematically-exact buffer area boundaries can contain circular arcs. To
 * represent these arcs using linear geometry they must be approximated with
 * line segments. The <code>quadrantSegments</code> argument allows
 * controlling the accuracy of the approximation by specifying the number of
 * line segments used to represent a quadrant of a circle
 * <p>
 * The end cap style specifies the buffer geometry that will be created at the
 * ends of linestrings. The styles provided are:
 * <ul>
 * <li><tt>BufferOp.CAP_ROUND</tt> - (default) a semi-circle
 * <li><tt>BufferOp.CAP_BUTT</tt> - a straight line perpendicular to the end
 * segment
 * <li><tt>BufferOp.CAP_SQUARE</tt> - a half-square
 * </ul>
 * <p>
 * The buffer operation always returns a polygonal result. The negative or
 * zero-distance buffer of lines and points is always an empty {@link Polygon}.
 * This is also the result for the buffers of degenerate (zero-area) polygons.
 *
 * @param {number}
 *          distance the width of the buffer (may be positive, negative or 0).
 * @param {number}
 *          quadrantSegments the number of line segments used to represent a
 *          quadrant of a circle.
 * @param {number}
 *          endCapStyle the end cap style to use.
 * @return {Geometry} a polygonal geometry representing the buffer region (which
 *         may be empty).
 *
 * @throws TopologyException
 *           if a robustness error occurs
 *
 * @see #buffer(double)
 * @see #buffer(double, int)
 * @see BufferOp
 */
jsts.geom.Geometry.prototype.buffer = function(distance, quadrantSegments, endCapStyle) {
  var params = new jsts.operation.buffer.BufferParameters(quadrantSegments, endCapStyle)
  return jsts.operation.buffer.BufferOp.bufferOp2(this, distance, params);
};


/**
 * Computes the smallest convex <code>Polygon</code> that contains all the
 * points in the <code>Geometry</code>. This obviously applies only to
 * <code>Geometry</code> s which contain 3 or more points; the results for
 * degenerate cases are specified as follows: <TABLE>
 * <TR>
 * <TH> Number of <code>Point</code>s in argument <code>Geometry</code>
 * </TH>
 * <TH> <code>Geometry</code> class of result </TH>
 * </TR>
 * <TR>
 * <TD> 0 </TD>
 * <TD> empty <code>GeometryCollection</code> </TD>
 * </TR>
 * <TR>
 * <TD> 1 </TD>
 * <TD> <code>Point</code> </TD>
 * </TR>
 * <TR>
 * <TD> 2 </TD>
 * <TD> <code>LineString</code> </TD>
 * </TR>
 * <TR>
 * <TD> 3 or more </TD>
 * <TD> <code>Polygon</code> </TD>
 * </TR>
 * </TABLE>
 *
 * @return {Geometry} the minimum-area convex polygon containing this
 *         <code>Geometry</code>' s points.
 */
jsts.geom.Geometry.prototype.convexHull = function() {
  return new jsts.algorithm.ConvexHull(this).getConvexHull();
};


/**
 * Computes a <code>Geometry</code> representing the points shared by this
 * <code>Geometry</code> and <code>other</code>. {@link GeometryCollection}s
 * support intersection with homogeneous collection types, with the semantics
 * that the result is a {@link GeometryCollection} of the intersection of each
 * element of the target with the argument.
 *
 * @param {Geometry}
 *          other the <code>Geometry</code> with which to compute the
 *          intersection.
 * @return {Geometry} the points common to the two <code>Geometry</code>s.
 * @throws TopologyException
 *           if a robustness error occurs
 * @throws IllegalArgumentException
 *           if the argument is a non-empty GeometryCollection
 */
jsts.geom.Geometry.prototype.intersection = function(other) {
  /**
   * TODO: MD - add optimization for P-A case using Point-In-Polygon
   */
  // special case: if one input is empty ==> empty
  if (this.isEmpty()) {
    return this.getFactory().createGeometryCollection(null);
  }
  if (other.isEmpty()) {
    return this.getFactory().createGeometryCollection(null);
  }

  // compute for GCs
  if (this.isGeometryCollection(this)) {
    var g2 = other;
    // TODO: probably not straightforward to port...
    /*
     * return GeometryCollectionMapper.map(this, new
     * GeometryCollectionMapper.MapOp() { public Geometry map(Geometry g) {
     * return g.intersection(g2); } });
     */
  }

  this.checkNotGeometryCollection(this);
  this.checkNotGeometryCollection(other);
  return jsts.operation.overlay.snap.SnapIfNeededOverlayOp.overlayOp(this,
      other, jsts.operation.overlay.OverlayOp.INTERSECTION);
};


/**
 * Computes a <code>Geometry</code> representing all the points in this
 * <code>Geometry</code> and <code>other</code>.
 *
 * Or without arguments:
 *
 * Computes the union of all the elements of this geometry. Heterogeneous
 * {@link GeometryCollection}s are fully supported.
 *
 * The result obeys the following contract:
 * <ul>
 * <li>Unioning a set of {@link LineString}s has the effect of fully noding
 * and dissolving the linework.
 * <li>Unioning a set of {@link Polygon}s will always return a
 * {@link Polygonal} geometry (unlike {link #union(Geometry)}, which may return
 * geometrys of lower dimension if a topology collapse occurred.
 * </ul>
 *
 * @param {Geometry}
 *          other the <code>Geometry</code> with which to compute the union.
 * @return {Geometry} a set combining the points of this <code>Geometry</code>
 *         and the points of <code>other.</code>
 * @throws TopologyException
 *           if a robustness error occurs
 * @throws IllegalArgumentException
 *           if either input is a non-empty GeometryCollection
 */
jsts.geom.Geometry.prototype.union = function(other) {
  if (arguments.length === 0) {
    return jsts.operation.union.UnaryUnionOp.union(this);
  }

  // special case: if either input is empty ==> other input
  if (this.isEmpty()) {
    return other.clone();
  }
  if (other.isEmpty()) {
    return this.clone();
  }

  // TODO: optimize if envelopes of geometries do not intersect

  this.checkNotGeometryCollection(this);
  this.checkNotGeometryCollection(other);
  return jsts.operation.overlay.snap.SnapIfNeededOverlayOp.overlayOp(this,
      other, jsts.operation.overlay.OverlayOp.UNION);
};


/**
 * Computes a <code>Geometry</code> representing the points making up this
 * <code>Geometry</code> that do not make up <code>other</code>. This
 * method returns the closure of the resultant <code>Geometry</code>.
 *
 * @param {Geometry}
 *          other the <code>Geometry</code> with which to compute the
 *          difference.
 * @return {Geometry} the point set difference of this <code>Geometry</code>
 *         with <code>other.</code>
 * @throws TopologyException
 *           if a robustness error occurs
 * @throws IllegalArgumentException
 *           if either input is a non-empty GeometryCollection
 */
jsts.geom.Geometry.prototype.difference = function(other) {
  // mod to handle empty cases better - return type of input
  // if (this.isEmpty() || other.isEmpty()) return (Geometry) clone();

  // special case: if A.isEmpty ==> empty; if B.isEmpty ==> A
  if (this.isEmpty()) {
    return this.getFactory().createGeometryCollection(null);
  }
  if (other.isEmpty()) {
    return this.clone();
  }

  this.checkNotGeometryCollection(this);
  this.checkNotGeometryCollection(other);
  return jsts.operation.overlay.snap.SnapIfNeededOverlayOp.overlayOp(this,
      other, jsts.operation.overlay.OverlayOp.DIFFERENCE);
};


/**
 * Returns a set combining the points in this <code>Geometry</code> not in
 * <code>other</code>, and the points in <code>other</code> not in this
 * <code>Geometry</code>. This method returns the closure of the resultant
 * <code>Geometry</code>.
 *
 * @param {Geometry}
 *          other the <code>Geometry</code> with which to compute the
 *          symmetric difference.
 * @return {Geometry} the point set symmetric difference of this
 *         <code>Geometry</code> with <code>other.</code>
 * @throws TopologyException
 *           if a robustness error occurs
 * @throws IllegalArgumentException
 *           if either input is a non-empty GeometryCollection
 */
jsts.geom.Geometry.prototype.symDifference = function(other) {
  // special case: if either input is empty ==> other input
  if (this.isEmpty()) {
    return other.clone();
  }
  if (other.isEmpty()) {
    return this.clone();
  }

  this.checkNotGeometryCollection(this);
  this.checkNotGeometryCollection(other);
  return jsts.operation.overlay.snap.SnapIfNeededOverlayOp.overlayOp(this,
      other, jsts.operation.overlay.OverlayOp.SYMDIFFERENCE);
};

/**
 * Returns true if the two <code>Geometry</code>s are exactly equal, up to a
 * specified distance tolerance. Two Geometries are exactly equal within a
 * distance tolerance if and only if:
 * <ul>
 * <li>they have the same class
 * <li>they have the same values for their vertices, within the given tolerance
 * distance, in exactly the same order.
 * </ul>
 * If this and the other <code>Geometry</code>s are composites and any
 * children are not <code>Geometry</code>s, returns <code>false</code>.
 *
 * @param {Geometry}
 *          other the <code>Geometry</code> with which to compare this
 *          <code>Geometry.</code>
 * @param {number}
 *          tolerance distance at or below which two <code>Coordinate</code>s
 *          are considered equal.
 * @return {boolean}
 */
jsts.geom.Geometry.prototype.equalsExact = function(other, tolerance) {
  throw new jsts.error.AbstractMethodInvocationError();
};

/**
 * Tests whether two geometries are exactly equal in their normalized forms.
 * This is a convenience method which creates normalized versions of both
 * geometries before computing {@link #equalsExact(Geometry)}. This method is
 * relatively expensive to compute. For maximum performance, the client should
 * instead perform normalization itself at an appropriate point during
 * execution.
 *
 * @param {Geometry}
 *          g a Geometry.
 * @return {boolean} true if the input geometries are exactly equal in their
 *         normalized form.
 */
jsts.geom.Geometry.prototype.equalsNorm = function(g) {
  if (g === null || g === undefined)
    return false;
  return this.norm().equalsExact(g.norm());
};


/**
 * Performs an operation with or on this <code>Geometry</code> and its
 * subelement <code>Geometry</code>s (if any). Only GeometryCollections and
 * subclasses have subelement Geometry's.
 *
 * @param filter
 *          the filter to apply to this <code>Geometry</code> (and its
 *          children, if it is a <code>GeometryCollection</code>).
 */
jsts.geom.Geometry.prototype.apply = function(filter) {
  throw new jsts.error.AbstractMethodInvocationError();
};

/**
 * Creates and returns a full copy of this {@link Geometry} object (including
 * all coordinates contained by it). Subclasses are responsible for overriding
 * this method and copying their internal data. Overrides should call this
 * method first.
 *
 * @return a clone of this instance.
 */
jsts.geom.Geometry.prototype.clone = function() {
  throw new jsts.error.AbstractMethodInvocationError();
};


/**
 * Converts this <code>Geometry</code> to <b>normal form</b> (or <b>
 * canonical form</b> ). Normal form is a unique representation for
 * <code>Geometry</code> s. It can be used to test whether two
 * <code>Geometry</code>s are equal in a way that is independent of the
 * ordering of the coordinates within them. Normal form equality is a stronger
 * condition than topological equality, but weaker than pointwise equality. The
 * definitions for normal form use the standard lexicographical ordering for
 * coordinates. "Sorted in order of coordinates" means the obvious extension of
 * this ordering to sequences of coordinates.
 */
jsts.geom.Geometry.prototype.normalize = function() {
  throw new jsts.error.AbstractMethodInvocationError();
};

/**
 * Creates a new Geometry which is a normalized copy of this Geometry.
 *
 * @return a normalized copy of this geometry.
 * @see #normalize()
 */
jsts.geom.Geometry.prototype.norm = function() {
  var copy = this.clone();
  copy.normalize();
  return copy;
};


/**
 * Returns whether this <code>Geometry</code> is greater than, equal to, or
 * less than another <code>Geometry</code>.
 * <P>
 *
 * If their classes are different, they are compared using the following
 * ordering:
 * <UL>
 * <LI> Point (lowest)
 * <LI> MultiPoint
 * <LI> LineString
 * <LI> LinearRing
 * <LI> MultiLineString
 * <LI> Polygon
 * <LI> MultiPolygon
 * <LI> GeometryCollection (highest)
 * </UL>
 * If the two <code>Geometry</code>s have the same class, their first
 * elements are compared. If those are the same, the second elements are
 * compared, etc.
 *
 * @param {Geometry}
 *          other a <code>Geometry</code> with which to compare this
 *          <code>Geometry.</code>
 * @return {number} a positive number, 0, or a negative number, depending on
 *         whether this object is greater than, equal to, or less than
 *         <code>o</code>, as defined in "Normal Form For Geometry" in the
 *         JTS Technical Specifications.
 */
jsts.geom.Geometry.prototype.compareTo = function(o) {
  var other = o;
  if (this.getClassSortIndex() !== other.getClassSortIndex()) {
    return this.getClassSortIndex() - other.getClassSortIndex();
  }
  if (this.isEmpty() && other.isEmpty()) {
    return 0;
  }
  if (this.isEmpty()) {
    return -1;
  }
  if (other.isEmpty()) {
    return 1;
  }
  return this.compareToSameClass(o);
};

/**
 * Returns whether the two <code>Geometry</code>s are equal, from the point
 * of view of the <code>equalsExact</code> method. Called by
 * <code>equalsExact</code> . In general, two <code>Geometry</code> classes
 * are considered to be "equivalent" only if they are the same class. An
 * exception is <code>LineString</code> , which is considered to be equivalent
 * to its subclasses.
 *
 * @param {Geometry}
 *          other the <code>Geometry</code> with which to compare this
 *          <code>Geometry</code> for equality.
 * @return {boolean} <code>true</code> if the classes of the two
 *         <code>Geometry</code> s are considered to be equal by the
 *         <code>equalsExact</code> method.
 */
jsts.geom.Geometry.prototype.isEquivalentClass = function(other) {
  if (this instanceof jsts.geom.Point && other instanceof jsts.geom.Point) {
    return true;
  } else if (this instanceof jsts.geom.LineString &&
      (other instanceof jsts.geom.LineString | other instanceof jsts.geom.LinearRing)) {
    return true;
  } else if (this instanceof jsts.geom.LinearRing &&
      (other instanceof jsts.geom.LineString | other instanceof jsts.geom.LinearRing)) {
    return true;
  } else if (this instanceof jsts.geom.Polygon &&
      (other instanceof jsts.geom.Polygon)) {
    return true;
  } else if (this instanceof jsts.geom.MultiPoint &&
      (other instanceof jsts.geom.MultiPoint)) {
    return true;
  } else if (this instanceof jsts.geom.MultiLineString &&
      (other instanceof jsts.geom.MultiLineString)) {
    return true;
  } else if (this instanceof jsts.geom.MultiPolygon &&
      (other instanceof jsts.geom.MultiPolygon)) {
    return true;
  } else if (this instanceof jsts.geom.GeometryCollection &&
      (other instanceof jsts.geom.GeometryCollection)) {
    return true;
  }

  return false;
};



/**
 * Throws an exception if <code>g</code>'s class is
 * <code>GeometryCollection</code> . (Its subclasses do not trigger an
 * exception).
 *
 * @param {Geometry}
 *          g the <code>Geometry</code> to check.
 * @throws Error
 *           if <code>g</code> is a <code>GeometryCollection</code> but not
 *           one of its subclasses
 */
jsts.geom.Geometry.prototype.checkNotGeometryCollection = function(g) {
  if (g.isGeometryCollectionBase()) {
    throw new jsts.error.IllegalArgumentError(
        'This method does not support GeometryCollection');
  }
};


/**
 *
 * @return {boolean} true if this is a GeometryCollection.
 */
jsts.geom.Geometry.prototype.isGeometryCollection = function() {
  return (this instanceof jsts.geom.GeometryCollection);
};

/**
 *
 * @return {boolean} true if this is a GeometryCollection but not subclass.
 */
jsts.geom.Geometry.prototype.isGeometryCollectionBase = function() {
  return (this.CLASS_NAME === 'jsts.geom.GeometryCollection');
};


/**
 * Returns the minimum and maximum x and y values in this <code>Geometry</code>,
 * or a null <code>Envelope</code> if this <code>Geometry</code> is empty.
 * Unlike <code>getEnvelopeInternal</code>, this method calculates the
 * <code>Envelope</code> each time it is called;
 * <code>getEnvelopeInternal</code> caches the result of this method.
 *
 * @return {Envelope} this <code>Geometry</code>s bounding box; if the
 *         <code>Geometry</code> is empty, <code>Envelope#isNull</code> will
 *         return <code>true.</code>
 */
jsts.geom.Geometry.prototype.computeEnvelopeInternal = function() {
  throw new jsts.error.AbstractMethodInvocationError();
};


/**
 * Returns whether this <code>Geometry</code> is greater than, equal to, or
 * less than another <code>Geometry</code> having the same class.
 *
 * @param o
 *          a <code>Geometry</code> having the same class as this
 *          <code>Geometry.</code>
 * @return a positive number, 0, or a negative number, depending on whether this
 *         object is greater than, equal to, or less than <code>o</code>, as
 *         defined in "Normal Form For Geometry" in the JTS Technical
 *         Specifications.
 */
jsts.geom.Geometry.prototype.compareToSameClass = function(o) {
  throw new jsts.error.AbstractMethodInvocationError();
};

/**
 * Returns the first non-zero result of <code>compareTo</code> encountered as
 * the two <code>Collection</code>s are iterated over. If, by the time one of
 * the iterations is complete, no non-zero result has been encountered, returns
 * 0 if the other iteration is also complete. If <code>b</code> completes
 * before <code>a</code>, a positive number is returned; if a before b, a
 * negative number.
 *
 * @param {Array}
 *          a a <code>Collection</code> of <code>Comparable</code>s.
 * @param {Array}
 *          b a <code>Collection</code> of <code>Comparable</code>s.
 * @return {number} the first non-zero <code>compareTo</code> result, if any;
 *         otherwise, zero.
 */
jsts.geom.Geometry.prototype.compare = function(a, b) {
  var i = a.iterator();
  var j = b.iterator();
  while (i.hasNext() && j.hasNext()) {
    var aElement = i.next();
    var bElement = j.next();
    var comparison = aElement.compareTo(bElement);
    if (comparison !== 0) {
      return comparison;
    }
  }
  if (i.hasNext()) {
    return 1;
  }
  if (j.hasNext()) {
    return -1;
  }
  return 0;
};


/**
 * @param {jsts.geom.Coordinate}
 *          a first Coordinate to compare.
 * @param {jsts.geom.Coordinate}
 *          b second Coordinate to compare.
 * @param {number}
 *          tolerance tolerance when comparing.
 * @return {boolean} true if equal.
 */
jsts.geom.Geometry.prototype.equal = function(a, b, tolerance) {
  if (tolerance === undefined || tolerance === null || tolerance === 0) {
    return a.equals(b);
  }
  return a.distance(b) <= tolerance;
};

/**
 * @private
 */
jsts.geom.Geometry.prototype.getClassSortIndex = function() {
  var sortedClasses = [jsts.geom.Point, jsts.geom.MultiPoint,
      jsts.geom.LineString, jsts.geom.LinearRing, jsts.geom.MultiLineString,
      jsts.geom.Polygon, jsts.geom.MultiPolygon, jsts.geom.GeometryCollection];

  for (var i = 0; i < sortedClasses.length; i++) {
    if (this instanceof sortedClasses[i])
      return i;
  }
  jsts.util.Assert.shouldNeverReachHere('Class not supported: ' + this);
  return -1;
};

jsts.geom.Geometry.prototype.toString = function() {
  return new jsts.io.WKTWriter().write(this);
};

/**
 * @return {Point}
 * @private
 */
jsts.geom.Geometry.prototype.createPointFromInternalCoord = function(coord,
    exemplar) {
  exemplar.getPrecisionModel().makePrecise(coord);
  return exemplar.getFactory().createPoint(coord);
};




// JSTS Coordinate
jsts.geom.Coordinate = function(x, y) {
    if (typeof x === 'number') {
      this.x = x;
      this.y = y;
    } else if (x instanceof jsts.geom.Coordinate) {
      this.x = parseFloat(x.x);
      this.y = parseFloat(x.y);
    } else if (x === undefined || x === null) {
      this.x = 0;
      this.y = 0;
    } else if (typeof x === 'string') {
      this.x = parseFloat(x);
      this.y = parseFloat(y);
    }
  };

  /**
   * Sets this <code>Coordinate</code>s (x,y,z) values to that of
   * <code>other</code>.
   *
   * @param {Coordinate}
   *          other the <code>Coordinate</code> to copy.
   */
  jsts.geom.Coordinate.prototype.setCoordinate = function(other) {
    this.x = other.x;
    this.y = other.y;
  };


  /**
   * Clones this instance.
   *
   * @return {Coordinate} A point instance cloned from this.
   */
  jsts.geom.Coordinate.prototype.clone = function() {
    return new jsts.geom.Coordinate(this.x, this.y);
  };


  /**
   * Computes the 2-dimensional Euclidean distance to another location. The
   * Z-ordinate is ignored.
   *
   * @param {Coordinate}
   *          p a point.
   * @return {number} the 2-dimensional Euclidean distance between the
   *         locations.
   */
  jsts.geom.Coordinate.prototype.distance = function(p) {
    var dx = this.x - p.x;
    var dy = this.y - p.y;

    return Math.sqrt(dx * dx + dy * dy);
  };

  /**
   * Returns whether the planar projections of the two <code>Coordinate</code>s
   * are equal.
   *
   * @param {Coordinate}
   *          other a <code>Coordinate</code> with which to do the 2D
   *          comparison.
   * @return {boolean} <code>true</code> if the x- and y-coordinates are
   *         equal; the z-coordinates do not have to be equal.
   */
  jsts.geom.Coordinate.prototype.equals2D = function(other) {
    if (this.x !== other.x) {
      return false;
    }

    if (this.y !== other.y) {
      return false;
    }

    return true;
  };

  /**
   * Returns <code>true</code> if <code>other</code> has the same values for
   * the x and y ordinates. Since Coordinates are 2.5D, this routine ignores the
   * z value when making the comparison.
   *
   * @param {Coordinate}
   *          other a <code>Coordinate</code> with which to do the comparison.
   * @return {boolean} <code>true</code> if <code>other</code> is a
   *         <code>Coordinate</code> with the same values for the x and y
   *         ordinates.
   */
  jsts.geom.Coordinate.prototype.equals = function(other) {
    if (!other instanceof jsts.geom.Coordinate || other === undefined) {
      return false;
    }
    return this.equals2D(other);
  };

  /**
   * Compares this {@link Coordinate} with the specified {@link Coordinate} for
   * order. This method ignores the z value when making the comparison. Returns:
   * <UL>
   * <LI> -1 : this.x < other.x || ((this.x == other.x) && (this.y < other.y))
   * <LI> 0 : this.x == other.x && this.y = other.y
   * <LI> 1 : this.x > other.x || ((this.x == other.x) && (this.y > other.y))
   *
   * </UL>
   * Note: This method assumes that ordinate values are valid numbers. NaN
   * values are not handled correctly.
   *
   * @param {Coordinate}
   *          other the <code>Coordinate</code> with which this
   *          <code>Coordinate</code> is being compared.
   * @return {number} -1, zero, or 1 as explained above.
   */
  jsts.geom.Coordinate.prototype.compareTo = function(other) {
    if (this.x < other.x) {
      return -1;
    }
    if (this.x > other.x) {
      return 1;
    }
    if (this.y < other.y) {
      return -1;
    }
    if (this.y > other.y) {
      return 1;
    }

    return 0;
  };

  jsts.geom.Coordinate.prototype.toString = function() {
    return '(' + this.x + ', ' + this.y + ')';
  };


  //    JSTS Point

  jsts.geom.Point = function(coordinate, factory) {
  this.factory = factory;

  if (coordinate === undefined)
    return;

  this.coordinate = coordinate;
};

jsts.geom.Point.prototype = new jsts.geom.Geometry();
jsts.geom.Point.constructor = jsts.geom.Point;


jsts.geom.Point.CLASS_NAME = 'jsts.geom.Point';


jsts.geom.Point.prototype.coordinate = null;


/**
 * @return {number} x-axis value of this Point.
 */
jsts.geom.Point.prototype.getX = function() {
  return this.coordinate.x;
};


/**
 * @return {number} y-axis value of this Point.
 */
jsts.geom.Point.prototype.getY = function() {
  return this.coordinate.y;
};

jsts.geom.Point.prototype.getCoordinate = function() {
  return this.coordinate;
};


/**
 * @return {Coordinate[]} this Point as coordinate array.
 */
jsts.geom.Point.prototype.getCoordinates = function() {
  return this.isEmpty() ? [] : [this.coordinate];
};

jsts.geom.Point.prototype.isEmpty = function() {
  return this.coordinate === null;
};

jsts.geom.Point.prototype.equalsExact = function(other, tolerance) {
  if (!this.isEquivalentClass(other)) {
    return false;
  }
  if (this.isEmpty() && other.isEmpty()) {
    return true;
  }
  return this.equal(other.getCoordinate(), this.getCoordinate(), tolerance);
};


/**
 * @return {number} number of coordinates (0 or 1).
 */
jsts.geom.Point.prototype.getNumPoints = function() {
  return this.isEmpty() ? 0 : 1;
};


/**
 * @return {boolean} Point is always simple.
 */
jsts.geom.Point.prototype.isSimple = function() {
  return true;
};


/**
 * Gets the boundary of this geometry. Zero-dimensional geometries have no
 * boundary by definition, so an empty GeometryCollection is returned.
 *
 * @return {GeometryCollection} an empty GeometryCollection.
 * @see Geometry#getBoundary
 */
jsts.geom.Point.prototype.getBoundary = function() {
  return new jsts.geom.GeometryCollection(null);
};


/**
 * @return {Envelope} Envelope of this point.
 */
jsts.geom.Point.prototype.computeEnvelopeInternal = function() {
  if (this.isEmpty()) {
    return new jsts.geom.Envelope();
  }
  return new jsts.geom.Envelope(this.coordinate);
};

jsts.geom.Point.prototype.apply = function(filter) {
  if (filter instanceof jsts.geom.GeometryFilter || filter instanceof jsts.geom.GeometryComponentFilter) {
    filter.filter(this);
  } else if (filter instanceof jsts.geom.CoordinateFilter) {
    if (this.isEmpty()) { return; }
    filter.filter(this.getCoordinate());
  }

};

jsts.geom.Point.prototype.clone = function() {
  return new jsts.geom.Point(this.coordinate.clone(), this.factory);
};


/**
 * @return {number} Always 0.
 */
jsts.geom.Point.prototype.getDimension = function() {
  return 0;
};


/**
 * @return {number} Always Dimension.FALSE.
 */
jsts.geom.Point.prototype.getBoundaryDimension = function() {
  return jsts.geom.Dimension.FALSE;
};


/**
 * @return {Point} Reversed point is a cloned point.
 */
jsts.geom.Point.prototype.reverse = function() {
  return this.clone();
};


/**
 * A Point is valid iff:
 * <ul>
 * <li>the coordinate which defines it is a valid coordinate (i.e does not have
 * an NaN X or Y ordinate)
 * </ul>
 *
 * @return {boolean} true iff the Point is valid.
 */
jsts.geom.Point.prototype.isValid = function() {
  if (!jsts.operation.valid.IsValidOp.isValid(this.getCoordinate())) {
    return false;
  }
  return true;
};


/**
 *
 */
jsts.geom.Point.prototype.normalize = function() {
  // a Point is always in normalized form
};

jsts.geom.Point.prototype.compareToSameClass = function(other) {
  var point = other;
  return this.getCoordinate().compareTo(point.getCoordinate());
};

/**
 * @return {string} String representation of Point type.
 */
jsts.geom.Point.prototype.getGeometryType = function() {
  return 'Point';
};

jsts.geom.Point.prototype.hashCode = function() {
  return 'Point_' + this.coordinate.hashCode();
};

jsts.geom.Point.prototype.CLASS_NAME = 'jsts.geom.Point';




//    JSTS Polygon

/**
   * Represents a linear polygon, which may include holes. The shell and holes
   * of the polygon are represented by {@link LinearRing}s. In a valid polygon,
   * holes may touch the shell or other holes at a single point. However, no
   * sequence of touching holes may split the polygon into two pieces. The
   * orientation of the rings in the polygon does not matter.
   *
   * The shell and holes must conform to the assertions specified in the <A
   * HREF="http://www.opengis.org/techno/specs.htm">OpenGIS Simple Features
   * Specification for SQL</A>.
   */


  /**
   * @requires jsts/geom/Geometry.js
   */

  /**
   * @extends {jsts.geom.Geometry}
   * @constructor
   */
  jsts.geom.Polygon = function(shell, holes, factory) {
    this.shell = shell || factory.createLinearRing(null);
    this.holes = holes || [];
    this.factory = factory;
  };

  jsts.geom.Polygon.prototype = new jsts.geom.Geometry();
  jsts.geom.Polygon.constructor = jsts.geom.Polygon;

  jsts.geom.Polygon.prototype.getCoordinate = function() {
    return this.shell.getCoordinate();
  };

  jsts.geom.Polygon.prototype.getCoordinates = function() {
    if (this.isEmpty()) {
      return [];
    }
    var coordinates = [];
    var k = -1;
    var shellCoordinates = this.shell.getCoordinates();
    for (var x = 0; x < shellCoordinates.length; x++) {
      k++;
      coordinates[k] = shellCoordinates[x];
    }
    for (var i = 0; i < this.holes.length; i++) {
      var childCoordinates = this.holes[i].getCoordinates();
      for (var j = 0; j < childCoordinates.length; j++) {
        k++;
        coordinates[k] = childCoordinates[j];
      }
    }
    return coordinates;
  };

  /**
   * @return {number}
   */
  jsts.geom.Polygon.prototype.getNumPoints = function() {
    var numPoints = this.shell.getNumPoints();
    for (var i = 0; i < this.holes.length; i++) {
      numPoints += this.holes[i].getNumPoints();
    }
    return numPoints;
  };

  /**
   * @return {boolean}
   */
  jsts.geom.Polygon.prototype.isEmpty = function() {
    return this.shell.isEmpty();
  };

  jsts.geom.Polygon.prototype.getExteriorRing = function() {
    return this.shell;
  };

  jsts.geom.Polygon.prototype.getInteriorRingN = function(n) {
    return this.holes[n];
  };

  jsts.geom.Polygon.prototype.getNumInteriorRing = function() {
    return this.holes.length;
  };

  /**
   * Returns the area of this <code>Polygon</code>
   *
   * @return the area of the polygon.
   */
  jsts.geom.Polygon.prototype.getArea = function() {
    var area = 0.0;
    area += Math.abs(jsts.algorithm.CGAlgorithms.signedArea(this.shell
        .getCoordinateSequence()));
    for (var i = 0; i < this.holes.length; i++) {
      area -= Math.abs(jsts.algorithm.CGAlgorithms.signedArea(this.holes[i]
          .getCoordinateSequence()));
    }
    return area;
  };

  /**
   * Returns the perimeter of this <code>Polygon</code>
   *
   * @return the perimeter of the polygon.
   */
  jsts.geom.Polygon.prototype.getLength = function() {
    var len = 0.0;
    len += this.shell.getLength();
    for (var i = 0; i < this.holes.length; i++) {
      len += this.holes[i].getLength();
    }
    return len;
  };

  /**
   * Computes the boundary of this geometry
   *
   * @return {Geometry} a lineal geometry (which may be empty).
   * @see Geometry#getBoundary
   */
  jsts.geom.Polygon.prototype.getBoundary = function() {
    if (this.isEmpty()) {
      return this.getFactory().createMultiLineString(null);
    }
    var rings = [];
    rings[0] = this.shell.clone();
    for (var i = 0, len = this.holes.length; i < len; i++) {
      rings[i + 1] = this.holes[i].clone();
    }
    // create LineString or MultiLineString as appropriate
    if (rings.length <= 1)
      return rings[0];
    return this.getFactory().createMultiLineString(rings);
  };

  jsts.geom.Polygon.prototype.computeEnvelopeInternal = function() {
    return this.shell.getEnvelopeInternal();
  };

  jsts.geom.Polygon.prototype.getDimension = function() {
    return 2;
  };

  jsts.geom.Polygon.prototype.getBoundaryDimension = function() {
    return 1;
  };


  /**
   * @param {Geometry}
   *          other
   * @param {number}
   *          tolerance
   * @return {boolean}
   */
  jsts.geom.Polygon.prototype.equalsExact = function(other, tolerance) {
    if (!this.isEquivalentClass(other)) {
      return false;
    }
    if (this.isEmpty() && other.isEmpty()) {
      return true;
    }
    if (this.isEmpty() !== other.isEmpty()) {
      return false;
    }

    if (!this.shell.equalsExact(other.shell, tolerance)) {
      return false;
    }
    if (this.holes.length !== other.holes.length) {
      return false;
    }
    if (this.holes.length !== other.holes.length) {
      return false;
    }
    for (var i = 0; i < this.holes.length; i++) {
      if (!(this.holes[i]).equalsExact(other.holes[i], tolerance)) {
        return false;
      }
    }
    return true;
  };

  jsts.geom.Polygon.prototype.compareToSameClass = function(o) {
    return this.shell.compareToSameClass(o.shell);
  };

  jsts.geom.Polygon.prototype.apply = function(filter) {
    if (filter instanceof jsts.geom.GeometryComponentFilter) {
      filter.filter(this);
      this.shell.apply(filter);
      for (var i = 0, len = this.holes.length; i < len; i++) {
        this.holes[i].apply(filter);
      }
    } else if (filter instanceof jsts.geom.GeometryFilter) {
      filter.filter(this);
    } else if (filter instanceof jsts.geom.CoordinateFilter) {
      this.shell.apply(filter);
      for (var i = 0, len = this.holes.length; i < len; i++) {
        this.holes[i].apply(filter);
      }
    } else if (filter instanceof jsts.geom.CoordinateSequenceFilter) {
      this.apply2.apply(this, arguments);
    }
  };

  jsts.geom.Polygon.prototype.apply2 = function(filter) {
    this.shell.apply(filter);
    if (!filter.isDone()) {
      for (var i = 0; i < this.holes.length; i++) {
        this.holes[i].apply(filter);
        if (filter.isDone())
          break;
      }
    }
    if (filter.isGeometryChanged()) {
      // TODO: call this.geometryChanged(); when ported
    }
  };

  /**
   * Creates and returns a full copy of this {@link Polygon} object. (including
   * all coordinates contained by it).
   *
   * @return a clone of this instance.
   */
  jsts.geom.Polygon.prototype.clone = function() {
    var holes = [];

    for (var i = 0, len = this.holes.length; i < len; i++) {
      holes.push(this.holes[i].clone());
    }

    return this.factory.createPolygon(this.shell.clone(), holes);
  };

  jsts.geom.Polygon.prototype.normalize = function() {
    this.normalize2(this.shell, true);
    for (var i = 0, len = this.holes.length; i < len; i++) {
      this.normalize2(this.holes[i], false);
    }
    // TODO: might need to supply comparison function
    this.holes.sort();
  };

  /**
   * @private
   */
  jsts.geom.Polygon.prototype.normalize2 = function(ring, clockwise) {
    if (ring.isEmpty()) {
      return;
    }
    var uniqueCoordinates = ring.points.slice(0, ring.points.length - 1);
    var minCoordinate = jsts.geom.CoordinateArrays.minCoordinate(ring.points);
    jsts.geom.CoordinateArrays.scroll(uniqueCoordinates, minCoordinate);
    ring.points = uniqueCoordinates.concat();
    ring.points[uniqueCoordinates.length] = uniqueCoordinates[0];
    if (jsts.algorithm.CGAlgorithms.isCCW(ring.points) === clockwise) {
      ring.points.reverse();
    }
  };

  /**
   * @return {String} String representation of Polygon type.
   */
  jsts.geom.Polygon.prototype.getGeometryType = function() {
    return 'Polygon';
  };
  
  jsts.geom.Polygon.prototype.CLASS_NAME = 'jsts.geom.Polygon';





//    JSTS LineString


  /**
   * @requires jsts/geom/Geometry.js
   * @requires jsts/geom/Dimension.js
   */

  var Dimension = jsts.geom.Dimension;

  /**
   * @extends jsts.geom.Geometry
   * @constructor
   */
  jsts.geom.LineString = function(points, factory) {
    this.factory = factory;
    this.points = points || [];
  };

  jsts.geom.LineString.prototype = new jsts.geom.Geometry();
  jsts.geom.LineString.constructor = jsts.geom.LineString;

  /**
   * @type {jsts.geom.Coordinate[]}
   * @private
   */
  jsts.geom.LineString.prototype.points = null;

  /**
   * @return {jsts.geom.Coordinate[]} this LineString's internal coordinate
   *         array.
   */
  jsts.geom.LineString.prototype.getCoordinates = function() {
    return this.points;
  };

  jsts.geom.LineString.prototype.getCoordinateSequence = function() {
    return this.points;
  };


  /**
   * @return {jsts.geom.Coordinate} The n'th coordinate of this
   *         jsts.geom.LineString.
   * @param {int}
   *          n index.
   */
  jsts.geom.LineString.prototype.getCoordinateN = function(n) {
    return this.points[n];
  };


  /**
   * @return {jsts.geom.Coordinate} The first coordinate of this LineString or
   *         null if empty.
   */
  jsts.geom.LineString.prototype.getCoordinate = function() {
    if (this.isEmpty()) {
      return null;
    }
    return this.getCoordinateN(0);
  };


  /**
   * @return {number} LineStrings are always 1-dimensional.
   */
  jsts.geom.LineString.prototype.getDimension = function() {
    return 1;
  };


  /**
   * @return {number} dimension of the boundary of this jsts.geom.LineString.
   */
  jsts.geom.LineString.prototype.getBoundaryDimension = function() {
    if (this.isClosed()) {
      return Dimension.FALSE;
    }
    return 0;
  };


  /**
   * @return {Boolean} true if empty.
   */
  jsts.geom.LineString.prototype.isEmpty = function() {
    return this.points.length === 0;
  };

  jsts.geom.LineString.prototype.getNumPoints = function() {
    return this.points.length;
  };

  jsts.geom.LineString.prototype.getPointN = function(n) {
    return this.getFactory().createPoint(this.points[n]);
  };


  jsts.geom.LineString.prototype.getStartPoint = function() {
    if (this.isEmpty()) {
      return null;
    }
    return this.getPointN(0);
  };

  jsts.geom.LineString.prototype.getEndPoint = function() {
    if (this.isEmpty()) {
      return null;
    }
    return this.getPointN(this.getNumPoints() - 1);
  };


  /**
   * @return {Boolean} true if LineString is Closed.
   */
  jsts.geom.LineString.prototype.isClosed = function() {
    if (this.isEmpty()) {
      return false;
    }
    return this.getCoordinateN(0).equals2D(
        this.getCoordinateN(this.points.length - 1));
  };


  /**
   * @return {Boolean} true if LineString is a Ring.
   */
  jsts.geom.LineString.prototype.isRing = function() {
    return this.isClosed() && this.isSimple();
  };


  /**
   * @return {String} String representation of LineString type.
   */
  jsts.geom.LineString.prototype.getGeometryType = function() {
    return 'LineString';
  };


  /**
   * Returns the length of this <code>LineString</code>
   *
   * @return the length of the linestring.
   */
  jsts.geom.LineString.prototype.getLength = function() {
    return jsts.algorithm.CGAlgorithms.computeLength(this.points);
  };

  /**
   * Gets the boundary of this geometry. The boundary of a lineal geometry is
   * always a zero-dimensional geometry (which may be empty).
   *
   * @return {Geometry} the boundary geometry.
   * @see Geometry#getBoundary
   */
  jsts.geom.LineString.prototype.getBoundary = function() {
    return (new jsts.operation.BoundaryOp(this)).getBoundary();
  };


  jsts.geom.LineString.prototype.computeEnvelopeInternal = function() {
    if (this.isEmpty()) {
      return new jsts.geom.Envelope();
    }

    var env = new jsts.geom.Envelope();
    this.points.forEach(function(component) {
      env.expandToInclude(component);
    });

    return env;
  };


  /**
   * @param {Geometry}
   *          other Geometry to compare this LineString to.
   * @param {double}
   *          tolerance Tolerance.
   * @return {Boolean} true if equal.
   */
  jsts.geom.LineString.prototype.equalsExact = function(other, tolerance) {
    if (!this.isEquivalentClass(other)) {
      return false;
    }

    if (this.points.length !== other.points.length) {
      return false;
    }

    if (this.isEmpty() && other.isEmpty()) {
      return true;
    }

    return this.points
        .reduce(function(equal, point, i) {
          return equal &&
              jsts.geom.Geometry.prototype.equal(point, other.points[i],
                  tolerance);
        });
  };

  jsts.geom.LineString.prototype.isEquivalentClass = function(other) {
    return other instanceof jsts.geom.LineString;
  };

  jsts.geom.LineString.prototype.compareToSameClass = function(o) {
    var line = o;
    // MD - optimized implementation
    var i = 0, il = this.points.length;
    var j = 0, jl = line.points.length;
    while (i < il && j < jl) {
      var comparison = this.points[i].compareTo(line.points[j]);
      if (comparison !== 0) {
        return comparison;
      }
      i++;
      j++;
    }
    if (i < il) {
      return 1;
    }
    if (j < jl) {
      return -1;
    }
    return 0;
  };

  jsts.geom.LineString.prototype.apply = function(filter) {
    if (filter instanceof jsts.geom.GeometryFilter ||
        filter instanceof jsts.geom.GeometryComponentFilter) {
      filter.filter(this);
    } else if (filter instanceof jsts.geom.CoordinateFilter) {
      for (var i = 0, len = this.points.length; i < len; i++) {
        filter.filter(this.points[i]);
      }
    } else if (filter instanceof jsts.geom.CoordinateSequenceFilter) {
      this.apply2.apply(this, arguments);
    }
  };

  jsts.geom.LineString.prototype.apply2 = function(filter) {
    if (this.points.length === 0)
      return;
    for (var i = 0; i < this.points.length; i++) {
      filter.filter(this.points, i);
      if (filter.isDone())
        break;
    }
    if (filter.isGeometryChanged()) {
      // TODO: call geometryChanged(); when ported
    }
  };

  jsts.geom.LineString.prototype.clone = function() {
    var points = [];

    for (var i = 0, len = this.points.length; i < len; i++) {
      points.push(this.points[i].clone());
    }

    return this.factory.createLineString(points);
  };

  /**
   * Normalizes a LineString. A normalized linestring has the first point which
   * is not equal to it's reflected point less than the reflected point.
   */
  jsts.geom.LineString.prototype.normalize = function() {
    var i, il, j, ci, cj, len;

    len = this.points.length;
    il = parseInt(len / 2);

    for (i = 0; i < il; i++) {
      j = len - 1 - i;
      // skip equal points on both ends
      ci = this.points[i];
      cj = this.points[j];
      if (!ci.equals(cj)) {
        if (ci.compareTo(cj) > 0) {
          this.points.reverse();
        }
        return;
      }
    }
  };

  jsts.geom.LineString.prototype.CLASS_NAME = 'jsts.geom.LineString';




//    JSTS





//    JSTS BUFFER PARAMETERS

/* Copyright (c) 2011 by The Authors.
 * Published under the LGPL 2.1 license.
 * See /license-notice.txt for the full text of the license notice.
 * See /license.txt for the full text of the license.
 */


/**
 * Contains the parameters which describe how a buffer should be constructed.
 *
 * @constructor
 */
jsts.operation.buffer.BufferParameters = function(quadrantSegments,
    endCapStyle, joinStyle, mitreLimit) {
  if (quadrantSegments)
    this.setQuadrantSegments(quadrantSegments);
  if (endCapStyle)
    this.setEndCapStyle(endCapStyle);
  if (joinStyle)
    this.setJoinStyle(joinStyle);
  if (mitreLimit)
    this.setMitreLimit(mitreLimit);
};


/**
 * Specifies a round line buffer end cap style.
 *
 * @type {int}
 */
jsts.operation.buffer.BufferParameters.CAP_ROUND = 1;


/**
 * Specifies a flat line buffer end cap style.
 *
 * @type {int}
 */
jsts.operation.buffer.BufferParameters.CAP_FLAT = 2;


/**
 * Specifies a square line buffer end cap style.
 *
 * @type {int}
 */
jsts.operation.buffer.BufferParameters.CAP_SQUARE = 3;


/**
 * Specifies a round join style.
 *
 * @type {int}
 */
jsts.operation.buffer.BufferParameters.JOIN_ROUND = 1;


/**
 * Specifies a mitre join style.
 */
jsts.operation.buffer.BufferParameters.JOIN_MITRE = 2;


/**
 * Specifies a bevel join style.
 *
 * @type {int}
 */
jsts.operation.buffer.BufferParameters.JOIN_BEVEL = 3;


/**
 * The default number of facets into which to divide a fillet of 90 degrees. A
 * value of 8 gives less than 2% max error in the buffer distance. For a max
 * error of < 1%, use QS = 12. For a max error of < 0.1%, use QS = 18.
 *
 * @type {int}
 */
jsts.operation.buffer.BufferParameters.DEFAULT_QUADRANT_SEGMENTS = 8;


/**
 * The default mitre limit Allows fairly pointy mitres.
 *
 * @type {double}
 */
jsts.operation.buffer.BufferParameters.DEFAULT_MITRE_LIMIT = 5.0;


/**
 * @type {int}
 * @private
 */
jsts.operation.buffer.BufferParameters.prototype.quadrantSegments = jsts.operation.buffer.BufferParameters.DEFAULT_QUADRANT_SEGMENTS;


/**
 * @type {int}
 * @private
 */
jsts.operation.buffer.BufferParameters.prototype.endCapStyle = jsts.operation.buffer.BufferParameters.CAP_ROUND;


/**
 * @type {int}
 * @private
 */
jsts.operation.buffer.BufferParameters.prototype.joinStyle = jsts.operation.buffer.BufferParameters.JOIN_ROUND;


/**
 * @type {double}
 * @private
 */
jsts.operation.buffer.BufferParameters.prototype.mitreLimit = jsts.operation.buffer.BufferParameters.DEFAULT_MITRE_LIMIT;

/**
 * @type {boolean}
 * @private
 */
jsts.operation.buffer.BufferParameters.prototype._isSingleSided = false;

/**
 * Gets the number of quadrant segments which will be used
 *
 * @return the number of quadrant segments.
 */
jsts.operation.buffer.BufferParameters.prototype.getQuadrantSegments = function() {
  return this.quadrantSegments;
};


/**
 * Sets the number of segments used to approximate a angle fillet
 *
 * @param {int}
 *          quadrantSegments the number of segments in a fillet for a quadrant.
 */
jsts.operation.buffer.BufferParameters.prototype.setQuadrantSegments = function(
    quadrantSegments) {
  this.quadrantSegments = quadrantSegments;
};


/**
 * Sets the number of line segments used to approximate an angle fillet.
 * <ul>
 * <li>If <tt>quadSegs</tt> >= 1, joins are round, and <tt>quadSegs</tt>
 * indicates the number of segments to use to approximate a quarter-circle.
 * <li>If <tt>quadSegs</tt> = 0, joins are bevelled (flat)
 * <li>If <tt>quadSegs</tt> < 0, joins are mitred, and the value of qs
 * indicates the mitre ration limit as
 *
 * <pre>
 * mitreLimit= |
 * <tt>
 * quadSegs
 * </tt>
 * |
 * </pre>
 *
 * </ul>
 * For round joins, <tt>quadSegs</tt> determines the maximum error in the
 * approximation to the true buffer curve. The default value of 8 gives less
 * than 2% max error in the buffer distance. For a max error of < 1%, use QS =
 * 12. For a max error of < 0.1%, use QS = 18. The error is always less than the
 * buffer distance (in other words, the computed buffer curve is always inside
 * the true curve).
 *
 * @param quadrantSegments
 *          the number of segments in a fillet for a quadrant.
 */
jsts.operation.buffer.BufferParameters.prototype.setQuadrantSegments = function(
    quadSegs) {
  this.quadrantSegments = quadSegs;

  /**
   * Indicates how to construct fillets. If qs >= 1, fillet is round, and qs
   * indicates number of segments to use to approximate a quarter-circle. If qs =
   * 0, fillet is bevelled flat (i.e. no filleting is performed) If qs < 0,
   * fillet is mitred, and absolute value of qs indicates maximum length of
   * mitre according to
   *
   * mitreLimit = |qs|
   */
  if (this.quadrantSegments === 0)
    this.joinStyle = jsts.operation.buffer.BufferParameters.JOIN_BEVEL;
  if (this.quadrantSegments < 0) {
    this.joinStyle = jsts.operation.buffer.BufferParameters.JOIN_MITRE;
    this.mitreLimit = Math.abs(this.quadrantSegments);
  }

  if (quadSegs <= 0) {
    this.quadrantSegments = 1;
  }

  /**
   * If join style was set by the quadSegs value, use the default for the actual
   * quadrantSegments value.
   */
  if (this.joinStyle !== jsts.operation.buffer.BufferParameters.JOIN_ROUND) {
    this.quadrantSegments = jsts.operation.buffer.BufferParameters.DEFAULT_QUADRANT_SEGMENTS;
  }
};


/**
 * Computes the maximum distance error due to a given level of approximation to
 * a true arc.
 *
 * @param quadSegs
 *          the number of segments used to approximate a quarter-circle.
 * @return the error of approximation.
 */
jsts.operation.buffer.BufferParameters.bufferDistanceError = function(quadSegs) {
  var alpha = Math.PI / 2.0 / quadSegs;
  return 1 - Math.cos(alpha / 2.0);
};


/**
 * Gets the end cap style.
 *
 * @return the end cap style.
 */
jsts.operation.buffer.BufferParameters.prototype.getEndCapStyle = function() {
  return this.endCapStyle;
};


/**
 * Specifies the end cap style of the generated buffer. The styles supported are
 * {@link #CAP_ROUND}, {@link #CAP_BUTT}, and {@link #CAP_SQUARE}. The
 * default is CAP_ROUND.
 *
 * @param {int}
 *          endCapStyle the end cap style to specify.
 */
jsts.operation.buffer.BufferParameters.prototype.setEndCapStyle = function(
    endCapStyle) {
  this.endCapStyle = endCapStyle;
};


/**
 * Gets the join style
 *
 * @return the join style code.
 */
jsts.operation.buffer.BufferParameters.prototype.getJoinStyle = function() {
  return this.joinStyle;
};


/**
 * Sets the join style for outside (reflex) corners between line segments.
 * Allowable values are {@link JOIN_ROUND} (which is the default),
 * {@link JOIN_MITRE} and {link JOIN_BEVEL}.
 *
 * @param joinStyle
 *          the code for the join style.
 */
jsts.operation.buffer.BufferParameters.prototype.setJoinStyle = function(
    joinStyle) {
  this.joinStyle = joinStyle;
};


/**
 * Gets the mitre ratio limit.
 *
 * @return the limit value.
 */
jsts.operation.buffer.BufferParameters.prototype.getMitreLimit = function() {
  return this.mitreLimit;
};


/**
 * Sets the limit on the mitre ratio used for very sharp corners. The mitre
 * ratio is the ratio of the distance from the corner to the end of the mitred
 * offset corner. When two line segments meet at a sharp angle, a miter join
 * will extend far beyond the original geometry. (and in the extreme case will
 * be infinitely far.) To prevent unreasonable geometry, the mitre limit allows
 * controlling the maximum length of the join corner. Corners with a ratio which
 * exceed the limit will be beveled.
 *
 * @param mitreLimit
 *          the mitre ratio limit.
 */
jsts.operation.buffer.BufferParameters.prototype.setMitreLimit = function(
    mitreLimit) {
  this.mitreLimit = mitreLimit;
};


/**
 * Sets whether the computed buffer should be single-sided. A single-sided
 * buffer is constructed on only one side of each input line.
 * <p>
 * The side used is determined by the sign of the buffer distance:
 * <ul>
 * <li>a positive distance indicates the left-hand side
 * <li>a negative distance indicates the right-hand side
 * </ul>
 * The single-sided buffer of point geometries is the same as the regular
 * buffer.
 * <p>
 * The End Cap Style for single-sided buffers is always ignored, and forced to
 * the equivalent of <tt>CAP_FLAT</tt>.
 *
 * @param isSingleSided
 *          true if a single-sided buffer should be constructed.
 */
jsts.operation.buffer.BufferParameters.prototype.setSingleSided = function(
    isSingleSided) {
  this._isSingleSided = isSingleSided;
};


/**
 * Tests whether the buffer is to be generated on a single side only.
 *
 * @return true if the generated buffer is to be single-sided.
 */
jsts.operation.buffer.BufferParameters.prototype.isSingleSided = function() {
  return this._isSingleSided;
};



//    JSTS BufferOp

/* Copyright (c) 2011 by The Authors.
 * Published under the LGPL 2.1 license.
 * See /license-notice.txt for the full text of the license notice.
 * See /license.txt for the full text of the license.
 */

/**
 * Computes the buffer of a geometry, for both positive and negative buffer
 * distances.
 *
 * In GIS, the positive buffer of a geometry is defined as
 * the Minkowski sum or difference of the geometry
 * with a circle of radius equal to the absolute value of the buffer distance.
 * In the CAD/CAM world buffers are known as </i>offset curves</i>.
 * In morphological analysis they are known as <i>erosion</i> and
 * <i>dilation</i>
 *
 * The buffer operation always returns a polygonal result.
 * The negative or zero-distance buffer of lines and points is always an empty
 * {@link Polygon}.
 *
 * Since true buffer curves may contain circular arcs,
 * computed buffer polygons can only be approximations to the true geometry.
 * The user can control the accuracy of the curve approximation by specifying
 * the number of linear segments used to approximate curves.
 *
 * The <b>end cap style</b> of a linear buffer may be specified. The
 * following end cap styles are supported:
 * <ul
 * <li>{@link #CAP_ROUND} - the usual round end caps
 * <li>{@link #CAP_BUTT} - end caps are truncated flat at the line ends
 * <li>{@link #CAP_SQUARE} - end caps are squared off at the buffer distance
 * beyond the line ends
 * </ul>
 *
 */



/**
 * Initializes a buffer computation for the given geometry with the given set of
 * parameters.
 *
 * @param {Geometry}
 *          g the geometry to buffer.
 * @param {BufferParameters}
 *          bufParams the buffer parameters to use.
 * @constructor
 */
jsts.operation.buffer.BufferOp = function(g, bufParams) {
  this.argGeom = g;
  this.bufParams = bufParams ? bufParams
      : new jsts.operation.buffer.BufferParameters();
};


/**
 * A number of digits of precision which leaves some computational "headroom"
 * for floating point operations.
 *
 * This value should be less than the decimal precision of double-precision
 * values (16).
 *
 * @type {int}
 */
jsts.operation.buffer.BufferOp.MAX_PRECISION_DIGITS = 12;


/**
 * Compute a scale factor to limit the precision of a given combination of
 * Geometry and buffer distance. The scale factor is determined by a combination
 * of the number of digits of precision in the (geometry + buffer distance),
 * limited by the supplied <code>maxPrecisionDigits</code> value.
 *
 * @param {Geometry}
 *          g the Geometry being buffered.
 * @param {double}
 *          distance the buffer distance.
 * @param {int}
 *          maxPrecisionDigits the max # of digits that should be allowed by the
 *          precision determined by the computed scale factor.
 *
 * @return {double} a scale factor for the buffer computation.
 */
jsts.operation.buffer.BufferOp.precisionScaleFactor = function(g, distance,
    maxPrecisionDigits) {
  var env = g.getEnvelopeInternal();
  var envSize = Math.max(env.getHeight(), env.getWidth());
  var expandByDistance = distance > 0.0 ? distance : 0.0;
  var bufEnvSize = envSize + 2 * expandByDistance;

  // the smallest power of 10 greater than the buffer envelope
  var bufEnvLog10 = (Math.log(bufEnvSize) / Math.log(10) + 1.0);
  var minUnitLog10 = bufEnvLog10 - maxPrecisionDigits;
  // scale factor is inverse of min Unit size, so flip sign of exponent
  var scaleFactor = Math.pow(10.0, -minUnitLog10);
  return scaleFactor;
};


/**
 * Computes the buffer of a geometry for a given buffer distance.
 *
 * @param {Geometry}
 *          g the geometry to buffer.
 * @param {double}
 *          distance the buffer distance.
 * @return {Geometry} the buffer of the input geometry.
 */
jsts.operation.buffer.BufferOp.bufferOp = function(g, distance) {
  if (arguments.length > 2) {
    return jsts.operation.buffer.BufferOp.bufferOp2.apply(this, arguments);
  }

  var gBuf = new jsts.operation.buffer.BufferOp(g);
  var geomBuf = gBuf.getResultGeometry(distance);
  return geomBuf;
};


/**
 * Computes the buffer for a geometry for a given buffer distance and accuracy
 * of approximation.
 *
 * @param {Geometry}
 *          g the geometry to buffer.
 * @param {double}
 *          distance the buffer distance.
 * @param {BufferParameters}
 *          params the buffer parameters to use.
 * @return {Geometry} the buffer of the input geometry.
 *
 */
jsts.operation.buffer.BufferOp.bufferOp2 = function(g, distance, params) {
  if (arguments.length > 3) {
    return jsts.operation.buffer.BufferOp.bufferOp3.apply(this, arguments);
  }

  var bufOp = new jsts.operation.buffer.BufferOp(g, params);
  var geomBuf = bufOp.getResultGeometry(distance);
  return geomBuf;
};


/**
 * Computes the buffer for a geometry for a given buffer distance and accuracy
 * of approximation.
 *
 * @param {Geometry}
 *          g the geometry to buffer.
 * @param {double}
 *          distance the buffer distance.
 * @param {int}
 *          quadrantSegments the number of segments used to approximate a
 *          quarter circle.
 * @return {Geometry} the buffer of the input geometry.
 *
 */
jsts.operation.buffer.BufferOp.bufferOp3 = function(g, distance,
    quadrantSegments) {
  if (arguments.length > 4) {
    return jsts.operation.buffer.BufferOp.bufferOp4.apply(this, arguments);
  }

  var bufOp = new jsts.operation.buffer.BufferOp(g);
  bufOp.setQuadrantSegments(quadrantSegments);
  var geomBuf = bufOp.getResultGeometry(distance);
  return geomBuf;
};


/**
 * Computes the buffer for a geometry for a given buffer distance and accuracy
 * of approximation.
 *
 * @param {Geometry}
 *          g the geometry to buffer.
 * @param {double}
 *          distance the buffer distance.
 * @param {int}
 *          quadrantSegments the number of segments used to approximate a
 *          quarter circle.
 * @param {int}
 *          endCapStyle the end cap style to use.
 * @return {Geometry} the buffer of the input geometry.
 *
 */
jsts.operation.buffer.BufferOp.bufferOp4 = function(g, distance,
    quadrantSegments, endCapStyle) {
  var bufOp = new jsts.operation.buffer.BufferOp(g);
  bufOp.setQuadrantSegments(quadrantSegments);
  bufOp.setEndCapStyle(endCapStyle);
  var geomBuf = bufOp.getResultGeometry(distance);
  return geomBuf;
};


/**
 * @type {Geometry}
 */
jsts.operation.buffer.BufferOp.prototype.argGeom = null;


/**
 * @type {double}
 */
jsts.operation.buffer.BufferOp.prototype.distance = null;


/**
 * @type {BufferParameters}
 */
jsts.operation.buffer.BufferOp.prototype.bufParams = null;


/**
 * @type {Geometry}
 */
jsts.operation.buffer.BufferOp.prototype.resultGeometry = null;


/**
 * Specifies the end cap style of the generated buffer. The styles supported are
 * {@link #CAP_ROUND}, {@link #CAP_BUTT}, and {@link #CAP_SQUARE}. The
 * default is CAP_ROUND.
 *
 * @param {int}
 *          endCapStyle the end cap style to specify.
 */
jsts.operation.buffer.BufferOp.prototype.setEndCapStyle = function(endCapStyle) {
  this.bufParams.setEndCapStyle(endCapStyle);
};


/**
 * Sets the number of segments used to approximate a angle fillet
 *
 * @param {int}
 *          quadrantSegments the number of segments in a fillet for a quadrant.
 */
jsts.operation.buffer.BufferOp.prototype.setQuadrantSegments = function(
    quadrantSegments) {
  this.bufParams.setQuadrantSegments(quadrantSegments);
};


/**
 * Returns the buffer computed for a geometry for a given buffer distance.
 *
 * @param {double}
 *          dist the buffer distance.
 * @return {Geometry} the buffer of the input geometry.
 */
jsts.operation.buffer.BufferOp.prototype.getResultGeometry = function(dist) {
  this.distance = dist;
  this.computeGeometry();
  return this.resultGeometry;
};

jsts.operation.buffer.BufferOp.prototype.computeGeometry = function() {
  this.bufferOriginalPrecision();
  if (this.resultGeometry !== null) {
    return;
  }

  var argPM = this.argGeom.getPrecisionModel();
  if (argPM.getType() === jsts.geom.PrecisionModel.FIXED) {
    this.bufferFixedPrecision(argPM);
  } else {
    this.bufferReducedPrecision();
  }
};


jsts.operation.buffer.BufferOp.prototype.bufferReducedPrecision = function() {
  var precDigits;
  var saveException = null;

  // try and compute with decreasing precision
  for (precDigits = jsts.operation.buffer.BufferOp.MAX_PRECISION_DIGITS; precDigits >= 0; precDigits--) {
    try {
      this.bufferReducedPrecision2(precDigits);
    } catch (/* TopologyException */ex) {
      saveException = ex;
      // don't propagate the exception - it will be detected by fact that
      // resultGeometry is null
    }
    if (this.resultGeometry !== null) {
      return;
    }
  }

  // tried everything - have to bail
  throw saveException;
};


jsts.operation.buffer.BufferOp.prototype.bufferOriginalPrecision = function() {
  try {
    // use fast noding by default
    var bufBuilder = new jsts.operation.buffer.BufferBuilder(this.bufParams);
    this.resultGeometry = bufBuilder.buffer(this.argGeom, this.distance);
  } catch (e) {
    // don't propagate the exception - it will be detected by fact that
    // resultGeometry is null
  }
};


/**
 * @param {int}
 *          precisionDigits
 */
jsts.operation.buffer.BufferOp.prototype.bufferReducedPrecision2 = function(
    precisionDigits) {

  var sizeBasedScaleFactor = jsts.operation.buffer.BufferOp
      .precisionScaleFactor(this.argGeom, this.distance, precisionDigits);

  var fixedPM = new jsts.geom.PrecisionModel(sizeBasedScaleFactor);
  this.bufferFixedPrecision(fixedPM);
};


/**
 * @param {PrecisionModel}
 *          fixedPM
 */
jsts.operation.buffer.BufferOp.prototype.bufferFixedPrecision = function(fixedPM) {
  var noder = new jsts.noding.ScaledNoder(
      new jsts.noding.snapround.MCIndexSnapRounder(
          new jsts.geom.PrecisionModel(1.0)), fixedPM.getScale());

  var bufBuilder = new jsts.operation.buffer.BufferBuilder(this.bufParams);
  bufBuilder.setWorkingPrecisionModel(fixedPM);
  bufBuilder.setNoder(noder);
  // this may throw an exception, if robustness errors are encountered
  this.resultGeometry = bufBuilder.buffer(this.argGeom, this.distance);
};




//    JSTS Envelope

/* Copyright (c) 2011 by The Authors.
 * Published under the LGPL 2.1 license.
 * See /license-notice.txt for the full text of the license notice.
 * See /license.txt for the full text of the license.
 */

/**
 * @requires jsts/geom/Coordinate.js
 */

/**
 * Defines a rectangular region of the 2D coordinate plane. It is often used to
 * represent the bounding box of a {@link Geometry}, e.g. the minimum and
 * maximum x and y values of the {@link Coordinate}s.
 * <p>
 * Note that Envelopes support infinite or half-infinite regions, by using the
 * values of <code>Double.POSITIVE_INFINITY</code> and
 * <code>Double.NEGATIVE_INFINITY</code>.
 * <p>
 * When Envelope objects are created or initialized, the supplies extent values
 * are automatically sorted into the correct order.
 */



/**
 * Creates an <code>Envelope</code> for a region defined by maximum and
 * minimum values.
 *
 * @constructor
 */
jsts.geom.Envelope = function() {
  jsts.geom.Envelope.prototype.init.apply(this, arguments);
};


/**
 * the minimum x-coordinate
 *
 * @type {?number}
 */
jsts.geom.Envelope.prototype.minx = null;


/**
 * the maximum x-coordinate
 *
 * @type {?number}
 */
jsts.geom.Envelope.prototype.maxx = null;


/**
 * the minimum y-coordinate
 *
 * @type {?number}
 */
jsts.geom.Envelope.prototype.miny = null;


/**
 * the maximum y-coordinate
 *
 * @type {?number}
 */
jsts.geom.Envelope.prototype.maxy = null;


/**
 * Creates an <code>Envelope</code> for a region defined by maximum and
 * minimum values.
 *
 * Will call appropriate init* method depending on arguments.
 */
jsts.geom.Envelope.prototype.init = function() {
  if (typeof arguments[0] === 'number' && arguments.length === 4) {
    this.initFromValues(arguments[0], arguments[1], arguments[2], arguments[3]);
  } else if (arguments[0] instanceof jsts.geom.Coordinate &&
      arguments.length === 1) {
    this.initFromCoordinate(arguments[0]);
  } else if (arguments[0] instanceof jsts.geom.Coordinate &&
      arguments.length === 2) {
    this.initFromCoordinates(arguments[0], arguments[1]);
  } else if (arguments[0] instanceof jsts.geom.Envelope &&
      arguments.length === 1) {
    this.initFromEnvelope(arguments[0]);
  } else {
    this.setToNull();
  }
};


/**
 * Initialize an <code>Envelope</code> for a region defined by maximum and
 * minimum values.
 *
 * @param {number}
 *          x1 the first x-value.
 * @param {number}
 *          x2 the second x-value.
 * @param {number}
 *          y1 the first y-value.
 * @param {number}
 *          y2 the second y-value.
 */
jsts.geom.Envelope.prototype.initFromValues = function(x1, x2, y1, y2) {
  if (x1 < x2) {
    this.minx = x1;
    this.maxx = x2;
  } else {
    this.minx = x2;
    this.maxx = x1;
  }
  if (y1 < y2) {
    this.miny = y1;
    this.maxy = y2;
  } else {
    this.miny = y2;
    this.maxy = y1;
  }
};


/**
 * Initialize an <code>Envelope</code> to a region defined by two Coordinates.
 *
 * @param {jsts.geom.Coordinate}
 *          p1 the first Coordinate.
 * @param {jsts.geom.Coordinate}
 *          p2 the second Coordinate.
 */
jsts.geom.Envelope.prototype.initFromCoordinates = function(p1, p2) {
  this.initFromValues(p1.x, p2.x, p1.y, p2.y);
};


/**
 * Initialize an <code>Envelope</code> to a region defined by a single
 * Coordinate.
 *
 * @param {jsts.geom.Coordinate}
 *          p the Coordinate.
 */
jsts.geom.Envelope.prototype.initFromCoordinate = function(p) {
  this.initFromValues(p.x, p.x, p.y, p.y);
};


/**
 * Initialize an <code>Envelope</code> from an existing Envelope.
 *
 * @param {jsts.geom.Envelope}
 *          env the Envelope to initialize from.
 */
jsts.geom.Envelope.prototype.initFromEnvelope = function(env) {
  this.minx = env.minx;
  this.maxx = env.maxx;
  this.miny = env.miny;
  this.maxy = env.maxy;
};


/**
 * Makes this <code>Envelope</code> a "null" envelope, that is, the envelope
 * of the empty geometry.
 */
jsts.geom.Envelope.prototype.setToNull = function() {
  this.minx = 0;
  this.maxx = -1;
  this.miny = 0;
  this.maxy = -1;
};


/**
 * Returns <code>true</code> if this <code>Envelope</code> is a "null"
 * envelope.
 *
 * @return {boolean} <code>true</code> if this <code>Envelope</code> is
 *         uninitialized or is the envelope of the empty geometry.
 */
jsts.geom.Envelope.prototype.isNull = function() {
  return this.maxx < this.minx;
};


/**
 * Returns the difference between the maximum and minimum y values.
 *
 * @return {number} max y - min y, or 0 if this is a null <code>Envelope.</code>
 */
jsts.geom.Envelope.prototype.getHeight = function() {
  if (this.isNull()) {
    return 0;
  }
  return this.maxy - this.miny;
};


/**
 * Returns the difference between the maximum and minimum x values.
 *
 * @return {number} max x - min x, or 0 if this is a null <code>Envelope.</code>
 */
jsts.geom.Envelope.prototype.getWidth = function() {
  if (this.isNull()) {
    return 0;
  }
  return this.maxx - this.minx;
};


/**
 * Returns the <code>Envelope</code>s minimum x-value. min x > max x
 * indicates that this is a null <code>Envelope</code>.
 *
 * @return {number} the minimum x-coordinate.
 */
jsts.geom.Envelope.prototype.getMinX = function() {
  return this.minx;
};


/**
 * Returns the <code>Envelope</code>s maximum x-value. min x > max x
 * indicates that this is a null <code>Envelope</code>.
 *
 * @return {number} the maximum x-coordinate.
 */
jsts.geom.Envelope.prototype.getMaxX = function() {
  return this.maxx;
};


/**
 * Returns the <code>Envelope</code>s minimum y-value. min y > max y
 * indicates that this is a null <code>Envelope</code>.
 *
 * @return {number} the minimum y-coordinate.
 */
jsts.geom.Envelope.prototype.getMinY = function() {
  return this.miny;
};


/**
 * Returns the <code>Envelope</code>s maximum y-value. min y > max y
 * indicates that this is a null <code>Envelope</code>.
 *
 * @return {number} the maximum y-coordinate.
 */
jsts.geom.Envelope.prototype.getMaxY = function() {
  return this.maxy;
};


/**
 * Gets the area of this envelope.
 *
 * @return {number} the area of the envelope, 0.0 if the envelope is null.
 */
jsts.geom.Envelope.prototype.getArea = function() {
  return this.getWidth() * this.getHeight();
};


/**
 * Enlarges this <code>Envelope</code>
 *
 * Will call appropriate expandToInclude* depending on arguments.
 */
jsts.geom.Envelope.prototype.expandToInclude = function() {
  if (arguments[0] instanceof jsts.geom.Coordinate) {
    this.expandToIncludeCoordinate(arguments[0]);
  } else if (arguments[0] instanceof jsts.geom.Envelope) {
    this.expandToIncludeEnvelope(arguments[0]);
  } else {
    this.expandToIncludeValues(arguments[0], arguments[1]);
  }
};


/**
 * Enlarges this <code>Envelope</code> so that it contains the given
 * {@link Coordinate}. Has no effect if the point is already on or within the
 * envelope.
 *
 * @param {jsts.geom.Coordinate}
 *          p the Coordinate to expand to include.
 */
jsts.geom.Envelope.prototype.expandToIncludeCoordinate = function(p) {
  this.expandToIncludeValues(p.x, p.y);
};


/**
 * Enlarges this <code>Envelope</code> so that it contains the given point.
 * Has no effect if the point is already on or within the envelope.
 *
 * @param {number}
 *          x the value to lower the minimum x to or to raise the maximum x to.
 * @param {number}
 *          y the value to lower the minimum y to or to raise the maximum y to.
 */
jsts.geom.Envelope.prototype.expandToIncludeValues = function(x, y) {
  if (this.isNull()) {
    this.minx = x;
    this.maxx = x;
    this.miny = y;
    this.maxy = y;
  } else {
    if (x < this.minx) {
      this.minx = x;
    }
    if (x > this.maxx) {
      this.maxx = x;
    }
    if (y < this.miny) {
      this.miny = y;
    }
    if (y > this.maxy) {
      this.maxy = y;
    }
  }
};


/**
 * Enlarges this <code>Envelope</code> so that it contains the
 * <code>other</code> Envelope. Has no effect if <code>other</code> is
 * wholly on or within the envelope.
 *
 * @param {jsts.geom.Envelope}
 *          other the <code>Envelope</code> to expand to include.
 */
jsts.geom.Envelope.prototype.expandToIncludeEnvelope = function(other) {
  if (other.isNull()) {
    return;
  }
  if (this.isNull()) {
    this.minx = other.getMinX();
    this.maxx = other.getMaxX();
    this.miny = other.getMinY();
    this.maxy = other.getMaxY();
  } else {
    if (other.minx < this.minx) {
      this.minx = other.minx;
    }
    if (other.maxx > this.maxx) {
      this.maxx = other.maxx;
    }
    if (other.miny < this.miny) {
      this.miny = other.miny;
    }
    if (other.maxy > this.maxy) {
      this.maxy = other.maxy;
    }
  }
};


/**
 * Enlarges this <code>Envelope</code>
 *
 * Will call appropriate expandBy* depending on arguments.
 */
jsts.geom.Envelope.prototype.expandBy = function() {
  if (arguments.length === 1) {
    this.expandByDistance(arguments[0]);
  } else {
    this.expandByDistances(arguments[0], arguments[1]);
  }
};


/**
 * Expands this envelope by a given distance in all directions. Both positive
 * and negative distances are supported.
 *
 * @param {number}
 *          distance the distance to expand the envelope.
 */
jsts.geom.Envelope.prototype.expandByDistance = function(distance) {
  this.expandByDistances(distance, distance);
};


/**
 * Expands this envelope by a given distance in all directions. Both positive
 * and negative distances are supported.
 *
 * @param {number}
 *          deltaX the distance to expand the envelope along the the X axis.
 * @param {number}
 *          deltaY the distance to expand the envelope along the the Y axis.
 */
jsts.geom.Envelope.prototype.expandByDistances = function(deltaX, deltaY) {
  if (this.isNull()) {
    return;
  }

  this.minx -= deltaX;
  this.maxx += deltaX;
  this.miny -= deltaY;
  this.maxy += deltaY;

  // check for envelope disappearing
  if (this.minx > this.maxx || this.miny > this.maxy) {
    this.setToNull();
  }
};


/**
 * Translates this envelope by given amounts in the X and Y direction.
 *
 * @param {number}
 *          transX the amount to translate along the X axis.
 * @param {number}
 *          transY the amount to translate along the Y axis.
 */
jsts.geom.Envelope.prototype.translate = function(transX, transY) {
  if (this.isNull()) {
    return;
  }
  this.init(this.minx + transX, this.maxx + transX, this.miny + transY,
      this.maxy + transY);
};


/**
 * Computes the coordinate of the centre of this envelope (as long as it is
 * non-null
 *
 * @return {jsts.geom.Coordinate} the centre coordinate of this envelope <code>null</code>
 *         if the envelope is null.
 */
jsts.geom.Envelope.prototype.centre = function() {
  if (this.isNull()) {
    return null;
  }
  return new jsts.geom.Coordinate((this.minx + this.maxx) / 2.0,
      (this.miny + this.maxy) / 2.0);
};


/**
 * Computes the intersection of two {@link Envelopes}
 *
 * @param {jsts.geom.Envelope}
 *          env the envelope to intersect with.
 * @return {jsts.geom.Envelope} a new Envelope representing the intersection of
 *         the envelopes (this will be the null envelope if either argument is
 *         null, or they do not intersect.
 */
jsts.geom.Envelope.prototype.intersection = function(env) {
  if (this.isNull() || env.isNull() || !this.intersects(env)) {
    return new jsts.geom.Envelope();
  }

  var intMinX = this.minx > env.minx ? this.minx : env.minx;
  var intMinY = this.miny > env.miny ? this.miny : env.miny;
  var intMaxX = this.maxx < env.maxx ? this.maxx : env.maxx;
  var intMaxY = this.maxy < env.maxy ? this.maxy : env.maxy;

  return new jsts.geom.Envelope(intMinX, intMaxX, intMinY, intMaxY);
};


/**
 * Check if the region defined by input overlaps (intersects) the region of this
 * <code>Envelope</code>.
 *
 * Will call appropriate intersects* depending on arguments.
 *
 * @return {boolean} <code>true</code> if an overlap is found.
 */
jsts.geom.Envelope.prototype.intersects = function() {
  if (arguments[0] instanceof jsts.geom.Envelope) {
    return this.intersectsEnvelope(arguments[0]);
  } else if (arguments[0] instanceof jsts.geom.Coordinate) {
    return this.intersectsCoordinate(arguments[0]);
  } else {
    return this.intersectsValues(arguments[0], arguments[1]);
  }
};


/**
 * Check if the region defined by <code>other</code> overlaps (intersects) the
 * region of this <code>Envelope</code>.
 *
 * @param {jsts.geom.Envelope}
 *          other the <code>Envelope</code> which this <code>Envelope</code>
 *          is being checked for overlapping.
 * @return {boolean} <code>true</code> if the <code>Envelope</code>s
 *         overlap.
 */
jsts.geom.Envelope.prototype.intersectsEnvelope = function(other) {
  if (this.isNull() || other.isNull()) {
    return false;
  }

  var result = !(other.minx > this.maxx || other.maxx < this.minx ||
      other.miny > this.maxy || other.maxy < this.miny);
  return result;
};


/**
 * Check if the point <code>p</code> overlaps (lies inside) the region of this
 * <code>Envelope</code>.
 *
 * @param {jsts.geom.Coordinate}
 *          p the <code>Coordinate</code> to be tested.
 * @return {boolean} <code>true</code> if the point overlaps this
 *         <code>Envelope.</code>
 */
jsts.geom.Envelope.prototype.intersectsCoordinate = function(p) {
  return this.intersectsValues(p.x, p.y);
};


/**
 * Check if the point <code>(x, y)</code> overlaps (lies inside) the region of
 * this <code>Envelope</code>.
 *
 * @param {number}
 *          x the x-ordinate of the point.
 * @param {number}
 *          y the y-ordinate of the point.
 * @return {boolean} <code>true</code> if the point overlaps this
 *         <code>Envelope.</code>
 */
jsts.geom.Envelope.prototype.intersectsValues = function(x, y) {
  if (this.isNull()) {
    return false;
  }

  return !(x > this.maxx || x < this.minx || y > this.maxy || y < this.miny);
};


/**
 * Tests if the input lies wholely inside this <code>Envelope</code>
 * (inclusive of the boundary).
 *
 * Will call appropriate contains* depending on arguments.
 *
 * @return {boolean} true if input is contained in this <code>Envelope.</code>
 */
jsts.geom.Envelope.prototype.contains = function() {
  if (arguments[0] instanceof jsts.geom.Envelope) {
    return this.containsEnvelope(arguments[0]);
  } else if (arguments[0] instanceof jsts.geom.Coordinate) {
    return this.containsCoordinate(arguments[0]);
  } else {
    return this.containsValues(arguments[0], arguments[1]);
  }
};


/**
 * Tests if the <code>Envelope other</code> lies wholely inside this
 * <code>Envelope</code> (inclusive of the boundary).
 * <p>
 * Note that this is <b>not</b> the same definition as the SFS
 * <tt>contains</tt>, which would exclude the envelope boundary.
 *
 * @param {jsts.geom.Envelope}
 *          other the <code>Envelope</code> to check.
 * @return {boolean} true if <code>other</code> is contained in this
 *         <code>Envelope.</code>
 *
 * @see covers(Envelope)
 */
jsts.geom.Envelope.prototype.containsEnvelope = function(other) {
  return this.coversEnvelope(other);
};


/**
 * Tests if the given point lies in or on the envelope.
 * <p>
 * Note that this is <b>not</b> the same definition as the SFS
 * <tt>contains</tt>, which would exclude the envelope boundary.
 *
 * @param {jsts.geom.Coordinate}
 *          p the point which this <code>Envelope</code> is being checked for
 *          containing.
 * @return {boolean} <code>true</code> if the point lies in the interior or on
 *         the boundary of this <code>Envelope</code>.
 *
 * @see covers(Coordinate)
 */
jsts.geom.Envelope.prototype.containsCoordinate = function(p) {
  return this.coversCoordinate(p);
};


/**
 * Tests if the given point lies in or on the envelope.
 * <p>
 * Note that this is <b>not</b> the same definition as the SFS
 * <tt>contains</tt>, which would exclude the envelope boundary.
 *
 * @param {number}
 *          x the x-coordinate of the point which this <code>Envelope</code>
 *          is being checked for containing.
 * @param {number}
 *          y the y-coordinate of the point which this <code>Envelope</code>
 *          is being checked for containing.
 * @return {boolean} <code>true</code> if <code>(x, y)</code> lies in the
 *         interior or on the boundary of this <code>Envelope</code>.
 *
 * @see covers(double, double)
 */
jsts.geom.Envelope.prototype.containsValues = function(x, y) {
  return this.coversValues(x, y);
};


/**
 * Tests if the given point lies in or on the envelope.
 *
 * Will call appropriate contains* depending on arguments.
 */
jsts.geom.Envelope.prototype.covers = function() {
  if (arguments[0] instanceof jsts.geom.Envelope) {
    this.coversEnvelope(arguments[0]);
  } else if (arguments[0] instanceof jsts.geom.Coordinate) {
    this.coversCoordinate(arguments[0]);
  } else {
    this.coversValues(arguments[0], arguments[1]);
  }
};


/**
 * Tests if the given point lies in or on the envelope.
 *
 * @param {number}
 *          x the x-coordinate of the point which this <code>Envelope</code>
 *          is being checked for containing.
 * @param {number}
 *          y the y-coordinate of the point which this <code>Envelope</code>
 *          is being checked for containing.
 * @return {boolean} <code>true</code> if <code>(x, y)</code> lies in the
 *         interior or on the boundary of this <code>Envelope</code>.
 */
jsts.geom.Envelope.prototype.coversValues = function(x, y) {
  if (this.isNull()) {
    return false;
  }
  return x >= this.minx && x <= this.maxx && y >= this.miny && y <= this.maxy;
};


/**
 * Tests if the given point lies in or on the envelope.
 *
 * @param {jsts.geom.Coordinate}
 *          p the point which this <code>Envelope</code> is being checked for
 *          containing.
 * @return {boolean} <code>true</code> if the point lies in the interior or on
 *         the boundary of this <code>Envelope</code>.
 */
jsts.geom.Envelope.prototype.coversCoordinate = function(p) {
  return this.coversValues(p.x, p.y);
};


/**
 * Tests if the <code>Envelope other</code> lies wholely inside this
 * <code>Envelope</code> (inclusive of the boundary).
 *
 * @param {jsts.geom.Envelope}
 *          other the <code>Envelope</code> to check.
 * @return {boolean} true if this <code>Envelope</code> covers the
 *         <code>other.</code>
 */
jsts.geom.Envelope.prototype.coversEnvelope = function(other) {
  if (this.isNull() || other.isNull()) {
    return false;
  }
  return other.minx >= this.minx && other.maxx <= this.maxx &&
      other.miny >= this.miny && other.maxy <= this.maxy;
};


/**
 * Computes the distance between this and another <code>Envelope</code>.
 *
 * @param {jsts.geom.Envelope}
 *          env The <code>Envelope</code> to test this <code>Envelope</code>
 *          against.
 * @return {number} The distance between overlapping Envelopes is 0. Otherwise,
 *         the distance is the Euclidean distance between the closest points.
 */
jsts.geom.Envelope.prototype.distance = function(env) {
  if (this.intersects(env)) {
    return 0;
  }
  var dx = 0.0;
  if (this.maxx < env.minx) {
    dx = env.minx - this.maxx;
  }
  if (this.minx > env.maxx) {
    dx = this.minx - env.maxx;
  }

  var dy = 0.0;
  if (this.maxy < env.miny) {
    dy = env.miny - this.maxy;
  }
  if (this.miny > env.maxy) {
    dy = this.miny - env.maxy;
  }

  // if either is zero, the envelopes overlap either vertically or horizontally
  if (dx === 0.0) {
    return dy;
  }
  if (dy === 0.0) {
    return dx;
  }
  return Math.sqrt(dx * dx + dy * dy);
};


/**
 * @param {jsts.geom.Envelope}
 *          other the <code>Envelope</code> to check against.
 * @return {boolean} true if envelopes are equal.
 */
jsts.geom.Envelope.prototype.equals = function(other) {
  if (this.isNull()) {
    return other.isNull();
  }
  return this.maxx === other.maxx && this.maxy === other.maxy &&
      this.minx === other.minx && this.miny === other.miny;
};


/**
 * @return {string} String representation of this <code>Envelope.</code>
 */
jsts.geom.Envelope.prototype.toString = function() {
  return 'Env[' + this.minx + ' : ' + this.maxx + ', ' + this.miny + ' : ' +
      this.maxy + ']';
};


/**
 * Test the point q to see whether it intersects the Envelope defined by p1-p2
 *
 * NOTE: calls intersectsEnvelope if four arguments are given to simulate
 * overloaded function
 *
 * @param {jsts.geom.Coordinate}
 *          p1 one extremal point of the envelope.
 * @param {jsts.geom.Coordinate}
 *          p2 another extremal point of the envelope.
 * @param {jsts.geom.Coordinate}
 *          q the point to test for intersection.
 * @return {boolean} <code>true</code> if q intersects the envelope p1-p2.
 */
jsts.geom.Envelope.intersects = function(p1, p2, q) {
  if (arguments.length === 4) {
    return jsts.geom.Envelope.intersectsEnvelope(arguments[0], arguments[1],
        arguments[2], arguments[3]);
  }

  var xc1 = p1.x < p2.x ? p1.x : p2.x;
  var xc2 = p1.x > p2.x ? p1.x : p2.x;
  var yc1 = p1.y < p2.y ? p1.y : p2.y;
  var yc2 = p1.y > p2.y ? p1.y : p2.y;

  if (((q.x >= xc1) && (q.x <= xc2)) && ((q.y >= yc1) && (q.y <= yc2))) {
    return true;
  }
  return false;
};


/**
 * Test the envelope defined by p1-p2 for intersection with the envelope defined
 * by q1-q2
 *
 * @param {jsts.geom.Coordinate}
 *          p1 one extremal point of the envelope P.
 * @param {jsts.geom.Coordinate}
 *          p2 another extremal point of the envelope P.
 * @param {jsts.geom.Coordinate}
 *          q1 one extremal point of the envelope Q.
 * @param {jsts.geom.Coordinate}
 *          q2 another extremal point of the envelope Q.
 * @return {boolean} <code>true</code> if Q intersects P.
 */
jsts.geom.Envelope.intersectsEnvelope = function(p1, p2, q1, q2) {
  var minq = Math.min(q1.x, q2.x);
  var maxq = Math.max(q1.x, q2.x);
  var minp = Math.min(p1.x, p2.x);
  var maxp = Math.max(p1.x, p2.x);

  if (minp > maxq) {
    return false;
  }
  if (maxp < minq) {
    return false;
  }

  minq = Math.min(q1.y, q2.y);
  maxq = Math.max(q1.y, q2.y);
  minp = Math.min(p1.y, p2.y);
  maxp = Math.max(p1.y, p2.y);

  if (minp > maxq) {
    return false;
  }
  if (maxp < minq) {
    return false;
  }
  return true;
};


/**
 * @return {jsts.geom.Envelope} A new instance copied from this.
 */
jsts.geom.Envelope.prototype.clone = function() {
  return new jsts.geom.Envelope(this.minx, this.maxx, this.miny, this.maxy);
};



//  JSTS NODER

  /**
   * Computes all intersections between segments in a set of
   * {@link SegmentString}s. Intersections found are represented as
   * {@link SegmentNode}s and added to the {@link SegmentString}s in which
   * they occur. As a final step in the noding a new set of segment strings
   * split at the nodes may be returned.
   *
   * @interface
   */
  jsts.noding.Noder = function() {

  };


  /**
   * Computes the noding for a collection of {@link SegmentString}s. Some
   * Noders may add all these nodes to the input SegmentStrings; others may only
   * add some or none at all.
   *
   * @param {Array}
   *          segStrings a collection of {@link SegmentString}s to node.
   */
  jsts.noding.Noder.prototype.computeNodes = jsts.abstractFunc;

  /**
   * Returns a {@link Collection} of fully noded {@link SegmentString}s. The
   * SegmentStrings have the same context as their parent.
   *
   * @return {Array} a Collection of SegmentStrings.
   */
  jsts.noding.Noder.prototype.getNodedSubstrings = jsts.abstractFunc;






//    JSTS ScaledNoder

/* Copyright (c) 2011 by The Authors.
 * Published under the LGPL 2.1 license.
 * See /license-notice.txt for the full text of the license notice.
 * See /license.txt for the full text of the license.
 */

/**
 * Port source: /jts/jts/java/src/com/vividsolutions/jts/noding/ScaledNoder.java
 * Revision: 478
 */

jsts.noding.ScaledNoder = function(noder, scaleFactor, offsetX, offsetY) {
  this.offsetX = offsetX ? offsetX : 0;
  this.offsetY = offsetY ? offsetY : 0;

  this.noder = noder;
  this.scaleFactor = scaleFactor;

  // no need to scale if input precision is already integral
  this.isScaled = !this.isIntegerPrecision();
};

jsts.noding.ScaledNoder.prototype = new jsts.noding.Noder();
jsts.noding.ScaledNoder.constructor = jsts.noding.ScaledNoder;

jsts.noding.ScaledNoder.prototype.noder = null;
jsts.noding.ScaledNoder.prototype.scaleFactor = undefined;
jsts.noding.ScaledNoder.prototype.offsetX = undefined;
jsts.noding.ScaledNoder.prototype.offsetY = undefined;
jsts.noding.ScaledNoder.prototype.isScaled = false;

jsts.noding.ScaledNoder.prototype.isIntegerPrecision = function() {
  return this.scaleFactor === 1.0;
};

jsts.noding.ScaledNoder.prototype.getNodedSubstrings = function() {
  var splitSS = this.noder.getNodedSubstrings();
  if (this.isScaled)
    this.rescale(splitSS);
  return splitSS;
};

jsts.noding.ScaledNoder.prototype.computeNodes = function(inputSegStrings) {
  var intSegStrings = inputSegStrings;
  if (this.isScaled)
    intSegStrings = this.scale(inputSegStrings);
  this.noder.computeNodes(intSegStrings);
};

/**
 * @private
 */
jsts.noding.ScaledNoder.prototype.scale = function(segStrings) {
  if (segStrings instanceof Array) {
    return this.scale2(segStrings);
  }

  var transformed = new javascript.util.ArrayList();
  for (var i = segStrings.iterator(); i.hasNext();) {
    var ss = i.next();
    transformed.add(new jsts.noding.NodedSegmentString(this.scale(ss
        .getCoordinates()), ss.getData()));
  }

  return transformed;
};

/**
 * @private
 */
jsts.noding.ScaledNoder.prototype.scale2 = function(pts) {
  var roundPts = [];
  for (var i = 0; i < pts.length; i++) {
    roundPts[i] = new jsts.geom.Coordinate(Math
        .round((pts[i].x - this.offsetX) * this.scaleFactor), Math
        .round((pts[i].y - this.offsetY) * this.scaleFactor));
  }
  var roundPtsNoDup = jsts.geom.CoordinateArrays.removeRepeatedPoints(roundPts);
  return roundPtsNoDup;
};

/**
 * @private
 */
jsts.noding.ScaledNoder.prototype.rescale = function(segStrings) {
  if (segStrings instanceof Array) {
    this.rescale2(segStrings);
    return;
  }

  for (var i = segStrings.iterator(); i.hasNext();) {
    var ss = i.next();
    this.rescale(ss.getCoordinates());
  }
};

/**
 * @private
 */
jsts.noding.ScaledNoder.prototype.rescale2 = function(pts) {
  for (var i = 0; i < pts.length; i++) {
    pts[i].x = pts[i].x / this.scaleFactor + this.offsetX;
    pts[i].y = pts[i].y / this.scaleFactor + this.offsetY;
  }
};




//    JSTS LineIntersector

/* Copyright (c) 2011 by The Authors.
 * Published under the LGPL 2.1 license.
 * See /license-notice.txt for the full text of the license notice.
 * See /license.txt for the full text of the license.
 */



/**
 * A LineIntersector is an algorithm that can both test whether two line
 * segments intersect and compute the intersection point if they do. The
 * intersection point may be computed in a precise or non-precise manner.
 * Computing it precisely involves rounding it to an integer. (This assumes that
 * the input coordinates have been made precise by scaling them to an integer
 * grid.)
 *
 * @constructor
 */
jsts.algorithm.LineIntersector = function() {
  this.inputLines = [[], []];
  this.intPt = [null, null];
  // alias the intersection points for ease of reference
  this.pa = this.intPt[0];
  this.pb = this.intPt[1];
  this.result = jsts.algorithm.LineIntersector.NO_INTERSECTION;
};


/**
 * Indicates that line segments do not intersect
 *
 * @type {int}
 */
jsts.algorithm.LineIntersector.NO_INTERSECTION = 0;


/**
 * Indicates that line segments intersect in a single point
 *
 * @type {int}
 */
jsts.algorithm.LineIntersector.POINT_INTERSECTION = 1;


/**
 * Indicates that line segments intersect in a line segment
 *
 * @type {int}
 */
jsts.algorithm.LineIntersector.COLLINEAR_INTERSECTION = 2;


/**
 * Force computed intersection to be rounded to a given precision model. No
 * getter is provided, because the precision model is not required to be
 * specified.
 *
 * @param precisionModel
 */
jsts.algorithm.LineIntersector.prototype.setPrecisionModel = function(
    precisionModel) {
  this.precisionModel = precisionModel;
};


/**
 * Gets an endpoint of an input segment.
 *
 * @param segmentIndex
 *          the index of the input segment (0 or 1).
 * @param ptIndex
 *          the index of the endpoint (0 or 1).
 * @return the specified endpoint.
 */
jsts.algorithm.LineIntersector.prototype.getEndpoint = function(segmentIndex,
    ptIndex) {
  return this.inputLines[segmentIndex][ptIndex];
};


/**
 * Computes the "edge distance" of an intersection point p along a segment. The
 * edge distance is a metric of the point along the edge. The metric used is a
 * robust and easy to compute metric function. It is <b>not</b> equivalent to
 * the usual Euclidean metric. It relies on the fact that either the x or the y
 * ordinates of the points in the edge are unique, depending on whether the edge
 * is longer in the horizontal or vertical direction.
 * <p>
 * NOTE: This function may produce incorrect distances for inputs where p is not
 * precisely on p1-p2 (E.g. p = (139,9) p1 = (139,10), p2 = (280,1) produces
 * distanct 0.0, which is incorrect.
 * <p>
 * My hypothesis is that the function is safe to use for points which are the
 * result of <b>rounding</b> points which lie on the line, but not safe to use
 * for <b>truncated</b> points.
 *
 * @param {Coordinate}
 *          p
 * @param {Coordinate}
 *          p0
 * @param {Coordinate}
 *          p1
 * @return {double}
 */
jsts.algorithm.LineIntersector.computeEdgeDistance = function(p, p0, p1) {
  var dx = Math.abs(p1.x - p0.x);
  var dy = Math.abs(p1.y - p0.y);

  var dist = -1.0; // sentinel value
  if (p.equals(p0)) {
    dist = 0.0;
  } else if (p.equals(p1)) {
    if (dx > dy) {
      dist = dx;
    } else {
      dist = dy;
    }
  } else {
    var pdx = Math.abs(p.x - p0.x);
    var pdy = Math.abs(p.y - p0.y);
    if (dx > dy) {
      dist = pdx;
    } else {
      dist = pdy;
    }
    // <FIX>
    // hack to ensure that non-endpoints always have a non-zero distance
    if (dist === 0.0 && !p.equals(p0)) {
      dist = Math.max(pdx, pdy);
    }
  }
  if (dist === 0.0 && !p.equals(p0)) {
    throw new jsts.error.IllegalArgumentError('Bad distance calculation');
  }
  return dist;
};


/**
 * This function is non-robust, since it may compute the square of large
 * numbers. Currently not sure how to improve this.
 *
 * @param {Coordinate}
 *          p
 * @param {Coordinate}
 *          p0
 * @param {Coordinate}
 *          p1
 * @return {double}
 */
jsts.algorithm.LineIntersector.nonRobustComputeEdgeDistance = function(p, p1,
    p2) {
  var dx = p.x - p1.x;
  var dy = p.y - p1.y;
  var dist = Math.sqrt(dx * dx + dy * dy); // dummy value
  if (!(dist === 0.0 && !p.equals(p1))) {
    throw new jsts.error.IllegalArgumentError('Invalid distance calculation');
  }
  return dist;
};


/**
 * @protected
 * @type {int}
 */
jsts.algorithm.LineIntersector.prototype.result = null;


/**
 * @protected
 * @type {Coordinate[][] }
 */
jsts.algorithm.LineIntersector.prototype.inputLines = null;


/**
 * @protected
 * @type {Coordinate[]}
 */
jsts.algorithm.LineIntersector.prototype.intPt = null;


/**
 * The indexes of the endpoints of the intersection lines, in order along the
 * corresponding line
 */
/**
 * @protected
 * @type {int[][]}
 */
jsts.algorithm.LineIntersector.prototype.intLineIndex = null;


/**
 * @protected
 * @type {boolean}
 */
jsts.algorithm.LineIntersector.prototype._isProper = null;


/**
 * @protected
 * @type {Coordinate}
 */
jsts.algorithm.LineIntersector.prototype.pa = null;


/**
 * @protected
 * @type {Coordinate}
 */
jsts.algorithm.LineIntersector.prototype.pb = null;


/**
 * @protected
 * @type {PrecisionModel}
 */
jsts.algorithm.LineIntersector.prototype.precisionModel = null;


/**
 * Compute the intersection of a point p and the line p1-p2. This function
 * computes the boolean value of the hasIntersection test. The actual value of
 * the intersection (if there is one) is equal to the value of <code>p</code>.
 *
 * @param {Coordinate}
 *          p
 * @param {Coordinate}
 *          p1
 * @param {Coordinate}
 *          p2
 */
jsts.algorithm.LineIntersector.prototype.computeIntersection = function(p, p1, p2) {
  throw new jsts.error.AbstractMethodInvocationError();
};


/**
 * @return {boolean}
 * @protected
 */
jsts.algorithm.LineIntersector.prototype.isCollinear = function() {
  return this.result === jsts.algorithm.LineIntersector.COLLINEAR_INTERSECTION;
};


/**
 * Computes the intersection of the lines p1-p2 and p3-p4. This function
 * computes both the boolean value of the hasIntersection test and the
 * (approximate) value of the intersection point itself (if there is one).
 *
 * @param {Coordinate}
 *          p1
 * @param {Coordinate}
 *          p2
 * @param {Coordinate}
 *          p3
 * @param {Coordinate}
 *          p4
 */
jsts.algorithm.LineIntersector.prototype.computeIntersection = function(p1, p2,
    p3, p4) {
  this.inputLines[0][0] = p1;
  this.inputLines[0][1] = p2;
  this.inputLines[1][0] = p3;
  this.inputLines[1][1] = p4;
  this.result = this.computeIntersect(p1, p2, p3, p4);
};


/**
 * @param {Coordinate}
 *          p1
 * @param {Coordinate}
 *          p2
 * @param {Coordinate}
 *          q1
 * @param {Coordinate}
 *          q2
 * @return {int}
 * @protected
 */
jsts.algorithm.LineIntersector.prototype.computeIntersect = function(p1, p2,
    q1, q2) {
  throw new jsts.error.AbstractMethodInvocationError();
};


/**
 * @return {boolean}
 * @protected
 */
jsts.algorithm.LineIntersector.prototype.isEndPoint = function() {
  return this.hasIntersection() && !this._isProper;
};


/**
 * Tests whether the input geometries intersect.
 *
 * @return {boolean} true if the input geometries intersect.
 */
jsts.algorithm.LineIntersector.prototype.hasIntersection = function() {
  return this.result !== jsts.algorithm.LineIntersector.NO_INTERSECTION;
};


/**
 * Returns the number of intersection points found. This will be either 0, 1 or
 * 2.
 *
 * @return {int}
 */
jsts.algorithm.LineIntersector.prototype.getIntersectionNum = function() {
  return this.result;
};


/**
 * Returns the intIndex'th intersection point
 *
 * @param {int}
 *          intIndex is 0 or 1.
 *
 * @return {Coordinate} the intIndex'th intersection point.
 */
jsts.algorithm.LineIntersector.prototype.getIntersection = function(intIndex) {
  return this.intPt[intIndex];
};


/**
 * @protected
 */
jsts.algorithm.LineIntersector.prototype.computeIntLineIndex = function() {
  if (this.intLineIndex === null) {
    this.intLineIndex = [[], []];
    this.computeIntLineIndex(0);
    this.computeIntLineIndex(1);
  }
};


/**
 * Test whether a point is a intersection point of two line segments. Note that
 * if the intersection is a line segment, this method only tests for equality
 * with the endpoints of the intersection segment. It does <b>not</b> return
 * true if the input point is internal to the intersection segment.
 *
 * @param {Coordinate}
 *          pt
 * @return {boolean} true if the input point is one of the intersection points.
 */
jsts.algorithm.LineIntersector.prototype.isIntersection = function(pt) {
  var i;
  for (i = 0; i < this.result; i++) {
    if (this.intPt[i].equals2D(pt)) {
      return true;
    }
  }
  return false;
};


/**
 * Tests whether either intersection point is an interior point of one of the
 * input segments.
 *
 * @return {boolean} <code>true</code> if either intersection point is in the
 *         interior of one of the input segments.
 */
jsts.algorithm.LineIntersector.prototype.isInteriorIntersection = function() {
  if (arguments.length === 1) {
    return this.isInteriorIntersection2.apply(this, arguments);
  }

  if (this.isInteriorIntersection(0)) {
    return true;
  }
  if (this.isInteriorIntersection(1)) {
    return true;
  }
  return false;
};


/**
 * Tests whether either intersection point is an interior point of the specified
 * input segment.
 *
 * @param {[]} inputLineIndex
 * @return {boolean} <code>true</code> if either intersection point is in the
 *         interior of the input segment.
 */
jsts.algorithm.LineIntersector.prototype.isInteriorIntersection2 = function(
    inputLineIndex) {
  var i;
  for (i = 0; i < this.result; i++) {
    if (!(this.intPt[i].equals2D(this.inputLines[inputLineIndex][0]) || this.intPt[i]
        .equals2D(this.inputLines[inputLineIndex][1]))) {
      return true;
    }
  }
  return false;
};


/**
 * Tests whether an intersection is proper. <br>
 * The intersection between two line segments is considered proper if they
 * intersect in a single point in the interior of both segments (e.g. the
 * intersection is a single point and is not equal to any of the endpoints).
 * <p>
 * The intersection between a point and a line segment is considered proper if
 * the point lies in the interior of the segment (e.g. is not equal to either of
 * the endpoints).
 *
 * @return {boolean} true if the intersection is proper.
 */
jsts.algorithm.LineIntersector.prototype.isProper = function() {
  return this.hasIntersection() && this._isProper;
};


/**
 * Computes the intIndex'th intersection point in the direction of a specified
 * input line segment
 *
 * @param {int}
 *          segmentIndex is 0 or 1.
 * @param {int}
 *          intIndex is 0 or 1.
 *
 * @return {Coordinate} the intIndex'th intersection point in the direction of
 *         the specified input line segment.
 */
jsts.algorithm.LineIntersector.prototype.getIntersectionAlongSegment = function(
    segmentIndex, intIndex) {
  // lazily compute int line array
  this.computeIntLineIndex();
  return this.intPt[intLineIndex[segmentIndex][intIndex]];
};


/**
 * Computes the index of the intIndex'th intersection point in the direction of
 * a specified input line segment
 *
 * @param {int}
 *          segmentIndex is 0 or 1.
 * @param {int}
 *          intIndex is 0 or 1.
 *
 * @return {int} the index of the intersection point along the segment (0 or 1).
 */
jsts.algorithm.LineIntersector.prototype.getIndexAlongSegment = function(
    segmentIndex, intIndex) {
  this.computeIntLineIndex();
  return this.intLineIndex[segmentIndex][intIndex];
};


/**
 * @param {int}
 *          segmentIndex
 * @protected
 */
jsts.algorithm.LineIntersector.prototype.computeIntLineIndex = function(
    segmentIndex) {
  var dist0 = this.getEdgeDistance(segmentIndex, 0);
  var dist1 = this.getEdgeDistance(segmentIndex, 1);
  if (dist0 > dist1) {
    this.intLineIndex[segmentIndex][0] = 0;
    this.intLineIndex[segmentIndex][1] = 1;
  } else {
    this.intLineIndex[segmentIndex][0] = 1;
    this.intLineIndex[segmentIndex][1] = 0;
  }
};


/**
 * Computes the "edge distance" of an intersection point along the specified
 * input line segment.
 *
 * @param {int}
 *          segmentIndex is 0 or 1.
 * @param {int}
 *          intIndex is 0 or 1.
 *
 * @return {double} the edge distance of the intersection point.
 */
jsts.algorithm.LineIntersector.prototype.getEdgeDistance = function(
    segmentIndex, intIndex) {
  var dist = jsts.algorithm.LineIntersector.computeEdgeDistance(
      this.intPt[intIndex], this.inputLines[segmentIndex][0],
      this.inputLines[segmentIndex][1]);
  return dist;
};




//    JSTS RobustLineIntersector

/* Copyright (c) 2011 by The Authors.
 * Published under the LGPL 2.1 license.
 * See /license-notice.txt for the full text of the license notice.
 * See /license.txt for the full text of the license.
 */


/**
 * @requires jsts/algorithm/LineIntersector.js
 */


/**
 * A robust version of {@LineIntersector}.
 *
 * @constructor
 * @augments jsts.algorithm.LineIntersector
 */
jsts.algorithm.RobustLineIntersector = function() {
  jsts.algorithm.RobustLineIntersector.prototype.constructor.call(this);
};

jsts.algorithm.RobustLineIntersector.prototype = new jsts.algorithm.LineIntersector();


/**
 * @param {Coordinate}
 *          p
 * @param {Coordinate}
 *          p1
 * @param {Coordinate}
 *          p2
 */
jsts.algorithm.RobustLineIntersector.prototype.computeIntersection = function(
    p, p1, p2) {

  if (arguments.length === 4) {
    jsts.algorithm.LineIntersector.prototype.computeIntersection.apply(this, arguments);
    return;
  }

  this._isProper = false;
  // do between check first, since it is faster than the orientation test
  if (jsts.geom.Envelope.intersects(p1, p2, p)) {
    if ((jsts.algorithm.CGAlgorithms.orientationIndex(p1, p2, p) === 0) &&
        (jsts.algorithm.CGAlgorithms.orientationIndex(p2, p1, p) === 0)) {
      this._isProper = true;
      if (p.equals(p1) || p.equals(p2)) {
        this._isProper = false;
      }
      this.result = jsts.algorithm.LineIntersector.POINT_INTERSECTION;
      return;
    }
  }
  this.result = jsts.algorithm.LineIntersector.NO_INTERSECTION;
};


/**
 * @param {Coordinate}
 *          p1
 * @param {Coordinate}
 *          p2
 * @param {Coordinate}
 *          q1
 * @param {Coordinate}
 *          q2
 * @return {Number}
 * @protected
 */
jsts.algorithm.RobustLineIntersector.prototype.computeIntersect = function(p1,
    p2, q1, q2) {
  this._isProper = false;

  // first try a fast test to see if the envelopes of the lines intersect
  if (!jsts.geom.Envelope.intersects(p1, p2, q1, q2)) {
    return jsts.algorithm.LineIntersector.NO_INTERSECTION;
  }

  // for each endpoint, compute which side of the other segment it lies
  // if both endpoints lie on the same side of the other segment,
  // the segments do not intersect
  var Pq1 = jsts.algorithm.CGAlgorithms.orientationIndex(p1, p2, q1);
  var Pq2 = jsts.algorithm.CGAlgorithms.orientationIndex(p1, p2, q2);

  if ((Pq1 > 0 && Pq2 > 0) || (Pq1 < 0 && Pq2 < 0)) {
    return jsts.algorithm.LineIntersector.NO_INTERSECTION;
  }

  var Qp1 = jsts.algorithm.CGAlgorithms.orientationIndex(q1, q2, p1);
  var Qp2 = jsts.algorithm.CGAlgorithms.orientationIndex(q1, q2, p2);

  if ((Qp1 > 0 && Qp2 > 0) || (Qp1 < 0 && Qp2 < 0)) {
    return jsts.algorithm.LineIntersector.NO_INTERSECTION;
  }

  var collinear = Pq1 === 0 && Pq2 === 0 && Qp1 === 0 && Qp2 === 0;
  if (collinear) {
    return this.computeCollinearIntersection(p1, p2, q1, q2);
  }

  /**
   * At this point we know that there is a single intersection point (since the
   * lines are not collinear).
   */

  /**
   * Check if the intersection is an endpoint. If it is, copy the endpoint as
   * the intersection point. Copying the point rather than computing it ensures
   * the point has the exact value, which is important for robustness. It is
   * sufficient to simply check for an endpoint which is on the other line,
   * since at this point we know that the inputLines must intersect.
   */
  if (Pq1 === 0 || Pq2 === 0 || Qp1 === 0 || Qp2 === 0) {
    this._isProper = false;

    /**
     * Check for two equal endpoints. This is done explicitly rather than by the
     * orientation tests below in order to improve robustness.
     *
     * [An example where the orientation tests fail to be consistent is the
     * following (where the true intersection is at the shared endpoint POINT
     * (19.850257749638203 46.29709338043669)
     *
     * LINESTRING ( 19.850257749638203 46.29709338043669, 20.31970698357233
     * 46.76654261437082 ) and LINESTRING ( -48.51001596420236
     * -22.063180333403878, 19.850257749638203 46.29709338043669 )
     *
     * which used to produce the INCORRECT result: (20.31970698357233,
     * 46.76654261437082, NaN)
     *
     */
    if (p1.equals2D(q1) || p1.equals2D(q2)) {
      this.intPt[0] = p1;
    } else if (p2.equals2D(q1) || p2.equals2D(q2)) {
      this.intPt[0] = p2;
    }

    /**
     * Now check to see if any endpoint lies on the interior of the other
     * segment.
     */
    else if (Pq1 === 0) {
      this.intPt[0] = new jsts.geom.Coordinate(q1);
    } else if (Pq2 === 0) {
      this.intPt[0] = new jsts.geom.Coordinate(q2);
    } else if (Qp1 === 0) {
      this.intPt[0] = new jsts.geom.Coordinate(p1);
    } else if (Qp2 === 0) {
      this.intPt[0] = new jsts.geom.Coordinate(p2);
    }
  } else {
    this._isProper = true;
    this.intPt[0] = this.intersection(p1, p2, q1, q2);
  }
  return jsts.algorithm.LineIntersector.POINT_INTERSECTION;
};


/**
 * @param {Coordinate}
 *          p1
 * @param {Coordinate}
 *          p2
 * @param {Coordinate}
 *          q1
 * @param {Coordinate}
 *          q2
 * @return {Number}
 * @private
 */
jsts.algorithm.RobustLineIntersector.prototype.computeCollinearIntersection = function(
    p1, p2, q1, q2) {
  var p1q1p2 = jsts.geom.Envelope.intersects(p1, p2, q1);
  var p1q2p2 = jsts.geom.Envelope.intersects(p1, p2, q2);
  var q1p1q2 = jsts.geom.Envelope.intersects(q1, q2, p1);
  var q1p2q2 = jsts.geom.Envelope.intersects(q1, q2, p2);

  if (p1q1p2 && p1q2p2) {
    this.intPt[0] = q1;
    this.intPt[1] = q2;
    return jsts.algorithm.LineIntersector.COLLINEAR_INTERSECTION;
  }
  if (q1p1q2 && q1p2q2) {
    this.intPt[0] = p1;
    this.intPt[1] = p2;
    return jsts.algorithm.LineIntersector.COLLINEAR_INTERSECTION;
  }
  if (p1q1p2 && q1p1q2) {
    this.intPt[0] = q1;
    this.intPt[1] = p1;
    return q1.equals(p1) && !p1q2p2 && !q1p2q2 ? jsts.algorithm.LineIntersector.POINT_INTERSECTION
        : jsts.algorithm.LineIntersector.COLLINEAR_INTERSECTION;
  }
  if (p1q1p2 && q1p2q2) {
    this.intPt[0] = q1;
    this.intPt[1] = p2;
    return q1.equals(p2) && !p1q2p2 && !q1p1q2 ? jsts.algorithm.LineIntersector.POINT_INTERSECTION
        : jsts.algorithm.LineIntersector.COLLINEAR_INTERSECTION;
  }
  if (p1q2p2 && q1p1q2) {
    this.intPt[0] = q2;
    this.intPt[1] = p1;
    return q2.equals(p1) && !p1q1p2 && !q1p2q2 ? jsts.algorithm.LineIntersector.POINT_INTERSECTION
        : jsts.algorithm.LineIntersector.COLLINEAR_INTERSECTION;
  }
  if (p1q2p2 && q1p2q2) {
    this.intPt[0] = q2;
    this.intPt[1] = p2;
    return q2.equals(p2) && !p1q1p2 && !q1p1q2 ? jsts.algorithm.LineIntersector.POINT_INTERSECTION
        : jsts.algorithm.LineIntersector.COLLINEAR_INTERSECTION;
  }
  return jsts.algorithm.LineIntersector.NO_INTERSECTION;
};


/**
 * This method computes the actual value of the intersection point. To obtain
 * the maximum precision from the intersection calculation, the coordinates are
 * normalized by subtracting the minimum ordinate values (in absolute value).
 * This has the effect of removing common significant digits from the
 * calculation to maintain more bits of precision.
 *
 * @param {Coordinate}
 *          p1
 * @param {Coordinate}
 *          p2
 * @param {Coordinate}
 *          q1
 * @param {Coordinate}
 *          q2
 * @return {Coordinate}
 * @private
 */
jsts.algorithm.RobustLineIntersector.prototype.intersection = function(p1, p2,
    q1, q2) {
  var intPt = this.intersectionWithNormalization(p1, p2, q1, q2);

  /**
   * Due to rounding it can happen that the computed intersection is outside the
   * envelopes of the input segments. Clearly this is inconsistent. This code
   * checks this condition and forces a more reasonable answer
   *
   * MD - May 4 2005 - This is still a problem. Here is a failure case:
   *
   * LINESTRING (2089426.5233462777 1180182.3877339689, 2085646.6891757075
   * 1195618.7333999649) LINESTRING (1889281.8148903656 1997547.0560044837,
   * 2259977.3672235999 483675.17050843034) int point =
   * (2097408.2633752143,1144595.8008114607)
   *
   * MD - Dec 14 2006 - This does not seem to be a failure case any longer
   */
  if (!this.isInSegmentEnvelopes(intPt)) {
    // System.out.println("Intersection outside segment envelopes: " + intPt);
    // System.out.println("Segments: " + this);
    // compute a safer result
    intPt = jsts.algorithm.CentralEndpointIntersector.getIntersection(p1, p2, q1, q2);
    // System.out.println("Snapped to " + intPt);
  }

  if (this.precisionModel !== null) {
    this.precisionModel.makePrecise(intPt);
  }

  return intPt;
};


/**
 * @param {Coordinate}
 *          p1
 * @param {Coordinate}
 *          p2
 * @param {Coordinate}
 *          q1
 * @param {Coordinate}
 *          q2
 * @return {Coordinate}
 * @private
 */
jsts.algorithm.RobustLineIntersector.prototype.intersectionWithNormalization = function(
    p1, p2, q1, q2) {
  var n1 = new jsts.geom.Coordinate(p1);
  var n2 = new jsts.geom.Coordinate(p2);
  var n3 = new jsts.geom.Coordinate(q1);
  var n4 = new jsts.geom.Coordinate(q2);
  var normPt = new jsts.geom.Coordinate();
  this.normalizeToEnvCentre(n1, n2, n3, n4, normPt);

  var intPt = this.safeHCoordinateIntersection(n1, n2, n3, n4);

  intPt.x += normPt.x;
  intPt.y += normPt.y;

  return intPt;
};


/**
 * Computes a segment intersection using homogeneous coordinates. Round-off
 * error can cause the raw computation to fail, (usually due to the segments
 * being approximately parallel). If this happens, a reasonable approximation is
 * computed instead.
 *
 * @param {Coordinate}
 *          p1 a segment endpoint.
 * @param {Coordinate}
 *          p2 a segment endpoint.
 * @param {Coordinate}
 *          q1 a segment endpoint.
 * @param {Coordinate}
 *          q2 a segment endpoint.
 * @return {Coordinate} the computed intersection point.
 * @private
 */
jsts.algorithm.RobustLineIntersector.prototype.safeHCoordinateIntersection = function(
    p1, p2, q1, q2) {
  var intPt = null;
  try {
    intPt = jsts.algorithm.HCoordinate.intersection(p1, p2, q1, q2);
  } catch (e) {
    if (e instanceof jsts.error.NotRepresentableError) {
      // System.out.println("Not calculable: " + this);
      // compute an approximate result
      intPt = jsts.algorithm.CentralEndpointIntersector.getIntersection(p1, p2,
          q1, q2);
      // System.out.println("Snapped to " + intPt);
    } else {
      throw e;
    }
  }

  return intPt;
};


/**
 * Normalize the supplied coordinates so that their minimum ordinate values lie
 * at the origin. NOTE: this normalization technique appears to cause large
 * errors in the position of the intersection point for some cases.
 *
 * @param {Coordinate}
 *          n1
 * @param {Coordinate}
 *          n2
 * @param {Coordinate}
 *          n3
 * @param {Coordinate}
 *          n4
 * @param {Coordinate}
 *          normPt
 */
jsts.algorithm.RobustLineIntersector.prototype.normalizeToMinimum = function(
    n1, n2, n3, n4, normPt) {
  normPt.x = this.smallestInAbsValue(n1.x, n2.x, n3.x, n4.x);
  normPt.y = this.smallestInAbsValue(n1.y, n2.y, n3.y, n4.y);
  n1.x -= normPt.x;
  n1.y -= normPt.y;
  n2.x -= normPt.x;
  n2.y -= normPt.y;
  n3.x -= normPt.x;
  n3.y -= normPt.y;
  n4.x -= normPt.x;
  n4.y -= normPt.y;
};


/**
 * Normalize the supplied coordinates to so that the midpoint of their
 * intersection envelope lies at the origin.
 *
 * @param {Coordinate}
 *          n00
 * @param {Coordinate}
 *          n01
 * @param {Coordinate}
 *          n10
 * @param {Coordinate}
 *          n11
 * @param {Coordinate}
 *          normPt
 */
jsts.algorithm.RobustLineIntersector.prototype.normalizeToEnvCentre = function(
    n00, n01, n10, n11, normPt) {
  var minX0 = n00.x < n01.x ? n00.x : n01.x;
  var minY0 = n00.y < n01.y ? n00.y : n01.y;
  var maxX0 = n00.x > n01.x ? n00.x : n01.x;
  var maxY0 = n00.y > n01.y ? n00.y : n01.y;

  var minX1 = n10.x < n11.x ? n10.x : n11.x;
  var minY1 = n10.y < n11.y ? n10.y : n11.y;
  var maxX1 = n10.x > n11.x ? n10.x : n11.x;
  var maxY1 = n10.y > n11.y ? n10.y : n11.y;

  var intMinX = minX0 > minX1 ? minX0 : minX1;
  var intMaxX = maxX0 < maxX1 ? maxX0 : maxX1;
  var intMinY = minY0 > minY1 ? minY0 : minY1;
  var intMaxY = maxY0 < maxY1 ? maxY0 : maxY1;

  var intMidX = (intMinX + intMaxX) / 2.0;
  var intMidY = (intMinY + intMaxY) / 2.0;
  normPt.x = intMidX;
  normPt.y = intMidY;

  /*
  // equilavalent code using more modular but slower method
  Envelope env0 = new Envelope(n00, n01);
  Envelope env1 = new Envelope(n10, n11);
  Envelope intEnv = env0.intersection(env1);
  Coordinate intMidPt = intEnv.centre();

  normPt.x = intMidPt.x;
  normPt.y = intMidPt.y;
  */

  n00.x -= normPt.x;
  n00.y -= normPt.y;
  n01.x -= normPt.x;
  n01.y -= normPt.y;
  n10.x -= normPt.x;
  n10.y -= normPt.y;
  n11.x -= normPt.x;
  n11.y -= normPt.y;
};


/**
 * @param {double}
 *          x1
 * @param {double}
 *          x2
 * @param {double}
 *          x3
 * @param {double}
 *          x4
 * @return {double}
 */
jsts.algorithm.RobustLineIntersector.prototype.smallestInAbsValue = function(
    x1, x2, x3, x4) {
  var x = x1;
  var xabs = Math.abs(x);
  if (Math.abs(x2) < xabs) {
    x = x2;
    xabs = Math.abs(x2);
  }
  if (Math.abs(x3) < xabs) {
    x = x3;
    xabs = Math.abs(x3);
  }
  if (Math.abs(x4) < xabs) {
    x = x4;
  }
  return x;
};


/**
 * Test whether a point lies in the envelopes of both input segments. A
 * correctly computed intersection point should return <code>true</code> for
 * this test. Since this test is for debugging purposes only, no attempt is made
 * to optimize the envelope test.
 *
 * @param {Coordinate}
 *          intPt
 * @return {boolean} <code>true</code> if the input point lies within both
 *         input segment envelopes.
 * @private
 */
jsts.algorithm.RobustLineIntersector.prototype.isInSegmentEnvelopes = function(
    intPt) {
  var env0 = new jsts.geom.Envelope(this.inputLines[0][0],
      this.inputLines[0][1]);
  var env1 = new jsts.geom.Envelope(this.inputLines[1][0],
      this.inputLines[1][1]);
  return env0.contains(intPt) && env1.contains(intPt);
};




//    JSTS MCIndexSnapRounder

/* Copyright (c) 2011 by The Authors.
 * Published under the LGPL 2.1 license.
 * See /license-notice.txt for the full text of the license notice.
 * See /license.txt for the full text of the license.
 */

/**
 * Port source: /jts/jts/java/src/com/vividsolutions/jts/noding/snapround/MCIndexSnapRounder.java
 * Revision: 486
 */

/**
 * @requires jsts/algorithm/RobustLineIntersector.js
 * @requires jsts/noding/Noder.js
 */

/**
 * Uses Snap Rounding to compute a rounded, fully noded arrangement from a set
 * of {@link SegmentString}s. Implements the Snap Rounding technique described
 * in papers by Hobby, Guibas & Marimont, and Goodrich et al. Snap Rounding
 * assumes that all vertices lie on a uniform grid; hence the precision model of
 * the input must be fixed precision, and all the input vertices must be rounded
 * to that precision.
 * <p>
 * This implementation uses a monotone chains and a spatial index to speed up
 * the intersection tests.
 * <p>
 * This implementation appears to be fully robust using an integer precision
 * model. It will function with non-integer precision models, but the results
 * are not 100% guaranteed to be correctly noded.
 */
jsts.noding.snapround.MCIndexSnapRounder = function(pm) {
  this.pm = pm;
  this.li = new jsts.algorithm.RobustLineIntersector();
  this.li.setPrecisionModel(pm);
  this.scaleFactor = pm.getScale();
};

jsts.noding.snapround.MCIndexSnapRounder.prototype = new jsts.noding.Noder();
jsts.noding.snapround.MCIndexSnapRounder.constructor = jsts.noding.snapround.MCIndexSnapRounder;


jsts.noding.snapround.MCIndexSnapRounder.prototype.pm = null;
jsts.noding.snapround.MCIndexSnapRounder.prototype.li = null;
jsts.noding.snapround.MCIndexSnapRounder.prototype.scaleFactor = null;
jsts.noding.snapround.MCIndexSnapRounder.prototype.noder = null;
jsts.noding.snapround.MCIndexSnapRounder.prototype.pointSnapper = null;
jsts.noding.snapround.MCIndexSnapRounder.prototype.nodedSegStrings = null;

jsts.noding.snapround.MCIndexSnapRounder.prototype.getNodedSubstrings = function() {
  return jsts.noding.NodedSegmentString
      .getNodedSubstrings(this.nodedSegStrings);
};

jsts.noding.snapround.MCIndexSnapRounder.prototype.computeNodes = function(
    inputSegmentStrings) {
  this.nodedSegStrings = inputSegmentStrings;
  this.noder = new jsts.noding.MCIndexNoder();
  this.pointSnapper = new jsts.noding.snapround.MCIndexPointSnapper(this.noder
      .getIndex());
  this.snapRound(inputSegmentStrings, this.li);
};

/**
 * @private
 */
jsts.noding.snapround.MCIndexSnapRounder.prototype.snapRound = function(
    segStrings, li) {
  var intersections = this.findInteriorIntersections(segStrings, li);
  this.computeIntersectionSnaps(intersections);
  this.computeVertexSnaps(segStrings);
};

/**
 * Computes all interior intersections in the collection of
 * {@link SegmentString}s, and returns their
 *
 * @link Coordinate}s.
 *
 * Does NOT node the segStrings.
 *
 * @return a list of Coordinates for the intersections.
 * @private
 */
jsts.noding.snapround.MCIndexSnapRounder.prototype.findInteriorIntersections = function(
    segStrings, li) {
  var intFinderAdder = new jsts.noding.IntersectionFinderAdder(li);
  this.noder.setSegmentIntersector(intFinderAdder);
  this.noder.computeNodes(segStrings);
  return intFinderAdder.getInteriorIntersections();
};

/**
 * Computes nodes introduced as a result of snapping segments to snap points
 * (hot pixels)
 *
 * @private
 */
jsts.noding.snapround.MCIndexSnapRounder.prototype.computeIntersectionSnaps = function(
    snapPts) {
  for (var it = snapPts.iterator(); it.hasNext();) {
    var snapPt = it.next();
    var hotPixel = new jsts.noding.snapround.HotPixel(snapPt, this.scaleFactor,
        this.li);
    this.pointSnapper.snap(hotPixel);
  }
};

/**
 * Computes nodes introduced as a result of snapping segments to vertices of
 * other segments
 *
 * @param edges
 *          the list of segment strings to snap together.
 */
jsts.noding.snapround.MCIndexSnapRounder.prototype.computeVertexSnaps = function(
    edges) {
  if (edges instanceof jsts.noding.NodedSegmentString) {
    this.computeVertexSnaps2.apply(this, arguments);
    return;
  }

  for (var i0 = edges.iterator(); i0.hasNext();) {
    var edge0 = i0.next();
    this.computeVertexSnaps(edge0);
  }
};

/**
 * Performs a brute-force comparison of every segment in each
 * {@link SegmentString}. This has n^2 performance.
 *
 * @private
 */
jsts.noding.snapround.MCIndexSnapRounder.prototype.computeVertexSnaps2 = function(
    e) {
  var pts0 = e.getCoordinates();
  for (var i = 0; i < pts0.length - 1; i++) {
    var hotPixel = new jsts.noding.snapround.HotPixel(pts0[i],
        this.scaleFactor, this.li);
    var isNodeAdded = this.pointSnapper.snap(hotPixel, e, i);
    // if a node is created for a vertex, that vertex must be noded too
    if (isNodeAdded) {
      e.addIntersection(pts0[i], i);
    }
  }
};




//    JSTS BufferBuilder

/* Copyright (c) 2011 by The Authors.
 * Published under the LGPL 2.1 license.
 * See /license-notice.txt for the full text of the license notice.
 * See /license.txt for the full text of the license.
 */



/**
 * Builds the buffer geometry for a given input geometry and precision model.
 * Allows setting the level of approximation for circular arcs, and the
 * precision model in which to carry out the computation.
 * <p>
 * When computing buffers in floating point double-precision it can happen that
 * the process of iterated noding can fail to converge (terminate). In this case
 * a TopologyException will be thrown. Retrying the computation in a fixed
 * precision can produce more robust results.
 *
 * @param {jsts.operation.buffer.BufferBuilder.BufferParameters}
 *          bufParams
 * @constructor
 */
jsts.operation.buffer.BufferBuilder = function(bufParams) {
  this.bufParams = bufParams;

  this.edgeList = new jsts.geomgraph.EdgeList();
};


/**
 * Compute the change in depth as an edge is crossed from R to L
 *
 * @param {Label}
 *          label
 * @return {Number}
 */
jsts.operation.buffer.BufferBuilder.depthDelta = function(label) {
  var lLoc = label.getLocation(0, jsts.geomgraph.Position.LEFT);
  var rLoc = label.getLocation(0, jsts.geomgraph.Position.RIGHT);
  if (lLoc === jsts.geom.Location.INTERIOR &&
      rLoc === jsts.geom.Location.EXTERIOR)
    return 1;
  else if (lLoc === jsts.geom.Location.EXTERIOR &&
      rLoc === jsts.geom.Location.INTERIOR)
    return -1;
  return 0;
};


/**
 * @type {BufferParameters}
 * @private
 */
jsts.operation.buffer.BufferBuilder.prototype.bufParams = null;


/**
 * @type {PrecisionModel}
 * @private
 */
jsts.operation.buffer.BufferBuilder.prototype.workingPrecisionModel = null;


/**
 * @type {Noder}
 * @private
 */
jsts.operation.buffer.BufferBuilder.prototype.workingNoder = null;


/**
 * @type {GeometryFactory}
 * @private
 */
jsts.operation.buffer.BufferBuilder.prototype.geomFact = null;


/**
 * @type {PlanarGraph}
 * @private
 */
jsts.operation.buffer.BufferBuilder.prototype.graph = null;


/**
 * @type {EdgeList}
 * @private
 */
jsts.operation.buffer.BufferBuilder.prototype.edgeList = null;


/**
 * Sets the precision model to use during the curve computation and noding, if
 * it is different to the precision model of the Geometry. If the precision
 * model is less than the precision of the Geometry precision model, the
 * Geometry must have previously been rounded to that precision.
 *
 * @param pm
 *          the precision model to use.
 */
jsts.operation.buffer.BufferBuilder.prototype.setWorkingPrecisionModel = function(
    pm) {
  this.workingPrecisionModel = pm;
};


/**
 * Sets the {@link Noder} to use during noding. This allows choosing fast but
 * non-robust noding, or slower but robust noding.
 *
 * @param noder
 *          the noder to use.
 */
jsts.operation.buffer.BufferBuilder.prototype.setNoder = function(noder) {
  this.workingNoder = noder;
};

jsts.operation.buffer.BufferBuilder.prototype.buffer = function(g, distance) {
  var precisionModel = this.workingPrecisionModel;
  if (precisionModel === null)
    precisionModel = g.getPrecisionModel();

  // factory must be the same as the one used by the input
  this.geomFact = g.getFactory();

  var curveBuilder = new jsts.operation.buffer.OffsetCurveBuilder(
      precisionModel, this.bufParams);

  var curveSetBuilder = new jsts.operation.buffer.OffsetCurveSetBuilder(g,
      distance, curveBuilder);

  var bufferSegStrList = curveSetBuilder.getCurves();

  // short-circuit test
  if (bufferSegStrList.size() <= 0) {
    return this.createEmptyResultGeometry();
  }

  this.computeNodedEdges(bufferSegStrList, precisionModel);
  this.graph = new jsts.geomgraph.PlanarGraph(
      new jsts.operation.overlay.OverlayNodeFactory());
  this.graph.addEdges(this.edgeList.getEdges());

  var subgraphList = this.createSubgraphs(this.graph);
  var polyBuilder = new jsts.operation.overlay.PolygonBuilder(this.geomFact);
  this.buildSubgraphs(subgraphList, polyBuilder);
  var resultPolyList = polyBuilder.getPolygons();

  // just in case...
  if (resultPolyList.size() <= 0) {
    return this.createEmptyResultGeometry();
  }

  var resultGeom = this.geomFact.buildGeometry(resultPolyList);
  return resultGeom;
};


/**
 * @private
 */
jsts.operation.buffer.BufferBuilder.prototype.getNoder = function(
    precisionModel) {
  if (this.workingNoder !== null)
    return this.workingNoder;

  // otherwise use a fast (but non-robust) noder
  var noder = new jsts.noding.MCIndexNoder();
  var li = new jsts.algorithm.RobustLineIntersector();
  li.setPrecisionModel(precisionModel);
  noder.setSegmentIntersector(new jsts.noding.IntersectionAdder(li));
  return noder;
};


/**
 * @private
 */
jsts.operation.buffer.BufferBuilder.prototype.computeNodedEdges = function(
    bufferSegStrList, precisionModel) {
  var noder = this.getNoder(precisionModel);
  noder.computeNodes(bufferSegStrList);
  var nodedSegStrings = noder.getNodedSubstrings();

  for (var i = nodedSegStrings.iterator(); i.hasNext();) {
    var segStr = i.next();
    var oldLabel = segStr.getData();
    var edge = new jsts.geomgraph.Edge(segStr.getCoordinates(),
        new jsts.geomgraph.Label(oldLabel));
    this.insertUniqueEdge(edge);
  }
};


/**
 * Inserted edges are checked to see if an identical edge already exists. If so,
 * the edge is not inserted, but its label is merged with the existing edge.
 *
 * @protected
 */
jsts.operation.buffer.BufferBuilder.prototype.insertUniqueEdge = function(e) {
  var existingEdge = this.edgeList.findEqualEdge(e);

  // If an identical edge already exists, simply update its label
  if (existingEdge != null) {
    var existingLabel = existingEdge.getLabel();

    var labelToMerge = e.getLabel();
    // check if new edge is in reverse direction to existing edge
    // if so, must flip the label before merging it
    if (!existingEdge.isPointwiseEqual(e)) {
      labelToMerge = new jsts.geomgraph.Label(e.getLabel());
      labelToMerge.flip();
    }
    existingLabel.merge(labelToMerge);

    // compute new depth delta of sum of edges
    var mergeDelta = jsts.operation.buffer.BufferBuilder
        .depthDelta(labelToMerge);
    var existingDelta = existingEdge.getDepthDelta();
    var newDelta = existingDelta + mergeDelta;
    existingEdge.setDepthDelta(newDelta);
  } else {
    // no matching existing edge was found
    // add this new edge to the list of edges in this graph
    this.edgeList.add(e);
    e.setDepthDelta(jsts.operation.buffer.BufferBuilder
        .depthDelta(e.getLabel()));
  }
};


/**
 * @param {PlanarGraph}
 *          graph
 * @private
 */
jsts.operation.buffer.BufferBuilder.prototype.createSubgraphs = function(graph) {
  var subgraphList = [];
  for (var i = graph.getNodes().iterator(); i.hasNext();) {
    var node = i.next();
    if (!node.isVisited()) {
      var subgraph = new jsts.operation.buffer.BufferSubgraph();
      subgraph.create(node);
      subgraphList.push(subgraph);
    }
  }
  /**
   * Sort the subgraphs in descending order of their rightmost coordinate. This
   * ensures that when the Polygons for the subgraphs are built, subgraphs for
   * shells will have been built before the subgraphs for any holes they
   * contain.
   */

  var compare = function(a, b) {
    return a.compareTo(b);
  };
  subgraphList.sort(compare);
  subgraphList.reverse();
  return subgraphList;
};


/**
 * Completes the building of the input subgraphs by depth-labelling them, and
 * adds them to the PolygonBuilder. The subgraph list must be sorted in
 * rightmost-coordinate order.
 *
 * @param {Array}
 *          subgraphList the subgraphs to build.
 * @param {PolygonBuilder}
 *          polyBuilder the PolygonBuilder which will build the final polygons.
 * @private
 */
jsts.operation.buffer.BufferBuilder.prototype.buildSubgraphs = function(
    subgraphList, polyBuilder) {
  var processedGraphs = [];
  for (var i = 0; i < subgraphList.length; i++) {
    var subgraph = subgraphList[i];
    var p = subgraph.getRightmostCoordinate();
    var locater = new jsts.operation.buffer.SubgraphDepthLocater(
        processedGraphs);
    var outsideDepth = locater.getDepth(p);
    subgraph.computeDepth(outsideDepth);
    subgraph.findResultEdges();
    processedGraphs.push(subgraph);
    polyBuilder.add(subgraph.getDirectedEdges(), subgraph.getNodes());
  }
};


/**
 * @private
 */
jsts.operation.buffer.BufferBuilder.convertSegStrings = function(it) {
  var fact = new jsts.geom.GeometryFactory();
  var lines = new javascript.util.ArrayList();
  while (it.hasNext()) {
    var ss = it.next();
    var line = fact.createLineString(ss.getCoordinates());
    lines.add(line);
  }
  return fact.buildGeometry(lines);
};


/**
 * Gets the standard result for an empty buffer. Since buffer always returns a
 * polygonal result, this is chosen to be an empty polygon.
 *
 * @return the empty result geometry.
 * @private
 */
jsts.operation.buffer.BufferBuilder.prototype.createEmptyResultGeometry = function() {
  var emptyGeom = this.geomFact.createPolygon(null, null);
  return emptyGeom;
};



//    JSTS EdgeList


  /**
   * @requires jsts/util/Assert.js
   */

  var ArrayList = javascript.util.ArrayList;
  var TreeMap = javascript.util.TreeMap;

  /**
   * A EdgeList is a list of Edges. It supports locating edges that are
   * pointwise equals to a target edge.
   *
   * @constructor
   */
  jsts.geomgraph.EdgeList = function() {
    this.edges = new javascript.util.ArrayList();
    this.ocaMap = new javascript.util.TreeMap();
  };


  /**
   * @type {javascript.util.ArrayList}
   * @private
   */
  jsts.geomgraph.EdgeList.prototype.edges = null;


  /**
   * An index of the edges, for fast lookup.
   *
   * @type {javascript.util.HashMap}
   * @private
   */
  jsts.geomgraph.EdgeList.prototype.ocaMap = null;


  /**
   * Insert an edge unless it is already in the list
   */
  jsts.geomgraph.EdgeList.prototype.add = function(e) {
    this.edges.add(e);
    var oca = new jsts.noding.OrientedCoordinateArray(e.getCoordinates());
    this.ocaMap.put(oca, e);
  };

  jsts.geomgraph.EdgeList.prototype.addAll = function(edgeColl) {
    for (var i = edgeColl.iterator(); i.hasNext();) {
      this.add(i.next());
    }
  };


  /**
   * @return {javascript.util.List}
   */
  jsts.geomgraph.EdgeList.prototype.getEdges = function() {
    return this.edges;
  };


  /**
   * If there is an edge equal to e already in the list, return it. Otherwise
   * return null.
   *
   * @param {Edge}
   *          e
   * @return {Edge} equal edge, if there is one already in the list null
   *         otherwise.
   */
  jsts.geomgraph.EdgeList.prototype.findEqualEdge = function(e) {
    var oca = new jsts.noding.OrientedCoordinateArray(e.getCoordinates());
    // will return null if no edge matches
    var matchEdge = this.ocaMap.get(oca);
    return matchEdge;
  };

  jsts.geomgraph.EdgeList.prototype.getEdges = function() {
    return this.edges;
  };

  jsts.geomgraph.EdgeList.prototype.iterator = function() {
    return this.edges.iterator();
  };

  jsts.geomgraph.EdgeList.prototype.get = function(i) {
    return this.edges.get(i);
  };


  /**
   * If the edge e is already in the list, return its index.
   *
   * @return {Number} index, if e is already in the list -1 otherwise.
   */
  jsts.geomgraph.EdgeList.prototype.findEdgeIndex = function(e) {
    for (var i = 0; i < this.edges.size(); i++) {
      if (this.edges.get(i).equals(e))
        return i;
    }
    return -1;
  };




//    JSTS OffsetCurveBuilder

/* Copyright (c) 2011 by The Authors.
 * Published under the LGPL 2.1 license.
 * See /license-notice.txt for the full text of the license notice.
 * See /license.txt for the full text of the license.
 */



/**
 * Computes the raw offset curve for a single {@link Geometry} component (ring,
 * line or point). A raw offset curve line is not noded - it may contain
 * self-intersections (and usually will). The final buffer polygon is computed
 * by forming a topological graph of all the noded raw curves and tracing
 * outside contours. The points in the raw curve are rounded to a given
 * {@link PrecisionModel}.
 *
 * @constructor
 */
jsts.operation.buffer.OffsetCurveBuilder = function(precisionModel, bufParams) {
  this.precisionModel = precisionModel;
  this.bufParams = bufParams;
};


/**
 * @type {double}
 * @private
 */
jsts.operation.buffer.OffsetCurveBuilder.prototype.distance = 0.0;


/**
 * @type {PrecisionModel}
 * @private
 */
jsts.operation.buffer.OffsetCurveBuilder.prototype.precisionModel = null;


/**
 * @type {BufferParameters}
 * @private
 */
jsts.operation.buffer.OffsetCurveBuilder.prototype.bufParams = null;


/**
 * Gets the buffer parameters being used to generate the curve.
 *
 * @return the buffer parameters being used.
 */
jsts.operation.buffer.OffsetCurveBuilder.prototype.getBufferParameters = function() {
  return this.bufParams;
};


/**
 * This method handles single points as well as LineStrings. LineStrings are
 * assumed <b>not</b> to be closed (the function will not fail for closed
 * lines, but will generate superfluous line caps).
 *
 * @param inputPts
 *          the vertices of the line to offset.
 * @param distance
 *          the offset distance.
 *
 * @return a Coordinate array representing the curve.
 * @return null if the curve is empty.
 */
jsts.operation.buffer.OffsetCurveBuilder.prototype.getLineCurve = function(
    inputPts, distance) {
  this.distance = distance;

  // a zero or negative width buffer of a line/point is empty
  if (this.distance < 0.0 && !this.bufParams.isSingleSided())
    return null;
  if (this.distance == 0.0)
    return null;

  var posDistance = Math.abs(this.distance);
  var segGen = this.getSegGen(posDistance);
  if (inputPts.length <= 1) {
    this.computePointCurve(inputPts[0], segGen);
  } else {
    if (this.bufParams.isSingleSided()) {
      var isRightSide = distance < 0.0;
      this.computeSingleSidedBufferCurve(inputPts, isRightSide, segGen);
    } else
      this.computeLineBufferCurve(inputPts, segGen);
  }

  var lineCoord = segGen.getCoordinates();
  return lineCoord;
};


/**
 * This method handles the degenerate cases of single points and lines, as well
 * as rings.
 *
 * @return a Coordinate array representing the curve.
 * @return null if the curve is empty.
 */
jsts.operation.buffer.OffsetCurveBuilder.prototype.getRingCurve = function(
    inputPts, side, distance) {
  this.distance = distance;
  if (inputPts.length <= 2)
    return this.getLineCurve(inputPts, distance);

  // optimize creating ring for for zero distance
  if (this.distance == 0.0) {
    return jsts.operation.buffer.OffsetCurveBuilder.copyCoordinates(inputPts);
  }
  var segGen = this.getSegGen(this.distance);
  this.computeRingBufferCurve(inputPts, side, segGen);
  return segGen.getCoordinates();
};

jsts.operation.buffer.OffsetCurveBuilder.prototype.getOffsetCurve = function(
    inputPts, distance) {
  this.distance = distance;

  // a zero width offset curve is empty
  if (this.distance === 0.0)
    return null;

  var isRightSide = this.distance < 0.0;
  var posDistance = Math.abs(this.distance);
  var segGen = this.getSegGen(posDistance);
  if (inputPts.length <= 1) {
    this.computePointCurve(inputPts[0], segGen);
  } else {
    this.computeOffsetCurve(inputPts, isRightSide, segGen);
  }
  var curvePts = segGen.getCoordinates();
  // for right side line is traversed in reverse direction, so have to reverse
  // generated line
  if (isRightSide)
    curvePts.reverse();
  return curvePts;
};


/**
 * @private
 */
jsts.operation.buffer.OffsetCurveBuilder.copyCoordinates = function(pts) {
  var copy = [];
  for (var i = 0; i < pts.length; i++) {
    copy.push(pts[i].clone());
  }
  return copy;
};


/**
 * @private
 */
jsts.operation.buffer.OffsetCurveBuilder.prototype.getSegGen = function(
    distance) {
  return new jsts.operation.buffer.OffsetSegmentGenerator(this.precisionModel,
      this.bufParams, distance);
};


/**
 * Use a value which results in a potential distance error which is
 * significantly less than the error due to the quadrant segment discretization.
 * For QS = 8 a value of 100 is reasonable. This should produce a maximum of 1%
 * distance error.
 *
 * @private
 */
jsts.operation.buffer.OffsetCurveBuilder.SIMPLIFY_FACTOR = 100.0;


/**
 * Computes the distance tolerance to use during input line simplification.
 *
 * @param distance
 *          the buffer distance.
 * @return the simplification tolerance.
 * @private
 */
jsts.operation.buffer.OffsetCurveBuilder.simplifyTolerance = function(
    bufDistance) {
  return bufDistance / jsts.operation.buffer.OffsetCurveBuilder.SIMPLIFY_FACTOR;
};


/**
 * @private
 */
jsts.operation.buffer.OffsetCurveBuilder.prototype.computePointCurve = function(
    pt, segGen) {
  switch (this.bufParams.getEndCapStyle()) {
  case jsts.operation.buffer.BufferParameters.CAP_ROUND:
    segGen.createCircle(pt);
    break;
  case jsts.operation.buffer.BufferParameters.CAP_SQUARE:
    segGen.createSquare(pt);
    break;
  // otherwise curve is empty (e.g. for a butt cap);
  }
};


/**
 * @private
 */
jsts.operation.buffer.OffsetCurveBuilder.prototype.computeLineBufferCurve = function(
    inputPts, segGen) {
  var distTol = jsts.operation.buffer.OffsetCurveBuilder
      .simplifyTolerance(this.distance);

  // --------- compute points for left side of line
  // Simplify the appropriate side of the line before generating
  var simp1 = jsts.operation.buffer.BufferInputLineSimplifier.simplify(
      inputPts, distTol);
  // MD - used for testing only (to eliminate simplification)
  // Coordinate[] simp1 = inputPts;

  var n1 = simp1.length - 1;
  segGen.initSideSegments(simp1[0], simp1[1], jsts.geomgraph.Position.LEFT);
  for (var i = 2; i <= n1; i++) {
    segGen.addNextSegment(simp1[i], true);
  }
  segGen.addLastSegment();
  // add line cap for end of line
  segGen.addLineEndCap(simp1[n1 - 1], simp1[n1]);

  // ---------- compute points for right side of line
  // Simplify the appropriate side of the line before generating
  var simp2 = jsts.operation.buffer.BufferInputLineSimplifier.simplify(
      inputPts, -distTol);
  // MD - used for testing only (to eliminate simplification)
  // Coordinate[] simp2 = inputPts;
  var n2 = simp2.length - 1;

  // since we are traversing line in opposite order, offset position is still
  // LEFT
  segGen.initSideSegments(simp2[n2], simp2[n2 - 1], jsts.geomgraph.Position.LEFT);
  for (var i = n2 - 2; i >= 0; i--) {
    segGen.addNextSegment(simp2[i], true);
  }
  segGen.addLastSegment();
  // add line cap for start of line
  segGen.addLineEndCap(simp2[1], simp2[0]);

  segGen.closeRing();
};


/**
 * @private
 */
jsts.operation.buffer.OffsetCurveBuilder.prototype.computeSingleSidedBufferCurve = function(
    inputPts, isRightSide, segGen) {
  var distTol = jsts.operation.buffer.OffsetCurveBuilder
      .simplifyTolerance(this.distance);

  if (isRightSide) {
    // add original line
    segGen.addSegments(inputPts, true);

    // ---------- compute points for right side of line
    // Simplify the appropriate side of the line before generating
    var simp2 = jsts.operation.buffer.BufferInputLineSimplifier.simplify(
        inputPts, -distTol);
    // MD - used for testing only (to eliminate simplification)
    // Coordinate[] simp2 = inputPts;
    var n2 = simp2.length - 1;

    // since we are traversing line in opposite order, offset position is still
    // LEFT
    segGen.initSideSegments(simp2[n2], simp2[n2 - 1],
        jsts.geomgraph.Position.LEFT);
    segGen.addFirstSegment();
    for (var i = n2 - 2; i >= 0; i--) {
      segGen.addNextSegment(simp2[i], true);
    }
  } else {
    // add original line
    segGen.addSegments(inputPts, false);

    // --------- compute points for left side of line
    // Simplify the appropriate side of the line before generating
    var simp1 = jsts.operation.buffer.BufferInputLineSimplifier.simplify(
        inputPts, distTol);
    // MD - used for testing only (to eliminate simplification)
    // Coordinate[] simp1 = inputPts;

    var n1 = simp1.length - 1;
    segGen.initSideSegments(simp1[0], simp1[1], jsts.geomgraph.Position.LEFT);
    segGen.addFirstSegment();
    for (var i = 2; i <= n1; i++) {
      segGen.addNextSegment(simp1[i], true);
    }
  }
  segGen.addLastSegment();
  segGen.closeRing();
};


/**
 * @private
 */
jsts.operation.buffer.OffsetCurveBuilder.prototype.computeOffsetCurve = function(
    inputPts, isRightSide, segGen) {
  var distTol = jsts.operation.buffer.OffsetCurveBuilder
      .simplifyTolerance(this.distance);

  if (isRightSide) {
    // ---------- compute points for right side of line
    // Simplify the appropriate side of the line before generating
    var simp2 = jsts.operation.buffer.BufferInputLineSimplifier.simplify(
        inputPts, -distTol);
    // MD - used for testing only (to eliminate simplification)
    // Coordinate[] simp2 = inputPts;
    var n2 = simp2.length - 1;

    // since we are traversing line in opposite order, offset position is still
    // LEFT
    segGen.initSideSegments(simp2[n2], simp2[n2 - 1],
        jsts.geomgraph.Position.LEFT);
    segGen.addFirstSegment();
    for (var i = n2 - 2; i >= 0; i--) {
      segGen.addNextSegment(simp2[i], true);
    }
  } else {
    // --------- compute points for left side of line
    // Simplify the appropriate side of the line before generating
    var simp1 = jsts.operation.buffer.BufferInputLineSimplifier.simplify(
        inputPts, distTol);
    // MD - used for testing only (to eliminate simplification)
    // Coordinate[] simp1 = inputPts;

    var n1 = simp1.length - 1;
    segGen.initSideSegments(simp1[0], simp1[1], jsts.geomgraph.Position.LEFT);
    segGen.addFirstSegment();
    for (var i = 2; i <= n1; i++) {
      segGen.addNextSegment(simp1[i], true);
    }
  }
  segGen.addLastSegment();
};


/**
 * @private
 */
jsts.operation.buffer.OffsetCurveBuilder.prototype.computeRingBufferCurve = function(
    inputPts, side, segGen) {
  // simplify input line to improve performance
  var distTol = jsts.operation.buffer.OffsetCurveBuilder
      .simplifyTolerance(this.distance);
  // ensure that correct side is simplified
  if (side === jsts.geomgraph.Position.RIGHT)
    distTol = -distTol;
  var simp = jsts.operation.buffer.BufferInputLineSimplifier.simplify(inputPts,
      distTol);

  var n = simp.length - 1;
  segGen.initSideSegments(simp[n - 1], simp[0], side);
  for (var i = 1; i <= n; i++) {
    var addStartPoint = i !== 1;
    segGen.addNextSegment(simp[i], addStartPoint);
  }
  segGen.closeRing();
};




//    JSTS OffsetCurveSetBuilder

/* Copyright (c) 2011 by The Authors.
 * Published under the LGPL 2.1 license.
 * See /license-notice.txt for the full text of the license notice.
 * See /license.txt for the full text of the license.
 */


/**
 * Creates all the raw offset curves for a buffer of a {@link Geometry}. Raw
 * curves need to be noded together and polygonized to form the final buffer
 * area.
 *
 * @constructor
 */
jsts.operation.buffer.OffsetCurveSetBuilder = function(inputGeom, distance,
    curveBuilder) {
  this.inputGeom = inputGeom;
  this.distance = distance;
  this.curveBuilder = curveBuilder;

  this.curveList = new javascript.util.ArrayList();
};


/**
 * @private
 */
jsts.operation.buffer.OffsetCurveSetBuilder.prototype.inputGeom = null;


/**
 * @private
 */
jsts.operation.buffer.OffsetCurveSetBuilder.prototype.distance = null;


/**
 * @private
 */
jsts.operation.buffer.OffsetCurveSetBuilder.prototype.curveBuilder = null;


/**
 * @private
 */
jsts.operation.buffer.OffsetCurveSetBuilder.prototype.curveList = null;


/**
 * Computes the set of raw offset curves for the buffer. Each offset curve has
 * an attached {@link Label} indicating its left and right location.
 *
 * @return a Collection of SegmentStrings representing the raw buffer curves.
 */
jsts.operation.buffer.OffsetCurveSetBuilder.prototype.getCurves = function() {
  this.add(this.inputGeom);
  return this.curveList;
};


/**
 * Creates a {@link SegmentString} for a coordinate list which is a raw offset
 * curve, and adds it to the list of buffer curves. The SegmentString is tagged
 * with a Label giving the topology of the curve. The curve may be oriented in
 * either direction. If the curve is oriented CW, the locations will be: <br>
 * Left: Location.EXTERIOR <br>
 * Right: Location.INTERIOR
 *
 * @private
 */
jsts.operation.buffer.OffsetCurveSetBuilder.prototype.addCurve = function(
    coord, leftLoc, rightLoc) {
  // don't add null or trivial curves
  if (coord == null || coord.length < 2)
    return;
  // add the edge for a coordinate list which is a raw offset curve
  var e = new jsts.noding.NodedSegmentString(coord, new jsts.geomgraph.Label(0,
      jsts.geom.Location.BOUNDARY, leftLoc, rightLoc));
  this.curveList.add(e);
};


/**
 * @private
 */
jsts.operation.buffer.OffsetCurveSetBuilder.prototype.add = function(g) {
  if (g.isEmpty())
    return;

  if (g instanceof jsts.geom.Polygon)
    this.addPolygon(g);
  // LineString also handles LinearRings
  else if (g instanceof jsts.geom.LineString)
    this.addLineString(g);
  else if (g instanceof jsts.geom.Point)
    this.addPoint(g);
  else if (g instanceof jsts.geom.MultiPoint)
    this.addCollection(g);
  else if (g instanceof jsts.geom.MultiLineString)
    this.addCollection(g);
  else if (g instanceof jsts.geom.MultiPolygon)
    this.addCollection(g);
  else if (g instanceof jsts.geom.GeometryCollection)
    this.addCollection(g);
  else
    throw new jsts.error.IllegalArgumentError();
};


/**
 * @private
 */
jsts.operation.buffer.OffsetCurveSetBuilder.prototype.addCollection = function(
    gc) {
  for (var i = 0; i < gc.getNumGeometries(); i++) {
    var g = gc.getGeometryN(i);
    this.add(g);
  }
};


/**
 * Add a Point to the graph.
 *
 * @private
 */
jsts.operation.buffer.OffsetCurveSetBuilder.prototype.addPoint = function(p) {
  // a zero or negative width buffer of a line/point is empty
  if (this.distance <= 0.0)
    return;
  var coord = p.getCoordinates();
  var curve = this.curveBuilder.getLineCurve(coord, this.distance);
  this.addCurve(curve, jsts.geom.Location.EXTERIOR, jsts.geom.Location.INTERIOR);
};


/**
 * @private
 */
jsts.operation.buffer.OffsetCurveSetBuilder.prototype.addLineString = function(
    line) {
  // a zero or negative width buffer of a line/point is empty
  if (this.distance <= 0.0 &&
      !this.curveBuilder.getBufferParameters().isSingleSided())
    return;
  var coord = jsts.geom.CoordinateArrays.removeRepeatedPoints(line
      .getCoordinates());
  var curve = this.curveBuilder.getLineCurve(coord, this.distance);
  this
      .addCurve(curve, jsts.geom.Location.EXTERIOR, jsts.geom.Location.INTERIOR);
};


/**
 * @private
 */
jsts.operation.buffer.OffsetCurveSetBuilder.prototype.addPolygon = function(p) {
  var offsetDistance = this.distance;
  var offsetSide = jsts.geomgraph.Position.LEFT;
  if (this.distance < 0.0) {
    offsetDistance = -this.distance;
    offsetSide = jsts.geomgraph.Position.RIGHT;
  }

  var shell = p.getExteriorRing();
  var shellCoord = jsts.geom.CoordinateArrays.removeRepeatedPoints(shell
      .getCoordinates());
  // optimization - don't bother computing buffer
  // if the polygon would be completely eroded
  if (this.distance < 0.0 && this.isErodedCompletely(shell, this.distance))
    return;
  // don't attemtp to buffer a polygon with too few distinct vertices
  if (this.distance <= 0.0 && shellCoord.length < 3)
    return;

  this.addPolygonRing(shellCoord, offsetDistance, offsetSide,
      jsts.geom.Location.EXTERIOR, jsts.geom.Location.INTERIOR);

  for (var i = 0; i < p.getNumInteriorRing(); i++) {

    var hole = p.getInteriorRingN(i);
    var holeCoord = jsts.geom.CoordinateArrays.removeRepeatedPoints(hole
        .getCoordinates());

    // optimization - don't bother computing buffer for this hole
    // if the hole would be completely covered
    if (this.distance > 0.0 && this.isErodedCompletely(hole, -this.distance))
      continue;

    // Holes are topologically labelled opposite to the shell, since
    // the interior of the polygon lies on their opposite side
    // (on the left, if the hole is oriented CCW)
    this.addPolygonRing(holeCoord, offsetDistance, jsts.geomgraph.Position
        .opposite(offsetSide), jsts.geom.Location.INTERIOR,
        jsts.geom.Location.EXTERIOR);
  }
};


/**
 * Adds an offset curve for a polygon ring. The side and left and right
 * topological location arguments assume that the ring is oriented CW. If the
 * ring is in the opposite orientation, the left and right locations must be
 * interchanged and the side flipped.
 *
 * @param coord
 *          the coordinates of the ring (must not contain repeated points).
 * @param offsetDistance
 *          the distance at which to create the buffer.
 * @param side
 *          the side of the ring on which to construct the buffer line.
 * @param cwLeftLoc
 *          the location on the L side of the ring (if it is CW).
 * @param cwRightLoc
 *          the location on the R side of the ring (if it is CW).
 */
/**
 * @private
 */
jsts.operation.buffer.OffsetCurveSetBuilder.prototype.addPolygonRing = function(
    coord, offsetDistance, side, cwLeftLoc, cwRightLoc) {
  // don't bother adding ring if it is "flat" and will disappear in the output
  if (offsetDistance == 0.0 &&
      coord.length < jsts.geom.LinearRing.MINIMUM_VALID_SIZE)
    return;

  var leftLoc = cwLeftLoc;
  var rightLoc = cwRightLoc;
  if (coord.length >= jsts.geom.LinearRing.MINIMUM_VALID_SIZE &&
      jsts.algorithm.CGAlgorithms.isCCW(coord)) {
    leftLoc = cwRightLoc;
    rightLoc = cwLeftLoc;
    side = jsts.geomgraph.Position.opposite(side);
  }
  var curve = this.curveBuilder.getRingCurve(coord, side, offsetDistance);
  this.addCurve(curve, leftLoc, rightLoc);
};


/**
 * The ringCoord is assumed to contain no repeated points. It may be degenerate
 * (i.e. contain only 1, 2, or 3 points). In this case it has no area, and hence
 * has a minimum diameter of 0.
 *
 * @param ringCoord
 * @param offsetDistance
 * @return
 */
/**
 * @private
 */
jsts.operation.buffer.OffsetCurveSetBuilder.prototype.isErodedCompletely = function(
    ring, bufferDistance) {
  var ringCoord = ring.getCoordinates();
  var minDiam = 0.0;
  // degenerate ring has no area
  if (ringCoord.length < 4)
    return bufferDistance < 0;

  // important test to eliminate inverted triangle bug
  // also optimizes erosion test for triangles
  if (ringCoord.length == 4)
    return this.isTriangleErodedCompletely(ringCoord, bufferDistance);

  // if envelope is narrower than twice the buffer distance, ring is eroded
  var env = ring.getEnvelopeInternal();
  var envMinDimension = Math.min(env.getHeight(), env.getWidth());
  if (bufferDistance < 0.0 && 2 * Math.abs(bufferDistance) > envMinDimension)
    return true;

  return false;
};


/**
 * Tests whether a triangular ring would be eroded completely by the given
 * buffer distance. This is a precise test. It uses the fact that the inner
 * buffer of a triangle converges on the inCentre of the triangle (the point
 * equidistant from all sides). If the buffer distance is greater than the
 * distance of the inCentre from a side, the triangle will be eroded completely.
 *
 * This test is important, since it removes a problematic case where the buffer
 * distance is slightly larger than the inCentre distance. In this case the
 * triangle buffer curve "inverts" with incorrect topology, producing an
 * incorrect hole in the buffer.
 *
 * @param triangleCoord
 * @param bufferDistance
 * @return
 *
 * @private
 */
jsts.operation.buffer.OffsetCurveSetBuilder.prototype.isTriangleErodedCompletely = function(
    triangleCoord, bufferDistance) {
  var tri = new jsts.geom.Triangle(triangleCoord[0], triangleCoord[1], triangleCoord[2]);
  var inCentre = tri.inCentre();
  var distToCentre = jsts.algorithm.CGAlgorithms.distancePointLine(inCentre,
      tri.p0, tri.p1);
  return distToCentre < Math.abs(bufferDistance);
};




//    OffsetSegmentGenerator

/* Copyright (c) 2011 by The Authors.
 * Published under the LGPL 2.1 license.
 * See /license-notice.txt for the full text of the license notice.
 * See /license.txt for the full text of the license.
 */


/**
 * Generates segments which form an offset curve. Supports all end cap and join
 * options provided for buffering. Implements various heuristics to produce
 * smoother, simpler curves which are still within a reasonable tolerance of the
 * true curve.
 * @constructor
 */
jsts.operation.buffer.OffsetSegmentGenerator = function(precisionModel,
    bufParams, distance) {
  this.seg0 = new jsts.geom.LineSegment();
  this.seg1 = new jsts.geom.LineSegment();
  this.offset0 = new jsts.geom.LineSegment();
  this.offset1 = new jsts.geom.LineSegment();

  this.precisionModel = precisionModel;
  this.bufParams = bufParams;

  // compute intersections in full precision, to provide accuracy
  // the points are rounded as they are inserted into the curve line
  this.li = new jsts.algorithm.RobustLineIntersector();
  this.filletAngleQuantum = Math.PI / 2.0 / bufParams.getQuadrantSegments();

  /**
   * Non-round joins cause issues with short closing segments, so don't use
   * them. In any case, non-round joins only really make sense for relatively
   * small buffer distances.
   */
  if (this.bufParams.getQuadrantSegments() >= 8 &&
      this.bufParams.getJoinStyle() === jsts.operation.buffer.BufferParameters.JOIN_ROUND) {
    this.closingSegLengthFactor = jsts.operation.buffer.OffsetSegmentGenerator.MAX_CLOSING_SEG_LEN_FACTOR;
  }
  this.init(distance);
};


/**
 * Factor which controls how close offset segments can be to skip adding a
 * filler or mitre.
 *
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.OFFSET_SEGMENT_SEPARATION_FACTOR = 1.0E-3;


/**
 * Factor which controls how close curve vertices on inside turns can be to be
 * snapped
 *
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.INSIDE_TURN_VERTEX_SNAP_DISTANCE_FACTOR = 1.0E-3;


/**
 * Factor which controls how close curve vertices can be to be snapped
 *
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.CURVE_VERTEX_SNAP_DISTANCE_FACTOR = 1.0E-6;


/**
 * Factor which determines how short closing segs can be for round buffers *
 *
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.MAX_CLOSING_SEG_LEN_FACTOR = 80;


/**
 * the max error of approximation (distance) between a quad segment and the true
 * fillet curve
 *
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.maxCurveSegmentError = 0.0;


/**
 * The angle quantum with which to approximate a fillet curve (based on the
 * input # of quadrant segments)
 *
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.filletAngleQuantum = null;


/**
 * The Closing Segment Length Factor controls how long "closing segments" are.
 * Closing segments are added at the middle of inside corners to ensure a
 * smoother boundary for the buffer offset curve. In some cases (particularly
 * for round joins with default-or-better quantization) the closing segments can
 * be made quite short. This substantially improves performance (due to fewer
 * intersections being created).
 *
 * A closingSegFactor of 0 results in lines to the corner vertex A
 * closingSegFactor of 1 results in lines halfway to the corner vertex A
 * closingSegFactor of 80 results in lines 1/81 of the way to the corner vertex
 * (this option is reasonable for the very common default situation of round
 * joins and quadrantSegs >= 8)
 *
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.closingSegLengthFactor = 1;


/**
 * @type {OffsetSegmentString}
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.segList = null;


/**
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.distance = 0.0;


/**
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.precisionModel = null;


/**
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.bufParams = null;


/**
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.li = null;


/**
 * @type {Coordinate}
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.s0 = null;


/**
 * @type {Coordinate}
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.s1 = null;


/**
 * @type {Coordinate}
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.s2 = null;


/**
 * @type {LineSegment}
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.seg0 = null;


/**
 * @type {LineSegment}
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.seg1 = null;


/**
 * @type {LineSegment}
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.offset0 = null;


/**
 * @type {LineSegment}
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.offset1 = null;


/**
 * @type {number}
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.side = 0;


/**
 * @type {boolean}
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.hasNarrowConcaveAngle = false;


/**
 * Tests whether the input has a narrow concave angle (relative to the offset
 * distance). In this case the generated offset curve will contain
 * self-intersections and heuristic closing segments. This is expected behaviour
 * in the case of buffer curves. For pure offset curves, the output needs to be
 * further treated before it can be used.
 *
 * @return true if the input has a narrow concave angle.
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.hasNarrowConcaveAngle = function() {
  return this.hasNarrowConcaveAngle;
};


/**
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.init = function(distance) {
  this.distance = distance;
  this.maxCurveSegmentError = this.distance *
      (1 - Math.cos(this.filletAngleQuantum / 2.0));
  this.segList = new jsts.operation.buffer.OffsetSegmentString();
  this.segList.setPrecisionModel(this.precisionModel);
  /**
   * Choose the min vertex separation as a small fraction of the offset
   * distance.
   */
  this.segList
      .setMinimumVertexDistance(this.distance *
          jsts.operation.buffer.OffsetSegmentGenerator.CURVE_VERTEX_SNAP_DISTANCE_FACTOR);
};


jsts.operation.buffer.OffsetSegmentGenerator.prototype.initSideSegments = function(
    s1, s2, side) {
  this.s1 = s1;
  this.s2 = s2;
  this.side = side;
  this.seg1.setCoordinates(this.s1, this.s2);
  this.computeOffsetSegment(this.seg1, this.side, this.distance, this.offset1);
};

jsts.operation.buffer.OffsetSegmentGenerator.prototype.getCoordinates = function() {
  return this.segList.getCoordinates();
};

jsts.operation.buffer.OffsetSegmentGenerator.prototype.closeRing = function() {
  this.segList.closeRing();
};

jsts.operation.buffer.OffsetSegmentGenerator.prototype.addSegments = function(
    pt, isForward) {
  this.segList.addPts(pt, isForward);
};

jsts.operation.buffer.OffsetSegmentGenerator.prototype.addFirstSegment = function() {
  this.segList.addPt(this.offset1.p0);
};


/**
 * Add last offset point
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.addLastSegment = function() {
  this.segList.addPt(this.offset1.p1);
};

jsts.operation.buffer.OffsetSegmentGenerator.prototype.addNextSegment = function(
    p, addStartPoint) {
  // s0-s1-s2 are the coordinates of the previous segment and the current one
  this.s0 = this.s1;
  this.s1 = this.s2;
  this.s2 = p;
  this.seg0.setCoordinates(this.s0, this.s1);
  this.computeOffsetSegment(this.seg0, this.side, this.distance, this.offset0);
  this.seg1.setCoordinates(this.s1, this.s2);
  this.computeOffsetSegment(this.seg1, this.side, this.distance, this.offset1);

  // do nothing if points are equal
  if (this.s1.equals(this.s2))
    return;

  var orientation = jsts.algorithm.CGAlgorithms.computeOrientation(this.s0,
      this.s1, this.s2);
  var outsideTurn = (orientation === jsts.algorithm.CGAlgorithms.CLOCKWISE && this.side === jsts.geomgraph.Position.LEFT) ||
      (orientation === jsts.algorithm.CGAlgorithms.COUNTERCLOCKWISE && this.side === jsts.geomgraph.Position.RIGHT);

  if (orientation == 0) { // lines are collinear
    this.addCollinear(addStartPoint);
  } else if (outsideTurn) {
    this.addOutsideTurn(orientation, addStartPoint);
  } else { // inside turn
    this.addInsideTurn(orientation, addStartPoint);
  }
};


/**
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.addCollinear = function(
    addStartPoint) {
  /**
   * This test could probably be done more efficiently, but the situation of
   * exact collinearity should be fairly rare.
   */
  this.li.computeIntersection(this.s0, this.s1, this.s1, this.s2);
  var numInt = this.li.getIntersectionNum();
  /**
   * if numInt is < 2, the lines are parallel and in the same direction. In this
   * case the point can be ignored, since the offset lines will also be
   * parallel.
   */
  if (numInt >= 2) {
    /**
     * segments are collinear but reversing. Add an "end-cap" fillet all the way
     * around to other direction This case should ONLY happen for LineStrings,
     * so the orientation is always CW. (Polygons can never have two consecutive
     * segments which are parallel but reversed, because that would be a self
     * intersection.
     *
     */
    if (this.bufParams.getJoinStyle() === jsts.operation.buffer.BufferParameters.JOIN_BEVEL ||
        this.bufParams.getJoinStyle() === jsts.operation.buffer.BufferParameters.JOIN_MITRE) {
      if (addStartPoint)
        this.segList.addPt(this.offset0.p1);
      this.segList.addPt(this.offset1.p0);
    } else {
      this.addFillet(this.s1, this.offset0.p1, this.offset1.p0,
          jsts.algorithm.CGAlgorithms.CLOCKWISE, this.distance);
    }
  }
};


/**
 * Adds the offset points for an outside (convex) turn
 *
 * @param orientation
 * @param addStartPoint
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.addOutsideTurn = function(
    orientation, addStartPoint) {
  /**
   * Heuristic: If offset endpoints are very close together, just use one of
   * them as the corner vertex. This avoids problems with computing mitre
   * corners in the case where the two segments are almost parallel (which is
   * hard to compute a robust intersection for).
   */
  if (this.offset0.p1.distance(this.offset1.p0) < this.distance *
      jsts.operation.buffer.OffsetSegmentGenerator.OFFSET_SEGMENT_SEPARATION_FACTOR) {
    this.segList.addPt(this.offset0.p1);
    return;
  }

  if (this.bufParams.getJoinStyle() === jsts.operation.buffer.BufferParameters.JOIN_MITRE) {
    this.addMitreJoin(this.s1, this.offset0, this.offset1, this.distance);
  } else if (this.bufParams.getJoinStyle() === jsts.operation.buffer.BufferParameters.JOIN_BEVEL) {
    this.addBevelJoin(this.offset0, this.offset1);
  } else {
    // add a circular fillet connecting the endpoints of the offset segments
    if (addStartPoint)
      this.segList.addPt(this.offset0.p1);
    // TESTING - comment out to produce beveled joins
    this.addFillet(this.s1, this.offset0.p1, this.offset1.p0, orientation,
        this.distance);
    this.segList.addPt(this.offset1.p0);
  }
};


/**
 * Adds the offset points for an inside (concave) turn.
 *
 * @param orientation
 * @param addStartPoint
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.addInsideTurn = function(
    orientation, addStartPoint) {
  /**
   * add intersection point of offset segments (if any)
   */
  this.li.computeIntersection(this.offset0.p0, this.offset0.p1, this.offset1.p0, this.offset1.p1);
  if (this.li.hasIntersection()) {
    this.segList.addPt(this.li.getIntersection(0));
  } else {
    /**
     * If no intersection is detected, it means the angle is so small and/or the
     * offset so large that the offsets segments don't intersect. In this case
     * we must add a "closing segment" to make sure the buffer curve is
     * continuous, fairly smooth (e.g. no sharp reversals in direction) and
     * tracks the buffer correctly around the corner. The curve connects the
     * endpoints of the segment offsets to points which lie toward the centre
     * point of the corner. The joining curve will not appear in the final
     * buffer outline, since it is completely internal to the buffer polygon.
     *
     * In complex buffer cases the closing segment may cut across many other
     * segments in the generated offset curve. In order to improve the
     * performance of the noding, the closing segment should be kept as short as
     * possible. (But not too short, since that would defeat its purpose). This
     * is the purpose of the closingSegFactor heuristic value.
     */

    /**
     * The intersection test above is vulnerable to robustness errors; i.e. it
     * may be that the offsets should intersect very close to their endpoints,
     * but aren't reported as such due to rounding. To handle this situation
     * appropriately, we use the following test: If the offset points are very
     * close, don't add closing segments but simply use one of the offset points
     */
    this.hasNarrowConcaveAngle = true;
    // System.out.println("NARROW ANGLE - distance = " + distance);
    if (this.offset0.p1.distance(this.offset1.p0) < this.distance *
        jsts.operation.buffer.OffsetSegmentGenerator.INSIDE_TURN_VERTEX_SNAP_DISTANCE_FACTOR) {
      this.segList.addPt(this.offset0.p1);
    } else {
      // add endpoint of this segment offset
      this.segList.addPt(this.offset0.p1);

      /**
       * Add "closing segment" of required length.
       */
      if (this.closingSegLengthFactor > 0) {
        var mid0 = new jsts.geom.Coordinate((this.closingSegLengthFactor *
            this.offset0.p1.x + this.s1.x) /
            (this.closingSegLengthFactor + 1), (this.closingSegLengthFactor *
            this.offset0.p1.y + this.s1.y) /
            (this.closingSegLengthFactor + 1));
        this.segList.addPt(mid0);
        var mid1 = new jsts.geom.Coordinate((this.closingSegLengthFactor *
            this.offset1.p0.x + this.s1.x) /
            (this.closingSegLengthFactor + 1), (this.closingSegLengthFactor *
            this.offset1.p0.y + this.s1.y) /
            (this.closingSegLengthFactor + 1));
        this.segList.addPt(mid1);
      } else {
        /**
         * This branch is not expected to be used except for testing purposes.
         * It is equivalent to the JTS 1.9 logic for closing segments (which
         * results in very poor performance for large buffer distances)
         */
        this.segList.addPt(this.s1);
      }

      // */
      // add start point of next segment offset
      this.segList.addPt(this.offset1.p0);
    }
  }
};


/**
 * Compute an offset segment for an input segment on a given side and at a given
 * distance. The offset points are computed in full double precision, for
 * accuracy.
 *
 * @param seg
 *          the segment to offset.
 * @param side
 *          the side of the segment ( {@link Position} ) the offset lies on.
 * @param distance
 *          the offset distance.
 * @param offset
 *          the points computed for the offset segment.
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.computeOffsetSegment = function(
    seg, side, distance, offset) {
  var sideSign = side === jsts.geomgraph.Position.LEFT ? 1 : -1;
  var dx = seg.p1.x - seg.p0.x;
  var dy = seg.p1.y - seg.p0.y;
  var len = Math.sqrt(dx * dx + dy * dy);
  // u is the vector that is the length of the offset, in the direction of the
  // segment
  var ux = sideSign * distance * dx / len;
  var uy = sideSign * distance * dy / len;
  offset.p0.x = seg.p0.x - uy;
  offset.p0.y = seg.p0.y + ux;
  offset.p1.x = seg.p1.x - uy;
  offset.p1.y = seg.p1.y + ux;
};


/**
 * Add an end cap around point p1, terminating a line segment coming from p0
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.addLineEndCap = function(
    p0, p1) {
  var seg = new jsts.geom.LineSegment(p0, p1);

  var offsetL = new jsts.geom.LineSegment();
  this.computeOffsetSegment(seg, jsts.geomgraph.Position.LEFT, this.distance,
      offsetL);
  var offsetR = new jsts.geom.LineSegment();
  this.computeOffsetSegment(seg, jsts.geomgraph.Position.RIGHT, this.distance,
      offsetR);

  var dx = p1.x - p0.x;
  var dy = p1.y - p0.y;
  var angle = Math.atan2(dy, dx);

  switch (this.bufParams.getEndCapStyle()) {
    case jsts.operation.buffer.BufferParameters.CAP_ROUND:
      // add offset seg points with a fillet between them
      this.segList.addPt(offsetL.p1);
      this.addFillet(p1, angle + Math.PI / 2, angle - Math.PI / 2,
          jsts.algorithm.CGAlgorithms.CLOCKWISE, this.distance);
      this.segList.addPt(offsetR.p1);
      break;
    case jsts.operation.buffer.BufferParameters.CAP_FLAT:
      // only offset segment points are added
      this.segList.addPt(offsetL.p1);
      this.segList.addPt(offsetR.p1);
      break;
    case jsts.operation.buffer.BufferParameters.CAP_SQUARE:
      // add a square defined by extensions of the offset segment endpoints
      var squareCapSideOffset = new jsts.geom.Coordinate();
      squareCapSideOffset.x = Math.abs(this.distance) * Math.cos(angle);
      squareCapSideOffset.y = Math.abs(this.distance) * Math.sin(angle);

      var squareCapLOffset = new jsts.geom.Coordinate(offsetL.p1.x +
          squareCapSideOffset.x, offsetL.p1.y + squareCapSideOffset.y);
      var squareCapROffset = new jsts.geom.Coordinate(offsetR.p1.x +
          squareCapSideOffset.x, offsetR.p1.y + squareCapSideOffset.y);
      this.segList.addPt(squareCapLOffset);
      this.segList.addPt(squareCapROffset);
      break;

  }
};


/**
 * Adds a mitre join connecting the two reflex offset segments. The mitre will
 * be beveled if it exceeds the mitre ratio limit.
 *
 * @param offset0
 *          the first offset segment.
 * @param offset1
 *          the second offset segment.
 * @param distance
 *          the offset distance.
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.addMitreJoin = function(
    p, offset0, offset1, distance) {
  var isMitreWithinLimit = true;
  var intPt = null;

  /**
   * This computation is unstable if the offset segments are nearly collinear.
   * Howver, this situation should have been eliminated earlier by the check for
   * whether the offset segment endpoints are almost coincident
   */
  try {
    intPt = jsts.algorithm.HCoordinate.intersection(offset0.p0, offset0.p1,
        offset1.p0, offset1.p1);

    var mitreRatio = distance <= 0.0 ? 1.0 : intPt.distance(p) /
        Math.abs(distance);

    if (mitreRatio > this.bufParams.getMitreLimit())
      this.isMitreWithinLimit = false;
  } catch (e) {
    if (e instanceof jsts.error.NotRepresentableError) {
      intPt = new jsts.geom.Coordinate(0, 0);
      this.isMitreWithinLimit = false;
    }
  }

  if (isMitreWithinLimit) {
    this.segList.addPt(intPt);
  } else {
    this.addLimitedMitreJoin(offset0, offset1, distance, bufParams
        .getMitreLimit());
    // addBevelJoin(offset0, offset1);
  }
};


/**
 * Adds a limited mitre join connecting the two reflex offset segments. A
 * limited mitre is a mitre which is beveled at the distance determined by the
 * mitre ratio limit.
 *
 * @param offset0
 *          the first offset segment.
 * @param offset1
 *          the second offset segment.
 * @param distance
 *          the offset distance.
 * @param mitreLimit
 *          the mitre limit ratio.
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.addLimitedMitreJoin = function(
    offset0, offset1, distance, mitreLimit) {
  var basePt = this.seg0.p1;

  var ang0 = jsts.algorithm.Angle.angle(basePt, this.seg0.p0);
  var ang1 = jsts.algorithm.Angle.angle(basePt, this.seg1.p1);

  // oriented angle between segments
  var angDiff = jsts.algorithm.Angle.angleBetweenOriented(this.seg0.p0, basePt,
      this.seg1.p1);
  // half of the interior angle
  var angDiffHalf = angDiff / 2;

  // angle for bisector of the interior angle between the segments
  var midAng = jsts.algorithm.Angle.normalize(ang0 + angDiffHalf);
  // rotating this by PI gives the bisector of the reflex angle
  var mitreMidAng = jsts.algorithm.Angle.normalize(midAng + Math.PI);

  // the miterLimit determines the distance to the mitre bevel
  var mitreDist = mitreLimit * distance;
  // the bevel delta is the difference between the buffer distance
  // and half of the length of the bevel segment
  var bevelDelta = mitreDist * Math.abs(Math.sin(angDiffHalf));
  var bevelHalfLen = distance - bevelDelta;

  // compute the midpoint of the bevel segment
  var bevelMidX = basePt.x + mitreDist * Math.cos(mitreMidAng);
  var bevelMidY = basePt.y + mitreDist * Math.sin(mitreMidAng);
  var bevelMidPt = new jsts.geom.Coordinate(bevelMidX, bevelMidY);

  // compute the mitre midline segment from the corner point to the bevel
  // segment midpoint
  var mitreMidLine = new jsts.geom.LineSegment(basePt, bevelMidPt);

  // finally the bevel segment endpoints are computed as offsets from
  // the mitre midline
  var bevelEndLeft = mitreMidLine.pointAlongOffset(1.0, bevelHalfLen);
  var bevelEndRight = mitreMidLine.pointAlongOffset(1.0, -bevelHalfLen);

  if (this.side == jsts.geomgraph.Position.LEFT) {
    this.segList.addPt(bevelEndLeft);
    this.segList.addPt(bevelEndRight);
  } else {
    this.segList.addPt(bevelEndRight);
    this.segList.addPt(bevelEndLeft);
  }
};


/**
 * Adds a bevel join connecting the two offset segments around a reflex corner.
 *
 * @param {LineSegment}
 *          offset0 the first offset segment.
 * @param {LineSegment}
 *          offset1 the second offset segment.
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.addBevelJoin = function(
    offset0, offset1) {
  this.segList.addPt(offset0.p1);
  this.segList.addPt(offset1.p0);
};


/**
 * Add points for a circular fillet around a reflex corner. Adds the start and
 * end points
 *
 * @param p
 *          base point of curve.
 * @param p0
 *          start point of fillet curve.
 * @param p1
 *          endpoint of fillet curve.
 * @param direction
 *          the orientation of the fillet.
 * @param radius
 *          the radius of the fillet.
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.addFillet = function(p,
    p0, p1, direction, radius) {
  if (!(p1 instanceof jsts.geom.Coordinate)) {
    this.addFillet2.apply(this, arguments);
    return;
  }

  var dx0 = p0.x - p.x;
  var dy0 = p0.y - p.y;
  var startAngle = Math.atan2(dy0, dx0);
  var dx1 = p1.x - p.x;
  var dy1 = p1.y - p.y;
  var endAngle = Math.atan2(dy1, dx1);

  if (direction === jsts.algorithm.CGAlgorithms.CLOCKWISE) {
    if (startAngle <= endAngle)
      startAngle += 2.0 * Math.PI;
  } else { // direction == COUNTERCLOCKWISE
    if (startAngle >= endAngle)
      startAngle -= 2.0 * Math.PI;
  }
  this.segList.addPt(p0);
  this.addFillet(p, startAngle, endAngle, direction, radius);
  this.segList.addPt(p1);
};


/**
 * Adds points for a circular fillet arc between two specified angles. The start
 * and end point for the fillet are not added - the caller must add them if
 * required.
 *
 * @param direction
 *          is -1 for a CW angle, 1 for a CCW angle.
 * @param radius
 *          the radius of the fillet.
 * @private
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.addFillet2 = function(p,
    startAngle, endAngle, direction, radius) {
  var directionFactor = direction === jsts.algorithm.CGAlgorithms.CLOCKWISE ? -1
      : 1;

  var totalAngle = Math.abs(startAngle - endAngle);
  var nSegs = parseInt((totalAngle / this.filletAngleQuantum + 0.5));

  if (nSegs < 1)
    return; // no segments because angle is less than increment - nothing to do!

  var initAngle, currAngleInc;

  // choose angle increment so that each segment has equal length
  initAngle = 0.0;
  currAngleInc = totalAngle / nSegs;

  var currAngle = initAngle;
  var pt = new jsts.geom.Coordinate();
  while (currAngle < totalAngle) {
    var angle = startAngle + directionFactor * currAngle;
    pt.x = p.x + radius * Math.cos(angle);
    pt.y = p.y + radius * Math.sin(angle);
    this.segList.addPt(pt);
    currAngle += currAngleInc;
  }
};


/**
 * Creates a CW circle around a point
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.createCircle = function(
    p) {
  // add start point
  var pt = new jsts.geom.Coordinate(p.x + this.distance, p.y);
  this.segList.addPt(pt);
  this.addFillet(p, 0.0, 2.0 * Math.PI, -1, this.distance);
  this.segList.closeRing();
};


/**
 * Creates a CW square around a point
 */
jsts.operation.buffer.OffsetSegmentGenerator.prototype.createSquare = function(
    p) {
  this.segList.addPt(new jsts.geom.Coordinate(p.x + distance, p.y + distance));
  this.segList.addPt(new jsts.geom.Coordinate(p.x + distance, p.y - distance));
  this.segList.addPt(new jsts.geom.Coordinate(p.x - distance, p.y - distance));
  this.segList.addPt(new jsts.geom.Coordinate(p.x - distance, p.y + distance));
  this.segList.closeRing();
};



//    JSTS LineSegment

/* Copyright (c) 2011 by The Authors.
 * Published under the LGPL 2.1 license.
 * See /license-notice.txt for the full text of the license notice.
 * See /license.txt for the full text of the license.
 */

/**
 * @requires jsts/geom/Coordinate.js
 * @requires jsts/algorithm/CGAlgorithms.js
 * @requires jsts/algorithm/RobustLineIntersector.js
 * @requires jsts/algorithm/HCoordinate.js
 */

/**
 * Represents a line segment defined by two {@link Coordinate}s. Provides
 * methods to compute various geometric properties and relationships of line
 * segments.
 * <p>
 * This class is designed to be easily mutable (to the extent of having its
 * contained points public). This supports a common pattern of reusing a single
 * LineSegment object as a way of computing segment properties on the segments
 * defined by arrays or lists of {@link Coordinate}s.
 *
 * @param {Coordinate}
 *          p0
 * @param {Coordinate}
 *          p1
 * @constructor
 */
jsts.geom.LineSegment = function () {
    if (arguments.length === 0) {
        this.p0 = new jsts.geom.Coordinate();
        this.p1 = new jsts.geom.Coordinate();
    } else if (arguments.length === 1) {
        this.p0 = arguments[0].p0;
        this.p1 = arguments[0].p1;
    } else if (arguments.length === 2) {
        this.p0 = arguments[0];
        this.p1 = arguments[1];
    } else if (arguments.length === 4) {
        this.p0 = new jsts.geom.Coordinate(arguments[0], arguments[1]);
        this.p1 = new jsts.geom.Coordinate(arguments[2], arguments[3]);
    }
};

/**
 * @type {Coordinate}
 */
jsts.geom.LineSegment.prototype.p0 = null;


/**
 * @type {Coordinate}
 */
jsts.geom.LineSegment.prototype.p1 = null;

/**
 * Computes the midpoint of a segment
 *
 * @param {jsts.geom.Coordinate} p0
 * @param {jsts.geom.Coordinate} p1
 * @return {jsts.geom.Coordinate} the midpoint of the segment
 */
jsts.geom.LineSegment.midPoint = function (p0, p1) {
    return new jsts.geom.Coordinate((p0.x + p1.x) / 2, (p0.y + p1.y) / 2);
};

/**
 * @param {number} i
 * @return {jsts.geom.Coordinate}
 */
jsts.geom.LineSegment.prototype.getCoordinate = function (i) {
    if (i === 0) return this.p0;
    return this.p1;
};

/**
 * Computes the length of the line segment.
 *
 * @return {number} the length of the line segment.
 */
jsts.geom.LineSegment.prototype.getLength = function () {
    return this.p0.distance(this.p1);
};

/**
 * Tests whether the segment is horizontal.
 *
 * @return {boolean} <code>true</code> if the segment is horizontal.
 */
jsts.geom.LineSegment.prototype.isHorizontal = function () {
    return this.p0.y === this.p1.y;
};
/**
 * Tests whether the segment is vertical.
 *
 * @return {boolean} <code>true</code> if the segment is vertical.
 */
jsts.geom.LineSegment.prototype.isVertical = function () {
    return this.p0.x === this.p1.x;
};

jsts.geom.LineSegment.prototype.orientationIndex = function (arg) {
    if (arg instanceof jsts.geom.LineSegment) {
        return this.orientationIndex1(arg);
    } else if (arg instanceof jsts.geom.Coordinate) {
        return this.orientationIndex2(arg);
    }
};

/**
  * Determines the orientation of a LineSegment relative to this segment.
  * The concept of orientation is specified as follows:
  * Given two line segments A and L,
  * <ul>
  * <li>A is to the left of a segment L if A lies wholly in the
  * closed half-plane lying to the left of L
  * <li>A is to the right of a segment L if A lies wholly in the
  * closed half-plane lying to the right of L
  * <li>otherwise, A has indeterminate orientation relative to L. This
  * happens if A is collinear with L or if A crosses the line determined by L.
  * </ul>
  *
  * @param {jsts.geom.LineSegment} seg the LineSegment to compare
  *
  * @return 1 if <code>seg</code> is to the left of this segment<br />
  * -1 if <code>seg</code> is to the right of this segment<br />
  * 0 if <code>seg</code> has indeterminate orientation relative to this segment
  */
jsts.geom.LineSegment.prototype.orientationIndex1 = function (seg) {
    var orient0 = jsts.algorithm.CGAlgorithms.orientationIndex(this.p0, this.p1, seg.p0);
    var orient1 = jsts.algorithm.CGAlgorithms.orientationIndex(this.p0, this.p1, seg.p1);
    // this handles the case where the points are L or collinear
    if (orient0 >= 0 && orient1 >= 0) {
        return Math.max(orient0, orient1);
    }
    // this handles the case where the points are R or collinear
    if (orient0 <= 0 && orient1 <= 0) {
        return Math.max(orient0, orient1);
    }
    // points lie on opposite sides ==> indeterminate orientation
    return 0;
};

/**
 * Determines the orientation index of a {@link Coordinate} relative to this segment.
 * The orientation index is as defined in {@link CGAlgorithms#computeOrientation}.
 *
 * @param {jsts.geom.Coordinate} p the coordinate to compare
 *
 * @return 1 (LEFT) if <code>p</code> is to the left of this segment
 * @return -1 (RIGHT) if <code>p</code> is to the right of this segment
 * @return 0 (COLLINEAR) if <code>p</code> is collinear with this segment
 * 
 * @see CGAlgorithms#computeOrientation(Coordinate, Coordinate, Coordinate)
 */
jsts.geom.LineSegment.prototype.orientationIndex2 = function (p) {
    return jsts.algorithm.CGAlgorithms.orientationIndex(this.p0, this.p1, p);
};

/**
 * Reverses the direction of the line segment.
 */
jsts.geom.LineSegment.prototype.reverse = function () {
    var temp = this.p0;
    this.p0 = this.p1;
    this.p1 = temp;
};

/**
 * Puts the line segment into a normalized form.
 * This is useful for using line segments in maps and indexes when
 * topological equality rather than exact equality is desired.
 * A segment in normalized form has the first point smaller
 * than the second (according to the standard ordering on {@link Coordinate}).
 */
jsts.geom.LineSegment.prototype.normalize = function () {
    if (this.p1.compareTo(this.p0) < 0) this.reverse();
};

/**
 * Computes the angle that the vector defined by this segment
 * makes with the X-axis.
 * The angle will be in the range [ -PI, PI ] radians.
 *
 * @return {number} the angle this segment makes with the X-axis (in radians)
 */
jsts.geom.LineSegment.prototype.angle = function () {
    return Math.atan2(this.p1.y - this.p0.y, this.p1.x - this.p0.x);
};

/**
 * Computes the midpoint of the segment
 *
 * @return {jsts.geom.Coordinate} the midpoint of the segment
 */
jsts.geom.LineSegment.prototype.midPoint = function () {
    return jsts.geom.LineSegment.midPoint(this.p0, this.p1);
};

jsts.geom.LineSegment.prototype.distance = function (arg) {
    if (arg instanceof jsts.geom.LineSegment) {
        return this.distance1(arg);
    } else if (arg instanceof jsts.geom.Coordinate) {
        return this.distance2(arg);
    }
};

/**
 * Computes the distance between this line segment and another segment.
 *
 * @param {jsts.geom.LineSegment} ls
 * @return {number} the distance to the other segment
 */
jsts.geom.LineSegment.prototype.distance1 = function (ls) {
    return jsts.algorithm.CGAlgorithms.distanceLineLine(this.p0, this.p1, ls.p0, ls.p1);
};

/**
 * Computes the distance between this line segment and a given point.
 *
 * @param {jsts.geom.Coordinate}
 *          p the coordinate.
 * @return {number}
 *          the distance from this segment to the given point.
 */
jsts.geom.LineSegment.prototype.distance2 = function (p) {
    return jsts.algorithm.CGAlgorithms.distancePointLine(p, this.p0, this.p1);
};

/**
 * Computes the {@link Coordinate} that lies a given
 * fraction along the line defined by this segment.
 * A fraction of <code>0.0</code> returns the start point of the segment;
 * a fraction of <code>1.0</code> returns the end point of the segment.
 * If the fraction is < 0.0 or > 1.0 the point returned 
 * will lie before the start or beyond the end of the segment. 
 *
 * @param {number} segmentLengthFraction the fraction of the segment length along the line
 * @return {jsts.geom.Coordinate} the point at that distance
 */
jsts.geom.LineSegment.prototype.pointAlong = function (segmentLengthFraction) {
    var coord = new jsts.geom.Coordinate();
    coord.x = this.p0.x + segmentLengthFraction * (this.p1.x - this.p0.x);
    coord.y = this.p0.y + segmentLengthFraction * (this.p1.y - this.p0.y);
    return coord;
};

/**
 * Computes the {@link Coordinate} that lies a given
 * fraction along the line defined by this segment and offset from 
 * the segment by a given distance.
 * A fraction of <code>0.0</code> offsets from the start point of the segment;
 * a fraction of <code>1.0</code> offsets from the end point of the segment.
 * The computed point is offset to the left of the line if the offset distance is
 * positive, to the right if negative.
 *
 * @param {number} segmentLengthFraction the fraction of the segment length along the line
 * @param {number} offsetDistance the distance the point is offset from the segment
 *    (positive is to the left, negative is to the right)
 * @return {jsts.geom.Coordinate} the point at that distance and offset
 */
jsts.geom.LineSegment.prototype.pointAlongOffset = function (segmentLengthFraction, offsetDistance) {
    // the point on the segment line
    var segx = this.p0.x + segmentLengthFraction * (this.p1.x - this.p0.x);
    var segy = this.p0.y + segmentLengthFraction * (this.p1.y - this.p0.y);

    var dx = this.p1.x - this.p0.x;
    var dy = this.p1.y - this.p0.y;
    var len = Math.sqrt(dx * dx + dy * dy);
    var ux = 0;
    var uy = 0;
    if (offsetDistance !== 0) {
        if (len <= 0) {
            throw "Cannot compute offset from zero-length line segment";
        }

        // u is the vector that is the length of the offset, in the direction of the segment
        ux = offsetDistance * dx / len;
        uy = offsetDistance * dy / len;
    }

    // the offset point is the seg point plus the offset vector rotated 90 degrees CCW
    var offsetx = segx - uy;
    var offsety = segy + ux;

    var coord = new jsts.geom.Coordinate(offsetx, offsety);
    return coord;
};

/**
 * Computes the Projection Factor for the projection of the point p onto this
 * LineSegment. The Projection Factor is the constant r by which the vector for
 * this segment must be multiplied to equal the vector for the projection of
 * <tt>p<//t> on the line
 * defined by this segment.
 * <p>
 * The projection factor returned will be in the range <tt>(-inf, +inf)</tt>.
 *
 * @param {Coordinate} p the point to compute the factor for.
 * @return {double} the projection factor for the point.
 */
jsts.geom.LineSegment.prototype.projectionFactor = function (p) {
    if (p.equals(this.p0))
        return 0.0;
    if (p.equals(this.p1))
        return 1.0;
    // Otherwise, use comp.graphics.algorithms Frequently Asked Questions method
    /*            AC dot AB
                   r = ---------
                         ||AB||^2
                r has the following meaning:
                r=0 P = A
                r=1 P = B
                r<0 P is on the backward extension of AB
                r>1 P is on the forward extension of AB
                0<r<1 P is interior to AB
        */
    var dx = this.p1.x - this.p0.x;
    var dy = this.p1.y - this.p0.y;
    var len2 = dx * dx + dy * dy;
    var r = ((p.x - this.p0.x) * dx + (p.y - this.p0.y) * dy) / len2;
    return r;
};

/**
 * Computes the fraction of distance (in <tt>[0.0, 1.0]</tt>) 
 * that the projection of a point occurs along this line segment.
 * If the point is beyond either ends of the line segment,
 * the closest fractional value (<tt>0.0</tt> or <tt>1.0</tt>) is returned.
 * <p>
 * Essentially, this is the {@link #projectionFactor} clamped to 
 * the range <tt>[0.0, 1.0]</tt>.
 * If the segment has zero length, 1.0 is returned.
 *  
 * @param {jsts.geom.Coordinate} inputPt the point
 * @return {number} the fraction along the line segment the projection of the point occurs
 */
jsts.geom.LineSegment.prototype.segmentFraction = function (inputPt) {
    var segFrac = this.projectionFactor(inputPt);
    if (segFrac < 0) {
        segFrac = 0;
    } else if (segFrac > 1 || isNaN(segFrac)) {
        segFrac = 1;
    }
    return segFrac;
};

jsts.geom.LineSegment.prototype.project = function (arg) {
    if (arg instanceof jsts.geom.Coordinate) {
        return this.project1(arg);
    } else if (arg instanceof jsts.geom.LineSegment) {
        return this.project2(arg);
    }
};

/**
 * Compute the projection of a point onto the line determined
 * by this line segment.
 * <p>
 * Note that the projected point
 * may lie outside the line segment.  If this is the case,
 * the projection factor will lie outside the range [0.0, 1.0].
 * @param {jsts.geom.Coordinate} p
 * @return {jsts.geom.Coordinate}
 */
jsts.geom.LineSegment.prototype.project1 = function (p) {
    if (p.equals(this.p0) || p.equals(this.p1)) {
        return new jsts.geom.Coordinate(p);
    }

    var r = this.projectionFactor(p);
    var coord = new jsts.geom.Coordinate();
    coord.x = this.p0.x + r * (this.p1.x - this.p0.x);
    coord.y = this.p0.y + r * (this.p1.y - this.p0.y);
    return coord;
};

/**
 * Project a line segment onto this line segment and return the resulting
 * line segment.  The returned line segment will be a subset of
 * the target line line segment.  This subset may be null, if
 * the segments are oriented in such a way that there is no projection.
 * <p>
 * Note that the returned line may have zero length (i.e. the same endpoints).
 * This can happen for instance if the lines are perpendicular to one another.
 *
 * @param {jsts.geom.LineSegment} seg the line segment to project
 * @return {jsts.geom.LineSegment} the projected line segment, or <code>null</code> if there is no overlap
 */
jsts.geom.LineSegment.prototype.project2 = function (seg) {
    var pf0 = this.projectionFactor(seg.p0);
    var pf1 = this.projectionFactor(seg.p1);
    // check if segment projects at all
    if (pf0 >= 1 && pf1 >= 1) return null;
    if (pf0 <= 0 && pf1 <= 0) return null;

    var newp0 = this.project(seg.p0);
    if (pf0 < 0) newp0 = p0;
    if (pf0 > 1) newp0 = p1;

    var newp1 = this.project(seg.p1);
    if (pf1 < 0.0) newp1 = p0;
    if (pf1 > 1.0) newp1 = p1;

    return new jsts.geom.LineSegment(newp0, newp1);
};

/**
 * Computes the closest point on this line segment to another point.
 *
 * @param {Coordinate}
 *          p the point to find the closest point to.
 * @return {Coordinate} a Coordinate which is the closest point on the line
 *         segment to the point p.
 */
jsts.geom.LineSegment.prototype.closestPoint = function (p) {
    var factor = this.projectionFactor(p);
    if (factor > 0 && factor < 1) {
        return this.project(p);
    }
    var dist0 = this.p0.distance(p);
    var dist1 = this.p1.distance(p);
    if (dist0 < dist1)
        return this.p0;
    return this.p1;
};


/**
 * Computes the closest points on two line segments.
 *
 * @param {LineSegment}
 *          line the segment to find the closest point to.
 * @return {[]} a pair of Coordinates which are the closest points on the line
 *         segments.
 */
jsts.geom.LineSegment.prototype.closestPoints = function (line) {
    // test for intersection
    var intPt = this.intersection(line);
    if (intPt !== null) {
        return [intPt, intPt];
    }

    /**
     * if no intersection closest pair contains at least one endpoint. Test each
     * endpoint in turn.
     */
    var closestPt = [];
    var minDistance = Number.MAX_VALUE;
    var dist;

    var close00 = this.closestPoint(line.p0);
    minDistance = close00.distance(line.p0);
    closestPt[0] = close00;
    closestPt[1] = line.p0;

    var close01 = this.closestPoint(line.p1);
    dist = close01.distance(line.p1);
    if (dist < minDistance) {
        minDistance = dist;
        closestPt[0] = close01;
        closestPt[1] = line.p1;
    }

    var close10 = line.closestPoint(this.p0);
    dist = close10.distance(this.p0);
    if (dist < minDistance) {
        minDistance = dist;
        closestPt[0] = this.p0;
        closestPt[1] = close10;
    }

    var close11 = line.closestPoint(this.p1);
    dist = close11.distance(this.p1);
    if (dist < minDistance) {
        minDistance = dist;
        closestPt[0] = this.p1;
        closestPt[1] = close11;
    }

    return closestPt;
};


/**
 * Computes an intersection point between two line segments, if there is one.
 * There may be 0, 1 or many intersection points between two segments. If there
 * are 0, null is returned. If there is 1 or more, exactly one of them is
 * returned (chosen at the discretion of the algorithm). If more information is
 * required about the details of the intersection, the
 * {@link RobustLineIntersector} class should be used.
 *
 * @param {LineSegment}
 *          line a line segment.
 * @return {Coordinate} an intersection point, or <code>null</code> if there
 *         is none.
 *
 * @see RobustLineIntersector
 */
jsts.geom.LineSegment.prototype.intersection = function (line) {
    var li = new jsts.algorithm.RobustLineIntersector();
    li.computeIntersection(this.p0, this.p1, line.p0, line.p1);
    if (li.hasIntersection())
        return li.getIntersection(0);
    return null;
};

jsts.geom.LineSegment.prototype.setCoordinates = function (ls) {
    if (ls instanceof jsts.geom.Coordinate) {
        this.setCoordinates2.apply(this, arguments);
        return;
    }

    this.setCoordinates2(ls.p0, ls.p1);
};

jsts.geom.LineSegment.prototype.setCoordinates2 = function (p0, p1) {
    this.p0.x = p0.x;
    this.p0.y = p0.y;
    this.p1.x = p1.x;
    this.p1.y = p1.y;
};

/**
 * Computes the perpendicular distance between the (infinite) line defined
 * by this line segment and a point.
 *
 * @param {jsts.geom.Coordinate} p the coordinate
 * @return {number} the perpendicular distance between the defined line and the given point
 */
jsts.geom.LineSegment.prototype.distancePerpendicular = function (p) {
    return jsts.algorithm.CGAlgorithms.distancePointLinePerpendicular(p, this.p0, this.p1);
};

/**
 * Computes the intersection point of the lines of infinite extent defined
 * by two line segments (if there is one).
 * There may be 0, 1 or an infinite number of intersection points 
 * between two lines.
 * If there is a unique intersection point, it is returned. 
 * Otherwise, <tt>null</tt> is returned.
 * If more information is required about the details of the intersection,
 * the {@link RobustLineIntersector} class should be used.
 *
 * @param {jsts.geom.LineSegment} line a line segment defining an straight line with infinite extent
 * @return {jsts.geom.Coordinate} an intersection point, 
 * or <code>null</code> if there is no point of intersection
 * or an infinite number of intersection points
 * 
 * @see RobustLineIntersector
 */
jsts.geom.LineSegment.prototype.lineIntersection = function (line) {
    try {
        var intPt = jsts.algorithm.HCoordinate.intersection(this.p0, this.p1, line.p0, line.p1);
        return intPt;
    } catch (ex) {
        // eat this exception, and return null;
    }
    return null;
};

/**
 * Creates a LineString with the same coordinates as this segment
 * 
 * @param {jsts.geom.GeometryFactory} geomFactory the geometery factory to use
 * @return {jsts.geom.LineString} a LineString with the same geometry as this segment
 */
jsts.geom.LineSegment.prototype.toGeometry = function (geomFactory) {
    return geomFactory.createLineString([this.p0, this.p1]);
};

/**
 *  Returns <code>true</code> if <code>other</code> has the same values for
 *  its points.
 *
 * @param {Object} o a <code>LineSegment</code> with which to do the comparison.
 * @return {boolean} <code>true</code> if <code>other</code> is a <code>LineSegment</code>
 *      with the same values for the x and y ordinates.
 */
jsts.geom.LineSegment.prototype.equals = function (o) {
    if (!(o instanceof jsts.geom.LineSegment)) {
        return false;
    }
    return this.p0.equals(o.p0) && this.p1.equals(o.p1);
};

/**
 *  Compares this object with the specified object for order.
 *  Uses the standard lexicographic ordering for the points in the LineSegment.
 *
 *@param {Object} o  the <code>LineSegment</code> with which this <code>LineSegment</code>
 *      is being compared
 *@return {number} a negative integer, zero, or a positive integer as this <code>LineSegment</code>
 *      is less than, equal to, or greater than the specified <code>LineSegment</code>
 */
jsts.geom.LineSegment.prototype.compareTo = function (o) {
    var comp0 = this.p0.compareTo(o.p0);
    if (comp0 !== 0) return comp0;
    return this.p1.compareTo(o.p1);
};

/**
 *  Returns <code>true</code> if <code>other</code> is
 *  topologically equal to this LineSegment (e.g. irrespective
 *  of orientation).
 *
 * @param {jsts.geom.LineSegment} other  a <code>LineSegment</code> with which to do the comparison.
 * @return {boolean} <code>true</code> if <code>other</code> is a <code>LineSegment</code>
 *      with the same values for the x and y ordinates.
 */
jsts.geom.LineSegment.prototype.equalsTopo = function (other) {
    return this.p0.equals(other.p0) && this.p1.equals(other.p1)
        || this.p0.equals(other.p1) && this.p1.equals(other.p0);
};

jsts.geom.LineSegment.prototype.toString = function () {
    return "LINESTRING(" +
        this.p0.x + " " + this.p0.y
        + ", " +
        this.p1.x + " " + this.p1.y + ")";
};





//    JSTS OffsetSegementString

/* Copyright (c) 2011 by The Authors.
 * Published under the LGPL 2.1 license.
 * See /license-notice.txt for the full text of the license notice.
 * See /license.txt for the full text of the license.
 */


/**
 * A dynamic list of the vertices in a constructed offset curve. Automatically
 * removes adjacent vertices which are closer than a given tolerance.
 * @constructor
 */
jsts.operation.buffer.OffsetSegmentString = function() {
  this.ptList = [];
};


/**
 * @private
 */
jsts.operation.buffer.OffsetSegmentString.prototype.ptList = null;


/**
 * @private
 */
jsts.operation.buffer.OffsetSegmentString.prototype.precisionModel = null;


/**
 * The distance below which two adjacent points on the curve are considered to
 * be coincident. This is chosen to be a small fraction of the offset distance.
 *
 * @private
 */
jsts.operation.buffer.OffsetSegmentString.prototype.minimimVertexDistance = 0.0;


jsts.operation.buffer.OffsetSegmentString.prototype.setPrecisionModel = function(
    precisionModel) {
  this.precisionModel = precisionModel;
};

jsts.operation.buffer.OffsetSegmentString.prototype.setMinimumVertexDistance = function(
    minimimVertexDistance) {
  this.minimimVertexDistance = minimimVertexDistance;
};

jsts.operation.buffer.OffsetSegmentString.prototype.addPt = function(pt) {
  var bufPt = new jsts.geom.Coordinate(pt);
  this.precisionModel.makePrecise(bufPt);
  // don't add duplicate (or near-duplicate) points
  if (this.isRedundant(bufPt))
    return;
  this.ptList.push(bufPt);
};

jsts.operation.buffer.OffsetSegmentString.prototype.addPts = function(pt,
    isForward) {
  if (isForward) {
    for (var i = 0; i < pt.length; i++) {
      this.addPt(pt[i]);
    }
  } else {
    for (var i = pt.length - 1; i >= 0; i--) {
      this.addPt(pt[i]);
    }
  }
};


/**
 * Tests whether the given point is redundant relative to the previous point in
 * the list (up to tolerance).
 *
 * @param pt
 * @return true if the point is redundant.
 * @private
 */
jsts.operation.buffer.OffsetSegmentString.prototype.isRedundant = function(pt) {
  if (this.ptList.length < 1)
    return false;
  var lastPt = this.ptList[this.ptList.length - 1];
  var ptDist = pt.distance(lastPt);
  if (ptDist < this.minimimVertexDistance)
    return true;
  return false;
};

jsts.operation.buffer.OffsetSegmentString.prototype.closeRing = function() {
  if (this.ptList.length < 1)
    return;
  var startPt = new jsts.geom.Coordinate(this.ptList[0]);
  var lastPt = this.ptList[this.ptList.length - 1];
  var last2Pt = null;
  if (this.ptList.length >= 2)
    last2Pt = this.ptList[this.ptList.length - 2];
  if (startPt.equals(lastPt))
    return;
  this.ptList.push(startPt);
};

jsts.operation.buffer.OffsetSegmentString.prototype.reverse = function() {

};

jsts.operation.buffer.OffsetSegmentString.prototype.getCoordinates = function() {
  return this.ptList;
};




//    JSTS CGAlgorithms

/* Copyright (c) 2011 by The Authors.
 * Published under the LGPL 2.1 license.
 * See /license-notice.txt for the full text of the license notice.
 * See /license.txt for the full text of the license.
 */



/**
 * Specifies and implements various fundamental Computational Geometric
 * algorithms. The algorithms supplied in this class are robust for
 * double-precision floating point.
 *
 * @constructor
 */
jsts.algorithm.CGAlgorithms = function() {

};


/**
 * A value that indicates an orientation of clockwise, or a right turn.
 */
jsts.algorithm.CGAlgorithms.CLOCKWISE = -1;


/**
 * A value that indicates an orientation of clockwise, or a right turn.
 */
jsts.algorithm.CGAlgorithms.RIGHT = jsts.algorithm.CGAlgorithms.CLOCKWISE;


/**
 * A value that indicates an orientation of counterclockwise, or a left turn.
 */
jsts.algorithm.CGAlgorithms.COUNTERCLOCKWISE = 1;


/**
 * A value that indicates an orientation of counterclockwise, or a left turn.
 */
jsts.algorithm.CGAlgorithms.LEFT = jsts.algorithm.CGAlgorithms.COUNTERCLOCKWISE;


/**
 * A value that indicates an orientation of collinear, or no turn (straight).
 */
jsts.algorithm.CGAlgorithms.COLLINEAR = 0;


/**
 * A value that indicates an orientation of collinear, or no turn (straight).
 */
jsts.algorithm.CGAlgorithms.STRAIGHT = jsts.algorithm.CGAlgorithms.COLLINEAR;


/**
 * Returns the index of the direction of the point <code>q</code> relative to
 * a vector specified by <code>p1-p2</code>.
 *
 * @param {jsts.geom.Coordinate}
 *          p1 the origin point of the vector.
 * @param {jsts.geom.Coordinate}
 *          p2 the final point of the vector.
 * @param {jsts.geom.Coordinate}
 *          q the point to compute the direction to.
 *
 * @return {Number} 1 if q is counter-clockwise (left) from p1-p2.
 * @return {Number} -1 if q is clockwise (right) from p1-p2.
 * @return {Number} 0 if q is collinear with p1-p2.
 */
jsts.algorithm.CGAlgorithms.orientationIndex = function(p1, p2, q) {
  /**
   * MD - 9 Aug 2010 It seems that the basic algorithm is slightly orientation
   * dependent, when computing the orientation of a point very close to a line.
   * This is possibly due to the arithmetic in the translation to the origin.
   *
   * For instance, the following situation produces identical results in spite
   * of the inverse orientation of the line segment:
   *
   * Coordinate p0 = new Coordinate(219.3649559090992, 140.84159161824724);
   * Coordinate p1 = new Coordinate(168.9018919682399, -5.713787599646864);
   *
   * Coordinate p = new Coordinate(186.80814046338352, 46.28973405831556); int
   * orient = orientationIndex(p0, p1, p); int orientInv = orientationIndex(p1,
   * p0, p);
   *
   * A way to force consistent results is to normalize the orientation of the
   * vector using the following code. However, this may make the results of
   * orientationIndex inconsistent through the triangle of points, so it's not
   * clear this is an appropriate patch.
   *
   */

  var dx1, dy1, dx2, dy2;
  dx1 = p2.x - p1.x;
  dy1 = p2.y - p1.y;
  dx2 = q.x - p2.x;
  dy2 = q.y - p2.y;

  return jsts.algorithm.RobustDeterminant.signOfDet2x2(dx1, dy1, dx2, dy2);
};


/**
 * Tests whether a point lies inside or on a ring. The ring may be oriented in
 * either direction. A point lying exactly on the ring boundary is considered to
 * be inside the ring.
 * <p>
 * This method does <i>not</i> first check the point against the envelope of
 * the ring.
 *
 * @param {jsts.geom.Coordinate}
 *          p point to check for ring inclusion.
 * @param {Array{jsts.geom.Coordinate}}
 *          ring an array of coordinates representing the ring (which must have
 *          first point identical to last point)
 * @return {Boolean} true if p is inside ring.
 *
 * @see locatePointInRing
 */
jsts.algorithm.CGAlgorithms.isPointInRing = function(p, ring) {
  return jsts.algorithm.CGAlgorithms.locatePointInRing(p, ring) !== jsts.geom.Location.EXTERIOR;
};


/**
 * Determines whether a point lies in the interior, on the boundary, or in the
 * exterior of a ring. The ring may be oriented in either direction.
 * <p>
 * This method does <i>not</i> first check the point against the envelope of
 * the ring.
 *
 * @param {jsts.geom.Coordinate}
 *          p point to check for ring inclusion.
 * @param {Array{jsts.geom.Coordinate}}
 *          ring an array of coordinates representing the ring (which must have
 *          first point identical to last point)
 * @return {jsts.geom.Location} the {@link Location} of p relative to the ring.
 */
jsts.algorithm.CGAlgorithms.locatePointInRing = function(p, ring) {
  return jsts.algorithm.RayCrossingCounter.locatePointInRing(p, ring);
};


/**
 * Tests whether a point lies on the line segments defined by a list of
 * coordinates.
 *
 * @param {jsts.geom.Coordinate}
 *          p the coordinate to test.
 * @param {Array{jsts.geom.Coordinate}}
 *          pt An array of coordinates defining line segments
 * @return {Boolean} true if the point is a vertex of the line or lies in the
 *         interior of a line segment in the linestring.
 */
jsts.algorithm.CGAlgorithms.isOnLine = function(p, pt) {
  var lineIntersector, i, il, p0, p1;
  lineIntersector = new jsts.algorithm.RobustLineIntersector();

  for (i = 1, il = pt.length; i < il; i++) {
    p0 = pt[i - 1];
    p1 = pt[i];
    lineIntersector.computeIntersection(p, p0, p1);

    if (lineIntersector.hasIntersection()) {
      return true;
    }
  }
  return false;
};


/**
 * Computes whether a ring defined by an array of {@link Coordinate}s is
 * oriented counter-clockwise.
 * <ul>
 * <li>The list of points is assumed to have the first and last points equal.
 * <li>This will handle coordinate lists which contain repeated points.
 * </ul>
 * This algorithm is <b>only</b> guaranteed to work with valid rings. If the
 * ring is invalid (e.g. self-crosses or touches), the computed result may not
 * be correct.
 *
 * @param {Array{jsts.geom.Coordinate}}
 *          ring an array of Coordinates forming a ring
 * @return {Boolean} true if the ring is oriented counter-clockwise.
 * @throws IllegalArgumentException
 *           if there are too few points to determine orientation (< 3)
 */
jsts.algorithm.CGAlgorithms.isCCW = function(ring) {
  var nPts, hiPt, hiIndex, p, iPrev, iNext, prev, next, i, disc, isCCW;

  // # of points without closing endpoint
  nPts = ring.length - 1;

  // sanity check
  if (nPts < 3) {
    throw new jsts.IllegalArgumentError(
        'Ring has fewer than 3 points, so orientation cannot be determined');
  }

  // find highets point
  hiPt = ring[0];
  hiIndex = 0;

  i = 1;
  for (i; i <= nPts; i++) {
    p = ring[i];
    if (p.y > hiPt.y) {
      hiPt = p;
      hiIndex = i;
    }
  }

  // find distinct point before highest point
  iPrev = hiIndex;
  do {
    iPrev = iPrev - 1;
    if (iPrev < 0) {
      iPrev = nPts;
    }
  } while (ring[iPrev].equals2D(hiPt) && iPrev !== hiIndex);

  // find distinct point after highest point
  iNext = hiIndex;
  do {
    iNext = (iNext + 1) % nPts;
  } while (ring[iNext].equals2D(hiPt) && iNext !== hiIndex);

  prev = ring[iPrev];
  next = ring[iNext];

  /**
   * This check catches cases where the ring contains an A-B-A configuration of
   * points. This can happen if the ring does not contain 3 distinct points
   * (including the case where the input array has fewer than 4 elements), or it
   * contains coincident line segments.
   */
  if (prev.equals2D(hiPt) || next.equals2D(hiPt) || prev.equals2D(next)) {
    return false;
  }

  disc = jsts.algorithm.CGAlgorithms.computeOrientation(prev, hiPt, next);

  /**
   * If disc is exactly 0, lines are collinear. There are two possible cases:
   * (1) the lines lie along the x axis in opposite directions (2) the lines lie
   * on top of one another
   *
   * (1) is handled by checking if next is left of prev ==> CCW (2) will never
   * happen if the ring is valid, so don't check for it (Might want to assert
   * this)
   */
  isCCW = false;
  if (disc === 0) {
    // poly is CCW if prev x is right of next x
    isCCW = (prev.x > next.x);
  } else {
    // if area is positive, points are ordered CCW
    isCCW = (disc > 0);
  }

  return isCCW;
};


/**
 * Computes the orientation of a point q to the directed line segment p1-p2. The
 * orientation of a point relative to a directed line segment indicates which
 * way you turn to get to q after travelling from p1 to p2.
 *
 * @param {jsts.geom.Coordinate}
 *          p1 First coordinate of the linesegment.
 * @param {jsts.geom.Coordinate}
 *          p2 Second coordinate of the linesegment.
 * @param {jsts.geom.Coordinate}
 *          q The point to calculate orientation of.
 *
 * @return {Number} 1 if q is counter-clockwise from p1-p2.
 * @return {Number} -1 if q is clockwise from p1-p2.
 * @return {Number} 0 if q is collinear with p1-p2.
 */
jsts.algorithm.CGAlgorithms.computeOrientation = function(p1, p2, q) {
  return jsts.algorithm.CGAlgorithms.orientationIndex(p1, p2, q);
};


/**
 * Computes the distance from a point p to a line segment AB
 *
 * Note: NON-ROBUST!
 *
 * @param {jsts.geom.Coordinate}
 *          p the point to compute the distance for.
 * @param {jsts.geom.Coordinate}
 *          A one point of the line.
 * @param {jsts.geom.Coordinate}
 *          B another point of the line (must be different to A).
 * @return {Number} the distance from p to line segment AB.
 */
jsts.algorithm.CGAlgorithms.distancePointLine = function(p, A, B) {
  if (!(A instanceof jsts.geom.Coordinate)) {
    jsts.algorithm.CGAlgorithms.distancePointLine2.apply(this, arguments);
  }

  // if start = end, then just compute distance to one of the endpoints
  if (A.x === B.x && A.y === B.y) {
    return p.distance(A);
  }
  // otherwise use comp.graphics.algorithms Frequently Asked Questions method
  /*(1)             AC dot AB
                   r = ---------
                         ||AB||^2
    r has the following meaning:
    r=0 P = A
    r=1 P = B
    r<0 P is on the backward extension of AB
    r>1 P is on the forward extension of AB
    0<r<1 P is interior to AB
  */
  var r, s;
  r = ((p.x - A.x) * (B.x - A.x) + (p.y - A.y) * (B.y - A.y)) /
      ((B.x - A.x) * (B.x - A.x) + (B.y - A.y) * (B.y - A.y));

  if (r <= 0.0) {
    return p.distance(A);
  }
  if (r >= 1.0) {
    return p.distance(B);
  }

  /*(2)
    (Ay-Cy)(Bx-Ax)-(Ax-Cx)(By-Ay)
  s = -----------------------------
             L^2

  Then the distance from C to P = |s|*L.
  */

  s = ((A.y - p.y) * (B.x - A.x) - (A.x - p.x) * (B.y - A.y)) /
      ((B.x - A.x) * (B.x - A.x) + (B.y - A.y) * (B.y - A.y));

  return Math.abs(s) *
      Math.sqrt(((B.x - A.x) * (B.x - A.x) + (B.y - A.y) * (B.y - A.y)));
};


/**
 * Computes the perpendicular distance from a point p to the (infinite) line
 * containing the points AB
 *
 * @param {jsts.geom.Coordinate}
 *          p the point to compute the distance for.
 * @param {jsts.geom.Coordinate}
 *          A one point of the line.
 * @param {jsts.geom.Coordinate}
 *          B another point of the line (must be different to A).
 * @return {Number} the distance from p to line AB.
 */
jsts.algorithm.CGAlgorithms.distancePointLinePerpendicular = function(p, A, B) {
  // use comp.graphics.algorithms Frequently Asked Questions method
  /*(2)
                   (Ay-Cy)(Bx-Ax)-(Ax-Cx)(By-Ay)
              s = -----------------------------
                                   L^2

              Then the distance from C to P = |s|*L.
  */
  var s = ((A.y - p.y) * (B.x - A.x) - (A.x - p.x) * (B.y - A.y)) /
      ((B.x - A.x) * (B.x - A.x) + (B.y - A.y) * (B.y - A.y));

  return Math.abs(s) *
      Math.sqrt(((B.x - A.x) * (B.x - A.x) + (B.y - A.y) * (B.y - A.y)));
};


/**
 * Computes the distance from a point to a sequence of line segments.
 *
 * @param {jsts.geom.Coordinate}
 *          p a point.
 * @param {Array{jsts.geom.Coordinate}}
 *          line a sequence of contiguous line segments defined by their
 *          vertices
 * @return {Number} the minimum distance between the point and the line
 *         segments.
 */
jsts.algorithm.CGAlgorithms.distancePointLine2 = function(p, line) {
  var minDistance, i, il, dist;
  if (line.length === 0) {
    throw new jsts.error.IllegalArgumentError(
        'Line array must contain at least one vertex');
  }
  minDistance = p.distance(line[0]);
  for (i = 0, il = line.length - 1; i < il; i++) {
    dist = jsts.algorithm.CGAlgorithms.distancePointLine(p, line[i],
        line[i + 1]);
    if (dist < minDistance) {
      minDistance = dist;
    }
  }
  return minDistance;
};

/**
 * Computes the distance from a line segment AB to a line segment CD
 *
 * Note: NON-ROBUST!
 *
 * @param {jsts.geom.Coordinate}
 *          A a point of one line.
 * @param {jsts.geom.Coordinate}
 *          B the second point of (must be different to A).
 * @param {jsts.geom.Coordinate}
 *          C one point of the line.
 * @param {jsts.geom.Coordinate}
 *          D another point of the line (must be different to A).
 * @return {Number} the distance.
 */

jsts.algorithm.CGAlgorithms.distanceLineLine = function(A, B, C, D) {
  // check for zero-length segments
  if (A.equals(B)) {
    return jsts.algorithm.CGAlgorithms.distancePointLine(A, C, D);
  }
  if (C.equals(D)) {
    return jsts.algorithm.CGAlgorithms.distancePointLine(D, A, B);
  }

  // AB and CD are line segments
  /* from comp.graphics.algo

  Solving the above for r and s yields
        (Ay-Cy)(Dx-Cx)-(Ax-Cx)(Dy-Cy)
             r = ----------------------------- (eqn 1)
        (Bx-Ax)(Dy-Cy)-(By-Ay)(Dx-Cx)

      (Ay-Cy)(Bx-Ax)-(Ax-Cx)(By-Ay)
    s = ----------------------------- (eqn 2)
      (Bx-Ax)(Dy-Cy)-(By-Ay)(Dx-Cx)
  Let P be the position vector of the intersection point, then
    P=A+r(B-A) or
    Px=Ax+r(Bx-Ax)
    Py=Ay+r(By-Ay)
  By examining the values of r & s, you can also determine some other
  limiting conditions:
    If 0<=r<=1 & 0<=s<=1, intersection exists
    r<0 or r>1 or s<0 or s>1 line segments do not intersect
    If the denominator in eqn 1 is zero, AB & CD are parallel
    If the numerator in eqn 1 is also zero, AB & CD are collinear.

  */
  var r_top, r_bot, s_top, s_bot, s, r;
  r_top = (A.y - C.y) * (D.x - C.x) - (A.x - C.x) * (D.y - C.y);
  r_bot = (B.x - A.x) * (D.y - C.y) - (B.y - A.y) * (D.x - C.x);

  s_top = (A.y - C.y) * (B.x - A.x) - (A.x - C.x) * (B.y - A.y);
  s_bot = (B.x - A.x) * (D.y - C.y) - (B.y - A.y) * (D.x - C.x);


  if ((r_bot === 0) || (s_bot === 0)) {
    return Math.min(jsts.algorithm.CGAlgorithms.distancePointLine(A, C, D),
        Math.min(jsts.algorithm.CGAlgorithms.distancePointLine(B, C, D), Math
            .min(jsts.algorithm.CGAlgorithms.distancePointLine(C, A, B),
                jsts.algorithm.CGAlgorithms.distancePointLine(D, A, B))));
  }

  s = s_top / s_bot;
  r = r_top / r_bot;
  if ((r < 0) || (r > 1) || (s < 0) || (s > 1)) {
    // no intersection
    return Math.min(jsts.algorithm.CGAlgorithms.distancePointLine(A, C, D),
        Math.min(jsts.algorithm.CGAlgorithms.distancePointLine(B, C, D), Math
            .min(jsts.algorithm.CGAlgorithms.distancePointLine(C, A, B),
                jsts.algorithm.CGAlgorithms.distancePointLine(D, A, B))));
  }

  return 0.0; // intersection exists
};


/**
 * Computes the signed area for a ring. The signed area is positive if the ring
 * is oriented CW, negative if the ring is oriented CCW, and zero if the ring is
 * degenerate or flat.
 *
 * @param {Array{jsts.geom.Coordinate}}
 *          ring the coordinates forming the ring
 * @return {Number} the signed area of the ring.
 */
jsts.algorithm.CGAlgorithms.signedArea = function(ring) {
  if (ring.length < 3) {
    return 0.0;
  }
  var sum, i, il, bx, by, cx, cy;

  sum = 0.0;

  for (i = 0, il = ring.length - 1; i < il; i++) {
    bx = ring[i].x;
    by = ring[i].y;
    cx = ring[i + 1].x;
    cy = ring[i + 1].y;
    sum += (bx + cx) * (cy - by);
  }

  return -sum / 2.0;
};


/**
 * Computes the signed area for a ring. The signed area is:
 * <ul>
 * <li>positive if the ring is oriented CW
 * <li>negative if the ring is oriented CCW
 * <li>zero if the ring is degenerate or flat
 * </ul>
 *
 * @param {Array{jsts.geom.Coordinate}}
 *          ring the coordinates forming the ring
 * @return {Number} the signed area of the ring.
 */
jsts.algorithm.CGAlgorithms.signedArea = function(ring) {
  var n, sum, p, bx, by, i, cx, cy;

  n = ring.length;
  if (n < 3) {
    return 0.0;
  }

  sum = 0.0;
  p = ring[0];

  bx = p.x;
  by = p.y;

  for (i = 1; i < n; i++) {
    p = ring[i];
    cx = p.x;
    cy = p.y;
    sum += (bx + cx) * (cy - by);
    bx = cx;
    by = cy;
  }

  return -sum / 2.0;
};

/**
 * Computes the length of a linestring specified by a sequence of points.
 *
 * NOTE: This is renamed from length() to computeLength() because 'length' is a
 * reserved keyword in javascript.
 *
 * @param {Array{jsts.geom.Coordinate}}
 *          pts the points specifying the linestring
 * @return {Number} the length of the linestring.
 */
jsts.algorithm.CGAlgorithms.computeLength = function(pts) {
  // optimized for processing CoordinateSequences
  var n = pts.length, len, x0, y0, x1, y1, dx, dy, p, i, il;
  if (n <= 1) {
    return 0.0;
  }

  len = 0.0;

  p = pts[0];

  x0 = p.x;
  y0 = p.y;

  i = 1, il = n;
  for (i; i < n; i++) {
    p = pts[i];

    x1 = p.x;
    y1 = p.y;
    dx = x1 - x0;
    dy = y1 - y0;

    len += Math.sqrt(dx * dx + dy * dy);

    x0 = x1;
    y0 = y1;
  }
  return len;
};

/**
 * @see {jsts.algorithm.CGAlgorithms.computeLength} Since 'length' is a reserved
 *      keyword in javascript this function does not act as a function. Please
 *      use 'computeLength' instead.
 */
jsts.algorithm.CGAlgorithms.length = function() {};





//    JSTS Location

/* Copyright (c) 2011 by The Authors.
 * Published under the LGPL 2.1 license.
 * See /license-notice.txt for the full text of the license notice.
 * See /license.txt for the full text of the license.
 */



/**
 * Constants representing the different topological locations which can occur in
 * a {@link Geometry}. The constants are also used as the row and column
 * indices of DE-9IM {@link IntersectionMatrix}es.
 *
 * @constructor
 */
jsts.geom.Location = function() {
};


/**
 * The location value for the interior of a geometry. Also, DE-9IM row index of
 * the interior of the first geometry and column index of the interior of the
 * second geometry.
 *
 * @const
 * @type {number}
 */
jsts.geom.Location.INTERIOR = 0;


/**
 * The location value for the boundary of a geometry. Also, DE-9IM row index of
 * the boundary of the first geometry and column index of the boundary of the
 * second geometry.
 *
 * @const
 * @type {number}
 */
jsts.geom.Location.BOUNDARY = 1;


/**
 * The location value for the exterior of a geometry. Also, DE-9IM row index of
 * the exterior of the first geometry and column index of the exterior of the
 * second geometry.
 *
 * @const
 * @type {number}
 */
jsts.geom.Location.EXTERIOR = 2;


/**
 * Used for uninitialized location values.
 *
 * @const
 * @type {number}
 */
jsts.geom.Location.NONE = -1;


/**
 * Converts the location value to a location symbol, for example,
 * <code>EXTERIOR => 'e'</code> .
 *
 * @param {number}
 *          locationValue either EXTERIOR, BOUNDARY, INTERIOR or NONE.
 * @return {string} either 'e', 'b', 'i' or '-'.
 */
jsts.geom.Location.toLocationSymbol = function(locationValue) {
  switch (locationValue) {
    case jsts.geom.Location.EXTERIOR:
      return 'e';
    case jsts.geom.Location.BOUNDARY:
      return 'b';
    case jsts.geom.Location.INTERIOR:
      return 'i';
    case jsts.geom.Location.NONE:
      return '-';
  }
  throw new jsts.IllegalArgumentError('Unknown location value: ' +
      locationValue);
};






//    JSTS Javascript Utils

(function(c){var d,b,a,h,e,f={".js":[],".json":[],".css":[],".html":[]};h=function(a){a=Error("Could not find module '"+a+"'");a.code="MODULE_NOT_FOUND";return a};e=function(a,l,b){var e,h;if("function"===typeof a[l+b])return l+b;for(e=0;h=f[b][e];++e)if("function"===typeof a[l+h])return l+h;return null};d=function(a,l,k,c,f){var g,n;k=k.split("/");g=k.pop();if("."===g||".."===g)k.push(g),g="";for(;null!=(n=k.shift());)if(n&&"."!==n&&(".."===n?a=l.pop():(l.push(a),a=a[n]),!a))throw h(c);g&&"function"!==
typeof a[g]&&((k=e(a,g,".js"))||(k=e(a,g,".json")),k||(k=e(a,g,".css")),k||(k=e(a,g,".html")),k?g=k:2!==f&&"object"===typeof a[g]&&(l.push(a),a=a[g],g=""));if(!g)return 1!==f&&a[":mainpath:"]?d(a,l,a[":mainpath:"],c,1):d(a,l,"index",c,2);f=a[g];if(!f)throw h(c);if(f.hasOwnProperty("module"))return f.module.exports;c={};f.module=g={exports:c};f.call(c,c,g,b(a,l));return g.exports};a=function(a,l,b){var e,f=b;e=b.charAt(0);var g=0;if("/"===e)f=f.slice(1),a=c["/"],l=[];else if("."!==e){e=f.split("/",
1)[0];a=c[e];if(!a)throw h(b);l=[];f=f.slice(e.length+1);f||((f=a[":mainpath:"])?g=1:(f="index",g=2))}return d(a,l,f,b,g)};b=function(b,l){return function(e){return a(b,[].concat(l),e)}};return b(c,[])})({"javascript.util":{lib:{"ArrayList.js":function(c,d,b){function a(a){this.array=[];a instanceof h&&this.addAll(a)}var h=b("./Collection");c=b("./List");var e=b("./OperationNotSupported"),f=b("./NoSuchElementException"),m=b("./IndexOutOfBoundsException");a.prototype=new c;a.prototype.array=null;a.prototype.add=
function(a){this.array.push(a);return!0};a.prototype.addAll=function(a){for(a=a.iterator();a.hasNext();)this.add(a.next());return!0};a.prototype.set=function(a,b){var e=this.array[a];this.array[a]=b;return e};a.prototype.iterator=function(){return new a.Iterator(this)};a.prototype.get=function(a){if(0>a||a>=this.size())throw new m;return this.array[a]};a.prototype.isEmpty=function(){return 0===this.array.length};a.prototype.size=function(){return this.array.length};a.prototype.toArray=function(){for(var a=
[],b=0,e=this.array.length;b<e;b++)a.push(this.array[b]);return a};a.prototype.remove=function(a){for(var b=!1,e=0,m=this.array.length;e<m;e++)if(this.array[e]===a){this.array.splice(e,1);b=!0;break}return b};a.Iterator=function(a){this.arrayList=a};a.Iterator.prototype.arrayList=null;a.Iterator.prototype.position=0;a.Iterator.prototype.next=function(){if(this.position===this.arrayList.size())throw new f;return this.arrayList.get(this.position++)};a.Iterator.prototype.hasNext=function(){return this.position<
this.arrayList.size()?!0:!1};a.Iterator.prototype.remove=function(){throw new e;};d.exports=a},"Arrays.js":function(c,d,b){function a(){}a.sort=function(){var a=arguments[0],b,f,m;if(1===arguments.length)a.sort();else if(2===arguments.length)f=arguments[1],m=function(a,b){return f.compare(a,b)},a.sort(m);else if(3===arguments.length)for(b=a.slice(arguments[1],arguments[2]),b.sort(),m=a.slice(0,arguments[1]).concat(b,a.slice(arguments[2],a.length)),a.splice(0,a.length),b=0;b<m.length;b++)a.push(m[b]);
else if(4===arguments.length)for(b=a.slice(arguments[1],arguments[2]),f=arguments[3],m=function(a,b){return f.compare(a,b)},b.sort(m),m=a.slice(0,arguments[1]).concat(b,a.slice(arguments[2],a.length)),a.splice(0,a.length),b=0;b<m.length;b++)a.push(m[b])};a.asList=function(a){for(var b=new javascript.util.ArrayList,f=0,m=a.length;f<m;f++)b.add(a[f]);return b};d.exports=a},"Collection.js":function(c,d,b){function a(){}b("./Iterator");a.prototype.add=function(a){};a.prototype.addAll=function(a){};a.prototype.isEmpty=
function(){};a.prototype.iterator=function(){};a.prototype.size=function(){};a.prototype.toArray=function(){};a.prototype.remove=function(a){};d.exports=a},"EmptyStackException.js":function(c,d,b){function a(a){this.message=a||""}a.prototype=Error();a.prototype.name="EmptyStackException";d.exports=a},"HashMap.js":function(c,d,b){function a(){this.object={}}c=b("./Map");var h=b("./ArrayList");a.prototype=new c;a.prototype.object=null;a.prototype.get=function(a){return this.object[a]||null};a.prototype.put=
function(a,b){return this.object[a]=b};a.prototype.values=function(){var a=new h,b;for(b in this.object)this.object.hasOwnProperty(b)&&a.add(this.object[b]);return a};a.prototype.size=function(){return this.values().size()};d.exports=a},"HashSet.js":function(c,d,b){function a(a){this.array=[];a instanceof h&&this.addAll(a)}var h=b("./Collection");c=b("./Set");var e=b("./OperationNotSupported"),f=b("./NoSuchElementException");a.prototype=new c;a.prototype.array=null;a.prototype.contains=function(a){for(var b=
0,e=this.array.length;b<e;b++)if(this.array[b]===a)return!0;return!1};a.prototype.add=function(a){if(this.contains(a))return!1;this.array.push(a);return!0};a.prototype.addAll=function(a){for(a=a.iterator();a.hasNext();)this.add(a.next());return!0};a.prototype.remove=function(a){throw new e;};a.prototype.size=function(){return this.array.length};a.prototype.isEmpty=function(){return 0===this.array.length};a.prototype.toArray=function(){for(var a=[],b=0,e=this.array.length;b<e;b++)a.push(this.array[b]);
return a};a.prototype.iterator=function(){return new a.Iterator(this)};a.Iterator=function(a){this.hashSet=a};a.Iterator.prototype.hashSet=null;a.Iterator.prototype.position=0;a.Iterator.prototype.next=function(){if(this.position===this.hashSet.size())throw new f;return this.hashSet.array[this.position++]};a.Iterator.prototype.hasNext=function(){return this.position<this.hashSet.size()?!0:!1};a.Iterator.prototype.remove=function(){throw new javascript.util.OperationNotSupported;};d.exports=a},"IndexOutOfBoundsException.js":function(c,
d,b){function a(a){this.message=a||""}a.prototype=Error();a.prototype.name="IndexOutOfBoundsException";d.exports=a},"Iterator.js":function(c,d,b){function a(){}a.prototype.hasNext=function(){};a.prototype.next=function(){};a.prototype.remove=function(){};d.exports=a},"List.js":function(c,d,b){function a(){}c=b("./Collection");a.prototype=new c;a.prototype.get=function(a){};a.prototype.set=function(a,b){};a.prototype.isEmpty=function(){};d.exports=a},"Map.js":function(c,d,b){function a(){}a.prototype.get=
function(a){};a.prototype.put=function(a,b){};a.prototype.size=function(){};a.prototype.values=function(){};d.exports=a},"NoSuchElementException.js":function(c,d,b){function a(a){this.message=a||""}a.prototype=Error();a.prototype.name="NoSuchElementException";d.exports=a},"OperationNotSupported.js":function(c,d,b){function a(a){this.message=a||""}a.prototype=Error();a.prototype.name="OperationNotSupported";d.exports=a},"Set.js":function(c,d,b){function a(){}c=b("./Collection");a.prototype=new c;a.prototype.contains=
function(a){};d.exports=a},"SortedMap.js":function(c,d,b){function a(){}c=b("./Map");a.prototype=new c;d.exports=a},"SortedSet.js":function(c,d,b){function a(){}c=b("./Set");a.prototype=new c;d.exports=a},"Stack.js":function(c,d,b){function a(){this.array=[]}c=b("./List");var h=b("./EmptyStackException");a.prototype=new c;a.prototype.array=null;a.prototype.push=function(a){this.array.push(a);return a};a.prototype.pop=function(a){if(0===this.array.length)throw new h;return this.array.pop()};a.prototype.peek=
function(){if(0===this.array.length)throw new h;return this.array[this.array.length-1]};a.prototype.empty=function(a){return 0===this.array.length?!0:!1};a.prototype.isEmpty=function(){return this.empty()};a.prototype.search=function(a){return this.array.indexOf(a)};a.prototype.size=function(){return this.array.length};a.prototype.toArray=function(){for(var a=[],b=0,c=this.array.length;b<c;b++)a.push(this.array[b]);return a};d.exports=a},"TreeMap.js":function(c,d,b){function a(){this.array=[]}c=b("./Map");
b("./SortedMap");var h=b("./ArrayList");a.prototype=new c;a.prototype.array=null;a.prototype.get=function(a){for(var b=0,c=this.array.length;b<c;b++){var d=this.array[b];if(0===d.key.compareTo(a))return d.value}return null};a.prototype.put=function(a,b){var c=this.get(a);if(c){var d=c.value;c.value=b;return d}for(var d={key:a,value:b},k=0,h=this.array.length;k<h;k++)if(c=this.array[k],1===c.key.compareTo(a))return this.array.splice(k,0,d),null;this.array.push({key:a,value:b});return null};a.prototype.values=
function(){for(var a=new h,b=0,c=this.array.length;b<c;b++)a.add(this.array[b].value);return a};a.prototype.size=function(){return this.values().size()};d.exports=a},"TreeSet.js":function(c,d,b){function a(a){this.array=[];a instanceof h&&this.addAll(a)}var h=b("./Collection");c=b("./SortedSet");var e=b("./OperationNotSupported"),f=b("./NoSuchElementException");a.prototype=new c;a.prototype.array=null;a.prototype.contains=function(a){for(var b=0,c=this.array.length;b<c;b++)if(0===this.array[b].compareTo(a))return!0;
return!1};a.prototype.add=function(a){if(this.contains(a))return!1;for(var b=0,c=this.array.length;b<c;b++)if(1===this.array[b].compareTo(a))return this.array.splice(b,0,a),!0;this.array.push(a);return!0};a.prototype.addAll=function(a){for(a=a.iterator();a.hasNext();)this.add(a.next());return!0};a.prototype.remove=function(a){throw new e;};a.prototype.size=function(){return this.array.length};a.prototype.isEmpty=function(){return 0===this.array.length};a.prototype.toArray=function(){for(var a=[],
b=0,c=this.array.length;b<c;b++)a.push(this.array[b]);return a};a.prototype.iterator=function(){return new a.Iterator(this)};a.Iterator=function(a){this.treeSet=a};a.Iterator.prototype.treeSet=null;a.Iterator.prototype.position=0;a.Iterator.prototype.next=function(){if(this.position===this.treeSet.size())throw new f;return this.treeSet.array[this.position++]};a.Iterator.prototype.hasNext=function(){return this.position<this.treeSet.size()?!0:!1};a.Iterator.prototype.remove=function(){throw new e;
};d.exports=a},"browser.js":function(c,d,b){javascript={util:b("./")}},"index.js":function(c,d,b){d.exports={ArrayList:b("./ArrayList"),Arrays:b("./Arrays"),Collection:b("./Collection"),HashMap:b("./HashMap"),HashSet:b("./HashSet"),Iterator:b("./Iterator"),List:b("./List"),Map:b("./Map"),Set:b("./Set"),SortedMap:b("./SortedMap"),SortedSet:b("./SortedSet"),Stack:b("./Stack"),TreeMap:b("./TreeMap"),TreeSet:b("./TreeSet")}}}}})("javascript.util/lib/browser");
