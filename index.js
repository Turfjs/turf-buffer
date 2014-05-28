//http://stackoverflow.com/questions/839899/how-do-i-calculate-a-point-on-a-circles-circumference
//radians = degrees * (pi/180)
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

var jsts = {}
jsts.io = {}
jsts.geom = {}

//  GEOJSON PARSER

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



// GEOJSON READER

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



//  JSTS GEOMETRY FACTORY

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



//  JSTS PrecisionModel

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



//   JSTS GEOMETRY

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


  //   JSTS Point

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




