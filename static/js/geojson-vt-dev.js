!function(e){if("object"==typeof exports&&"undefined"!=typeof module)module.exports=e();else if("function"==typeof define&&define.amd)define([],e);else{var o;"undefined"!=typeof window?o=window:"undefined"!=typeof global?o=global:"undefined"!=typeof self&&(o=self),o.geojsonvt=e()}}(function(){var define,module,exports;return (function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({1:[function(require,module,exports){
'use strict';

module.exports = clip;

/* clip features between two axis-parallel lines:
 *     |        |
 *  ___|___     |     /
 * /   |   \____|____/
 *     |        |
 */

function clip(features, scale, k1, k2, axis, intersect) {

    k1 /= scale;
    k2 /= scale;

    var clipped = [];

    for (var i = 0; i < features.length; i++) {

        var feature = features[i],
            geometry = feature.geometry,
            type = feature.type,
            min, max;

        if (feature.min) {
            min = feature.min[axis];
            max = feature.max[axis];

            if (min >= k1 && max <= k2) { // trivial accept
                clipped.push(feature);
                continue;
            } else if (min > k2 || max < k1) continue; // trivial reject
        }

        var slices = type === 1 ?
                clipPoints(geometry, k1, k2, axis) :
                clipGeometry(geometry, k1, k2, axis, intersect, type === 3);

        if (slices.length) {
            // if a feature got clipped, it will likely get clipped on the next zoom level as well,
            // so there's no need to recalculate bboxes
            clipped.push({
                geometry: slices,
                type: type,
                tags: features[i].tags || null
            });
        }
    }

    return clipped.length ? clipped : null;
}

function clipPoints(geometry, k1, k2, axis) {
    var slice = [];

    for (var i = 0; i < geometry.length; i++) {
        var a = geometry[i],
            ak = a[axis];

        if (ak >= k1 && ak <= k2) slice.push(a);
    }
    return slice;
}

function clipGeometry(geometry, k1, k2, axis, intersect, closed) {

    var slices = [];

    for (var i = 0; i < geometry.length; i++) {

        var ak = 0,
            bk = 0,
            b = null,
            points = geometry[i],
            area = points.area,
            dist = points.dist,
            len = points.length,
            a, j;

        var slice = [];

        for (j = 0; j < len - 1; j++) {
            a = b || points[j];
            b = points[j + 1];
            ak = bk || a[axis];
            bk = b[axis];

            if (ak < k1) {

                if ((bk > k2)) { // ---|-----|-->
                    slice.push(intersect(a, b, k1), intersect(a, b, k2));
                    if (!closed) slice = newSlice(slices, slice, area, dist);

                } else if (bk >= k1) slice.push(intersect(a, b, k1)); // ---|-->  |

            } else if (ak > k2) {

                if ((bk < k1)) { // <--|-----|---
                    slice.push(intersect(a, b, k2), intersect(a, b, k1));
                    if (!closed) slice = newSlice(slices, slice, area, dist);

                } else if (bk <= k2) slice.push(intersect(a, b, k2)); // |  <--|---

            } else {

                slice.push(a);

                if (bk < k1) { // <--|---  |
                    slice.push(intersect(a, b, k1));
                    if (!closed) slice = newSlice(slices, slice, area, dist);

                } else if (bk > k2) { // |  ---|-->
                    slice.push(intersect(a, b, k2));
                    if (!closed) slice = newSlice(slices, slice, area, dist);
                }
                // | --> |
            }
        }

        // add the last point
        a = points[len - 1];
        ak = a[axis];
        if (ak >= k1 && ak <= k2) slice.push(a);

        // close the polygon if its endpoints are not the same after clipping
        if (closed && slice[0] !== slice[slice.length - 1]) slice.push(slice[0]);

        // add the final slice
        newSlice(slices, slice, area, dist);
    }

    return slices;
}

function newSlice(slices, slice, area, dist) {
    if (slice.length) {
        // we don't recalculate the area/length of the unclipped geometry because the case where it goes
        // below the visibility threshold as a result of clipping is rare, so we avoid doing unnecessary work
        slice.area = area;
        slice.dist = dist;

        slices.push(slice);
    }
    return [];
}

},{}],2:[function(require,module,exports){
'use strict';

module.exports = convert;

var simplify = require('./simplify');

// converts GeoJSON feature into an intermediate projected JSON vector format with simplification data

function convert(data, tolerance) {
    var features = [];

    if (data.type === 'FeatureCollection') {
        for (var i = 0; i < data.features.length; i++) {
            convertFeature(features, data.features[i], tolerance);
        }
    } else if (data.type === 'Feature') {
        convertFeature(features, data, tolerance);

    } else {
        // single geometry or a geometry collection
        convertFeature(features, {geometry: data}, tolerance);
    }
    return features;
}

function convertFeature(features, feature, tolerance) {
    var geom = feature.geometry,
        type = geom.type,
        coords = geom.coordinates,
        tags = feature.properties,
        i, j, rings;

    if (type === 'Point') {
        features.push(create(tags, 1, [projectPoint(coords)]));

    } else if (type === 'MultiPoint') {
        features.push(create(tags, 1, project(coords)));

    } else if (type === 'LineString') {
        features.push(create(tags, 2, [project(coords, tolerance)]));

    } else if (type === 'MultiLineString' || type === 'Polygon') {
        rings = [];
        for (i = 0; i < coords.length; i++) {
            rings.push(project(coords[i], tolerance));
        }
        features.push(create(tags, type === 'Polygon' ? 3 : 2, rings));

    } else if (type === 'MultiPolygon') {
        rings = [];
        for (i = 0; i < coords.length; i++) {
            for (j = 0; j < coords[i].length; j++) {
                rings.push(project(coords[i][j], tolerance));
            }
        }
        features.push(create(tags, 3, rings));

    } else if (type === 'GeometryCollection') {
        for (i = 0; i < geom.geometries.length; i++) {
            convertFeature(features, {
                geometry: geom.geometries[i],
                properties: tags
            }, tolerance);
        }

    } else {
        console.warn('Unsupported GeoJSON type: ' + geom.type);
    }
}

function create(tags, type, geometry) {
    var feature = {
        geometry: geometry,
        type: type,
        tags: tags || null,
        min: [1, 1], // initial bbox values;
        max: [0, 0]  // note that all coords are in [0..1] range
    };
    calcBBox(feature);
    return feature;
}

function project(lonlats, tolerance) {
    var projected = [];
    for (var i = 0; i < lonlats.length; i++) {
        projected.push(projectPoint(lonlats[i]));
    }
    if (tolerance) {
        simplify(projected, tolerance);
        calcSize(projected);
    }
    return projected;
}

function projectPoint(p) {
    var sin = Math.sin(p[1] * Math.PI / 180),
        x = (p[0] / 360 + 0.5),
        y = (0.5 - 0.25 * Math.log((1 + sin) / (1 - sin)) / Math.PI);
    return [x, y, 0];
}

// calculate area and length of the poly
function calcSize(points) {
    var area = 0,
        dist = 0;

    for (var i = 0, a, b; i < points.length - 1; i++) {
        a = b || points[i];
        b = points[i + 1];

        area += a[0] * b[1] - b[0] * a[1];

        // use Manhattan distance instead of Euclidian one to avoid expensive square root computation
        dist += Math.abs(b[0] - a[0]) + Math.abs(b[1] - a[1]);
    }
    points.area = Math.abs(area / 2);
    points.dist = dist;
}

// calculate the feature bounding box for faster clipping later
function calcBBox(feature) {
    var geometry = feature.geometry,
        min = feature.min,
        max = feature.max;

    if (feature.type === 1) {
        calcRingBBox(min, max, geometry);
    } else {
        for (var i = 0; i < geometry.length; i++) {
            calcRingBBox(min, max, geometry[i]);
        }
    }
    return feature;
}

function calcRingBBox(min, max, points) {
    for (var i = 0, p; i < points.length; i++) {
        p = points[i];
        min[0] = Math.min(p[0], min[0]);
        max[0] = Math.max(p[0], max[0]);
        min[1] = Math.min(p[1], min[1]);
        max[1] = Math.max(p[1], max[1]);
    }
}

},{"./simplify":4}],3:[function(require,module,exports){
'use strict';

module.exports = geojsonvt;

var convert = require('./convert'), // GeoJSON conversion and preprocessing
    clip = require('./clip'),       // stripe clipping algorithm
    createTile = require('./tile'); // final simplified tile generation


function geojsonvt(data, options) {
    return new GeoJSONVT(data, options);
}

function GeoJSONVT(data, options) {
    options = this.options = extend(Object.create(this.options), options);

    var debug = options.debug;

    if (debug) console.time('preprocess data');

    var z2 = 1 << options.baseZoom, // 2^z
        features = convert(data, options.tolerance / (z2 * options.extent));

    this.tiles = {};

    if (debug) {
        console.timeEnd('preprocess data');
        console.time('generate tiles up to z' + options.maxZoom);
        this.stats = [];
        this.total = 0;
    }

    // start slicing from the top tile down
    this.splitTile(features, 0, 0, 0);

    if (debug) {
        console.log('features: %d, points: %d', this.tiles[0].numFeatures, this.tiles[0].numPoints);
        console.timeEnd('generate tiles up to z' + options.maxZoom);
        console.log('tiles generated:', this.total, this.stats);
    }
}

GeoJSONVT.prototype.options = {
    baseZoom: 14,   // max zoom to preserve detail on
    maxZoom: 4,     // zoom to slice down to on first pass
    maxPoints: 100, // stop slicing a tile below this number of points
    tolerance: 3,   // simplification tolerance (higher means simpler)
    extent: 4096,   // tile extent
    buffer: 64,     // tile buffer on each side
    debug: 0        // logging level (0, 1 or 2)
};

GeoJSONVT.prototype.splitTile = function (features, z, x, y, cz, cx, cy) {

    var stack = [features, z, x, y],
        options = this.options,
        debug = options.debug,
        extent = options.extent,
        buffer = options.buffer;

    // avoid recursion by using a processing queue
    while (stack.length) {
        features = stack.shift();
        z = stack.shift();
        x = stack.shift();
        y = stack.shift();

        var z2 = 1 << z,
            id = toID(z, x, y),
            tile = this.tiles[id],
            tileTolerance = z === options.baseZoom ? 0 : options.tolerance / (z2 * extent);

        if (!tile) {
            if (debug > 1) console.time('creation');

            tile = this.tiles[id] = createTile(features, z2, x, y, tileTolerance, extent, z === options.baseZoom);

            if (debug) {
                if (debug > 1) {
                    console.log('tile z%d-%d-%d (features: %d, points: %d, simplified: %d)',
                        z, x, y, tile.numFeatures, tile.numPoints, tile.numSimplified);
                    console.timeEnd('creation');
                }
                this.stats[z] = (this.stats[z] || 0) + 1;
                this.total++;
            }
        }

        if (!cz && (z === options.maxZoom || tile.numPoints <= options.maxPoints ||
                isClippedSquare(tile.features, extent, buffer)) || z === options.baseZoom || z === cz) {
            tile.source = features;
            continue; // stop tiling
        }

        if (cz) tile.source = features;
        else tile.source = null;

        if (debug > 1) console.time('clipping');

        // values we'll use for clipping
        var k1 = 0.5 * buffer / extent,
            k2 = 0.5 - k1,
            k3 = 0.5 + k1,
            k4 = 1 + k1,

            tl, bl, tr, br, left, right,
            m, goLeft, goTop;

        if (cz) { // if we have a specific tile to drill down to, calculate where to go
            m = 1 << (cz - z);
            goLeft = cx / m - x < 0.5;
            goTop = cy / m - y < 0.5;
        }

        tl = bl = tr = br = left = right = null;

        if (!cz ||  goLeft) left  = clip(features, z2, x - k1, x + k3, 0, intersectX);
        if (!cz || !goLeft) right = clip(features, z2, x + k2, x + k4, 0, intersectX);

        if (left) {
            if (!cz ||  goTop) tl = clip(left, z2, y - k1, y + k3, 1, intersectY);
            if (!cz || !goTop) bl = clip(left, z2, y + k2, y + k4, 1, intersectY);
        }

        if (right) {
            if (!cz ||  goTop) tr = clip(right, z2, y - k1, y + k3, 1, intersectY);
            if (!cz || !goTop) br = clip(right, z2, y + k2, y + k4, 1, intersectY);
        }

        if (debug > 1) console.timeEnd('clipping');

        if (tl) stack.push(tl, z + 1, x * 2,     y * 2);
        if (bl) stack.push(bl, z + 1, x * 2,     y * 2 + 1);
        if (tr) stack.push(tr, z + 1, x * 2 + 1, y * 2);
        if (br) stack.push(br, z + 1, x * 2 + 1, y * 2 + 1);
    }
};

GeoJSONVT.prototype.getTile = function (z, x, y) {
    var id = toID(z, x, y);
    if (this.tiles[id]) return this.tiles[id];

    var options = this.options,
        debug = options.debug;

    if (debug > 1) console.log('drilling down to z%d-%d-%d', z, x, y);

    var z0 = z,
        x0 = x,
        y0 = y,
        parent;

    while (!parent && z0 > 0) {
        z0--;
        x0 = Math.floor(x0 / 2);
        y0 = Math.floor(y0 / 2);
        parent = this.tiles[toID(z0, x0, y0)];
    }

    if (debug > 1) console.log('found parent tile z%d-%d-%d', z0, x0, y0);

    // if we found a parent tile containing the original geometry, we can drill down from it
    if (parent.source) {
        if (isClippedSquare(parent.features, options.extent, options.buffer)) return parent;

        if (debug) console.time('drilling down');
        this.splitTile(parent.source, z0, x0, y0, z, x, y);
        if (debug) console.timeEnd('drilling down');
    }

    return this.tiles[id];
};

// checks whether a tile is a whole-area fill after clipping; if it is, there's no sense slicing it further
function isClippedSquare(features, extent, buffer) {
    if (features.length !== 1) return false;

    var feature = features[0];
    if (feature.type !== 3 || feature.geometry.length > 1) return false;

    for (var i = 0; i < feature.geometry[0].length; i++) {
        var p = feature.geometry[0][i];
        if ((p[0] !== -buffer && p[0] !== extent + buffer) ||
            (p[1] !== -buffer && p[1] !== extent + buffer)) return false;
    }
    return true;
}

function toID(z, x, y) {
    return (((1 << z) * y + x) * 32) + z;
}

function intersectX(a, b, x) {
    return [x, (x - a[0]) * (b[1] - a[1]) / (b[0] - a[0]) + a[1], 1];
}
function intersectY(a, b, y) {
    return [(y - a[1]) * (b[0] - a[0]) / (b[1] - a[1]) + a[0], y, 1];
}

function extend(dest, src) {
    for (var i in src) {
        dest[i] = src[i];
    }
    return dest;
}

},{"./clip":1,"./convert":2,"./tile":5}],4:[function(require,module,exports){
'use strict';

module.exports = simplify;

// calculate simplification data using optimized Douglas-Peucker algorithm

function simplify(points, tolerance) {

    var sqTolerance = tolerance * tolerance,
        len = points.length,
        first = 0,
        last = len - 1,
        stack = [],
        i, maxSqDist, sqDist, index;

    // always retain the endpoints (1 is the max value)
    points[first][2] = 1;
    points[last][2] = 1;

    // avoid recursion by using a stack
    while (last) {

        maxSqDist = 0;

        for (i = first + 1; i < last; i++) {
            sqDist = getSqSegDist(points[i], points[first], points[last]);

            if (sqDist > maxSqDist) {
                index = i;
                maxSqDist = sqDist;
            }
        }

        if (maxSqDist > sqTolerance) {
            points[index][2] = maxSqDist; // save the point importance in squared pixels as a z coordinate
            stack.push(first, index, index, last);
        }

        last = stack.pop();
        first = stack.pop();
    }
}

// square distance from a point to a segment
function getSqSegDist(p, a, b) {

    var x = a[0], y = a[1],
        bx = b[0], by = b[1],
        px = p[0], py = p[1],
        dx = bx - x,
        dy = by - y;

    if (dx !== 0 || dy !== 0) {

        var t = ((px - x) * dx + (py - y) * dy) / (dx * dx + dy * dy);

        if (t > 1) {
            x = bx;
            y = by;

        } else if (t > 0) {
            x += dx * t;
            y += dy * t;
        }
    }

    dx = px - x;
    dy = py - y;

    return dx * dx + dy * dy;
}

},{}],5:[function(require,module,exports){
'use strict';

module.exports = createTile;

function createTile(features, z2, tx, ty, tolerance, extent, noSimplify) {
    var tile = {
        features: [],
        numPoints: 0,
        numSimplified: 0,
        numFeatures: 0,
        source: null
    };
    for (var i = 0; i < features.length; i++) {
        tile.numFeatures++;
        addFeature(tile, features[i], z2, tx, ty, tolerance, extent, noSimplify);
    }
    return tile;
}

function addFeature(tile, feature, z2, tx, ty, tolerance, extent, noSimplify) {

    var geom = feature.geometry,
        type = feature.type,
        transformed = [],
        sqTolerance = tolerance * tolerance,
        i, j, ring, p;

    if (type === 1) {
        for (i = 0; i < geom.length; i++) {
            transformed.push(transformPoint(geom[i], z2, tx, ty, extent));
            tile.numPoints++;
            tile.numSimplified++;
        }

    } else {

        // simplify and transform projected coordinates for tile geometry
        for (i = 0; i < geom.length; i++) {
            ring = geom[i];

            // filter out tiny polylines & polygons
            if (!noSimplify && ((type === 2 && ring.dist < tolerance) ||
                                (type === 3 && ring.area < sqTolerance))) {
                tile.numPoints += ring.length;
                continue;
            }

            var transformedRing = [];

            for (j = 0; j < ring.length; j++) {
                p = ring[j];
                // keep points with importance > tolerance
                if (noSimplify || p[2] > sqTolerance) {
                    transformedRing.push(transformPoint(p, z2, tx, ty, extent));
                    tile.numSimplified++;
                }
                tile.numPoints++;
            }

            transformed.push(transformedRing);
        }
    }

    if (transformed.length) {
        tile.features.push({
            geometry: transformed,
            type: type,
            tags: feature.tags || null
        });
    }
}

function transformPoint(p, z2, tx, ty, extent) {
    var x = Math.round(extent * (p[0] * z2 - tx)),
        y = Math.round(extent * (p[1] * z2 - ty));
    return [x, y];
}

},{}]},{},[3])(3)
});
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIm5vZGVfbW9kdWxlcy9icm93c2VyaWZ5L25vZGVfbW9kdWxlcy9icm93c2VyLXBhY2svX3ByZWx1ZGUuanMiLCJzcmMvY2xpcC5qcyIsInNyYy9jb252ZXJ0LmpzIiwic3JjL2luZGV4LmpzIiwic3JjL3NpbXBsaWZ5LmpzIiwic3JjL3RpbGUuanMiXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6IkFBQUE7QUNBQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDbEpBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ2hKQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzdNQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDdkVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSIsImZpbGUiOiJnZW5lcmF0ZWQuanMiLCJzb3VyY2VSb290IjoiIiwic291cmNlc0NvbnRlbnQiOlsiKGZ1bmN0aW9uIGUodCxuLHIpe2Z1bmN0aW9uIHMobyx1KXtpZighbltvXSl7aWYoIXRbb10pe3ZhciBhPXR5cGVvZiByZXF1aXJlPT1cImZ1bmN0aW9uXCImJnJlcXVpcmU7aWYoIXUmJmEpcmV0dXJuIGEobywhMCk7aWYoaSlyZXR1cm4gaShvLCEwKTt2YXIgZj1uZXcgRXJyb3IoXCJDYW5ub3QgZmluZCBtb2R1bGUgJ1wiK28rXCInXCIpO3Rocm93IGYuY29kZT1cIk1PRFVMRV9OT1RfRk9VTkRcIixmfXZhciBsPW5bb109e2V4cG9ydHM6e319O3Rbb11bMF0uY2FsbChsLmV4cG9ydHMsZnVuY3Rpb24oZSl7dmFyIG49dFtvXVsxXVtlXTtyZXR1cm4gcyhuP246ZSl9LGwsbC5leHBvcnRzLGUsdCxuLHIpfXJldHVybiBuW29dLmV4cG9ydHN9dmFyIGk9dHlwZW9mIHJlcXVpcmU9PVwiZnVuY3Rpb25cIiYmcmVxdWlyZTtmb3IodmFyIG89MDtvPHIubGVuZ3RoO28rKylzKHJbb10pO3JldHVybiBzfSkiLCIndXNlIHN0cmljdCc7XG5cbm1vZHVsZS5leHBvcnRzID0gY2xpcDtcblxuLyogY2xpcCBmZWF0dXJlcyBiZXR3ZWVuIHR3byBheGlzLXBhcmFsbGVsIGxpbmVzOlxuICogICAgIHwgICAgICAgIHxcbiAqICBfX198X19fICAgICB8ICAgICAvXG4gKiAvICAgfCAgIFxcX19fX3xfX19fL1xuICogICAgIHwgICAgICAgIHxcbiAqL1xuXG5mdW5jdGlvbiBjbGlwKGZlYXR1cmVzLCBzY2FsZSwgazEsIGsyLCBheGlzLCBpbnRlcnNlY3QpIHtcblxuICAgIGsxIC89IHNjYWxlO1xuICAgIGsyIC89IHNjYWxlO1xuXG4gICAgdmFyIGNsaXBwZWQgPSBbXTtcblxuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgZmVhdHVyZXMubGVuZ3RoOyBpKyspIHtcblxuICAgICAgICB2YXIgZmVhdHVyZSA9IGZlYXR1cmVzW2ldLFxuICAgICAgICAgICAgZ2VvbWV0cnkgPSBmZWF0dXJlLmdlb21ldHJ5LFxuICAgICAgICAgICAgdHlwZSA9IGZlYXR1cmUudHlwZSxcbiAgICAgICAgICAgIG1pbiwgbWF4O1xuXG4gICAgICAgIGlmIChmZWF0dXJlLm1pbikge1xuICAgICAgICAgICAgbWluID0gZmVhdHVyZS5taW5bYXhpc107XG4gICAgICAgICAgICBtYXggPSBmZWF0dXJlLm1heFtheGlzXTtcblxuICAgICAgICAgICAgaWYgKG1pbiA+PSBrMSAmJiBtYXggPD0gazIpIHsgLy8gdHJpdmlhbCBhY2NlcHRcbiAgICAgICAgICAgICAgICBjbGlwcGVkLnB1c2goZmVhdHVyZSk7XG4gICAgICAgICAgICAgICAgY29udGludWU7XG4gICAgICAgICAgICB9IGVsc2UgaWYgKG1pbiA+IGsyIHx8IG1heCA8IGsxKSBjb250aW51ZTsgLy8gdHJpdmlhbCByZWplY3RcbiAgICAgICAgfVxuXG4gICAgICAgIHZhciBzbGljZXMgPSB0eXBlID09PSAxID9cbiAgICAgICAgICAgICAgICBjbGlwUG9pbnRzKGdlb21ldHJ5LCBrMSwgazIsIGF4aXMpIDpcbiAgICAgICAgICAgICAgICBjbGlwR2VvbWV0cnkoZ2VvbWV0cnksIGsxLCBrMiwgYXhpcywgaW50ZXJzZWN0LCB0eXBlID09PSAzKTtcblxuICAgICAgICBpZiAoc2xpY2VzLmxlbmd0aCkge1xuICAgICAgICAgICAgLy8gaWYgYSBmZWF0dXJlIGdvdCBjbGlwcGVkLCBpdCB3aWxsIGxpa2VseSBnZXQgY2xpcHBlZCBvbiB0aGUgbmV4dCB6b29tIGxldmVsIGFzIHdlbGwsXG4gICAgICAgICAgICAvLyBzbyB0aGVyZSdzIG5vIG5lZWQgdG8gcmVjYWxjdWxhdGUgYmJveGVzXG4gICAgICAgICAgICBjbGlwcGVkLnB1c2goe1xuICAgICAgICAgICAgICAgIGdlb21ldHJ5OiBzbGljZXMsXG4gICAgICAgICAgICAgICAgdHlwZTogdHlwZSxcbiAgICAgICAgICAgICAgICB0YWdzOiBmZWF0dXJlc1tpXS50YWdzIHx8IG51bGxcbiAgICAgICAgICAgIH0pO1xuICAgICAgICB9XG4gICAgfVxuXG4gICAgcmV0dXJuIGNsaXBwZWQubGVuZ3RoID8gY2xpcHBlZCA6IG51bGw7XG59XG5cbmZ1bmN0aW9uIGNsaXBQb2ludHMoZ2VvbWV0cnksIGsxLCBrMiwgYXhpcykge1xuICAgIHZhciBzbGljZSA9IFtdO1xuXG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBnZW9tZXRyeS5sZW5ndGg7IGkrKykge1xuICAgICAgICB2YXIgYSA9IGdlb21ldHJ5W2ldLFxuICAgICAgICAgICAgYWsgPSBhW2F4aXNdO1xuXG4gICAgICAgIGlmIChhayA+PSBrMSAmJiBhayA8PSBrMikgc2xpY2UucHVzaChhKTtcbiAgICB9XG4gICAgcmV0dXJuIHNsaWNlO1xufVxuXG5mdW5jdGlvbiBjbGlwR2VvbWV0cnkoZ2VvbWV0cnksIGsxLCBrMiwgYXhpcywgaW50ZXJzZWN0LCBjbG9zZWQpIHtcblxuICAgIHZhciBzbGljZXMgPSBbXTtcblxuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgZ2VvbWV0cnkubGVuZ3RoOyBpKyspIHtcblxuICAgICAgICB2YXIgYWsgPSAwLFxuICAgICAgICAgICAgYmsgPSAwLFxuICAgICAgICAgICAgYiA9IG51bGwsXG4gICAgICAgICAgICBwb2ludHMgPSBnZW9tZXRyeVtpXSxcbiAgICAgICAgICAgIGFyZWEgPSBwb2ludHMuYXJlYSxcbiAgICAgICAgICAgIGRpc3QgPSBwb2ludHMuZGlzdCxcbiAgICAgICAgICAgIGxlbiA9IHBvaW50cy5sZW5ndGgsXG4gICAgICAgICAgICBhLCBqO1xuXG4gICAgICAgIHZhciBzbGljZSA9IFtdO1xuXG4gICAgICAgIGZvciAoaiA9IDA7IGogPCBsZW4gLSAxOyBqKyspIHtcbiAgICAgICAgICAgIGEgPSBiIHx8IHBvaW50c1tqXTtcbiAgICAgICAgICAgIGIgPSBwb2ludHNbaiArIDFdO1xuICAgICAgICAgICAgYWsgPSBiayB8fCBhW2F4aXNdO1xuICAgICAgICAgICAgYmsgPSBiW2F4aXNdO1xuXG4gICAgICAgICAgICBpZiAoYWsgPCBrMSkge1xuXG4gICAgICAgICAgICAgICAgaWYgKChiayA+IGsyKSkgeyAvLyAtLS18LS0tLS18LS0+XG4gICAgICAgICAgICAgICAgICAgIHNsaWNlLnB1c2goaW50ZXJzZWN0KGEsIGIsIGsxKSwgaW50ZXJzZWN0KGEsIGIsIGsyKSk7XG4gICAgICAgICAgICAgICAgICAgIGlmICghY2xvc2VkKSBzbGljZSA9IG5ld1NsaWNlKHNsaWNlcywgc2xpY2UsIGFyZWEsIGRpc3QpO1xuXG4gICAgICAgICAgICAgICAgfSBlbHNlIGlmIChiayA+PSBrMSkgc2xpY2UucHVzaChpbnRlcnNlY3QoYSwgYiwgazEpKTsgLy8gLS0tfC0tPiAgfFxuXG4gICAgICAgICAgICB9IGVsc2UgaWYgKGFrID4gazIpIHtcblxuICAgICAgICAgICAgICAgIGlmICgoYmsgPCBrMSkpIHsgLy8gPC0tfC0tLS0tfC0tLVxuICAgICAgICAgICAgICAgICAgICBzbGljZS5wdXNoKGludGVyc2VjdChhLCBiLCBrMiksIGludGVyc2VjdChhLCBiLCBrMSkpO1xuICAgICAgICAgICAgICAgICAgICBpZiAoIWNsb3NlZCkgc2xpY2UgPSBuZXdTbGljZShzbGljZXMsIHNsaWNlLCBhcmVhLCBkaXN0KTtcblxuICAgICAgICAgICAgICAgIH0gZWxzZSBpZiAoYmsgPD0gazIpIHNsaWNlLnB1c2goaW50ZXJzZWN0KGEsIGIsIGsyKSk7IC8vIHwgIDwtLXwtLS1cblxuICAgICAgICAgICAgfSBlbHNlIHtcblxuICAgICAgICAgICAgICAgIHNsaWNlLnB1c2goYSk7XG5cbiAgICAgICAgICAgICAgICBpZiAoYmsgPCBrMSkgeyAvLyA8LS18LS0tICB8XG4gICAgICAgICAgICAgICAgICAgIHNsaWNlLnB1c2goaW50ZXJzZWN0KGEsIGIsIGsxKSk7XG4gICAgICAgICAgICAgICAgICAgIGlmICghY2xvc2VkKSBzbGljZSA9IG5ld1NsaWNlKHNsaWNlcywgc2xpY2UsIGFyZWEsIGRpc3QpO1xuXG4gICAgICAgICAgICAgICAgfSBlbHNlIGlmIChiayA+IGsyKSB7IC8vIHwgIC0tLXwtLT5cbiAgICAgICAgICAgICAgICAgICAgc2xpY2UucHVzaChpbnRlcnNlY3QoYSwgYiwgazIpKTtcbiAgICAgICAgICAgICAgICAgICAgaWYgKCFjbG9zZWQpIHNsaWNlID0gbmV3U2xpY2Uoc2xpY2VzLCBzbGljZSwgYXJlYSwgZGlzdCk7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIC8vIHwgLS0+IHxcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuXG4gICAgICAgIC8vIGFkZCB0aGUgbGFzdCBwb2ludFxuICAgICAgICBhID0gcG9pbnRzW2xlbiAtIDFdO1xuICAgICAgICBhayA9IGFbYXhpc107XG4gICAgICAgIGlmIChhayA+PSBrMSAmJiBhayA8PSBrMikgc2xpY2UucHVzaChhKTtcblxuICAgICAgICAvLyBjbG9zZSB0aGUgcG9seWdvbiBpZiBpdHMgZW5kcG9pbnRzIGFyZSBub3QgdGhlIHNhbWUgYWZ0ZXIgY2xpcHBpbmdcbiAgICAgICAgaWYgKGNsb3NlZCAmJiBzbGljZVswXSAhPT0gc2xpY2Vbc2xpY2UubGVuZ3RoIC0gMV0pIHNsaWNlLnB1c2goc2xpY2VbMF0pO1xuXG4gICAgICAgIC8vIGFkZCB0aGUgZmluYWwgc2xpY2VcbiAgICAgICAgbmV3U2xpY2Uoc2xpY2VzLCBzbGljZSwgYXJlYSwgZGlzdCk7XG4gICAgfVxuXG4gICAgcmV0dXJuIHNsaWNlcztcbn1cblxuZnVuY3Rpb24gbmV3U2xpY2Uoc2xpY2VzLCBzbGljZSwgYXJlYSwgZGlzdCkge1xuICAgIGlmIChzbGljZS5sZW5ndGgpIHtcbiAgICAgICAgLy8gd2UgZG9uJ3QgcmVjYWxjdWxhdGUgdGhlIGFyZWEvbGVuZ3RoIG9mIHRoZSB1bmNsaXBwZWQgZ2VvbWV0cnkgYmVjYXVzZSB0aGUgY2FzZSB3aGVyZSBpdCBnb2VzXG4gICAgICAgIC8vIGJlbG93IHRoZSB2aXNpYmlsaXR5IHRocmVzaG9sZCBhcyBhIHJlc3VsdCBvZiBjbGlwcGluZyBpcyByYXJlLCBzbyB3ZSBhdm9pZCBkb2luZyB1bm5lY2Vzc2FyeSB3b3JrXG4gICAgICAgIHNsaWNlLmFyZWEgPSBhcmVhO1xuICAgICAgICBzbGljZS5kaXN0ID0gZGlzdDtcblxuICAgICAgICBzbGljZXMucHVzaChzbGljZSk7XG4gICAgfVxuICAgIHJldHVybiBbXTtcbn1cbiIsIid1c2Ugc3RyaWN0JztcblxubW9kdWxlLmV4cG9ydHMgPSBjb252ZXJ0O1xuXG52YXIgc2ltcGxpZnkgPSByZXF1aXJlKCcuL3NpbXBsaWZ5Jyk7XG5cbi8vIGNvbnZlcnRzIEdlb0pTT04gZmVhdHVyZSBpbnRvIGFuIGludGVybWVkaWF0ZSBwcm9qZWN0ZWQgSlNPTiB2ZWN0b3IgZm9ybWF0IHdpdGggc2ltcGxpZmljYXRpb24gZGF0YVxuXG5mdW5jdGlvbiBjb252ZXJ0KGRhdGEsIHRvbGVyYW5jZSkge1xuICAgIHZhciBmZWF0dXJlcyA9IFtdO1xuXG4gICAgaWYgKGRhdGEudHlwZSA9PT0gJ0ZlYXR1cmVDb2xsZWN0aW9uJykge1xuICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IGRhdGEuZmVhdHVyZXMubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgICAgIGNvbnZlcnRGZWF0dXJlKGZlYXR1cmVzLCBkYXRhLmZlYXR1cmVzW2ldLCB0b2xlcmFuY2UpO1xuICAgICAgICB9XG4gICAgfSBlbHNlIGlmIChkYXRhLnR5cGUgPT09ICdGZWF0dXJlJykge1xuICAgICAgICBjb252ZXJ0RmVhdHVyZShmZWF0dXJlcywgZGF0YSwgdG9sZXJhbmNlKTtcblxuICAgIH0gZWxzZSB7XG4gICAgICAgIC8vIHNpbmdsZSBnZW9tZXRyeSBvciBhIGdlb21ldHJ5IGNvbGxlY3Rpb25cbiAgICAgICAgY29udmVydEZlYXR1cmUoZmVhdHVyZXMsIHtnZW9tZXRyeTogZGF0YX0sIHRvbGVyYW5jZSk7XG4gICAgfVxuICAgIHJldHVybiBmZWF0dXJlcztcbn1cblxuZnVuY3Rpb24gY29udmVydEZlYXR1cmUoZmVhdHVyZXMsIGZlYXR1cmUsIHRvbGVyYW5jZSkge1xuICAgIHZhciBnZW9tID0gZmVhdHVyZS5nZW9tZXRyeSxcbiAgICAgICAgdHlwZSA9IGdlb20udHlwZSxcbiAgICAgICAgY29vcmRzID0gZ2VvbS5jb29yZGluYXRlcyxcbiAgICAgICAgdGFncyA9IGZlYXR1cmUucHJvcGVydGllcyxcbiAgICAgICAgaSwgaiwgcmluZ3M7XG5cbiAgICBpZiAodHlwZSA9PT0gJ1BvaW50Jykge1xuICAgICAgICBmZWF0dXJlcy5wdXNoKGNyZWF0ZSh0YWdzLCAxLCBbcHJvamVjdFBvaW50KGNvb3JkcyldKSk7XG5cbiAgICB9IGVsc2UgaWYgKHR5cGUgPT09ICdNdWx0aVBvaW50Jykge1xuICAgICAgICBmZWF0dXJlcy5wdXNoKGNyZWF0ZSh0YWdzLCAxLCBwcm9qZWN0KGNvb3JkcykpKTtcblxuICAgIH0gZWxzZSBpZiAodHlwZSA9PT0gJ0xpbmVTdHJpbmcnKSB7XG4gICAgICAgIGZlYXR1cmVzLnB1c2goY3JlYXRlKHRhZ3MsIDIsIFtwcm9qZWN0KGNvb3JkcywgdG9sZXJhbmNlKV0pKTtcblxuICAgIH0gZWxzZSBpZiAodHlwZSA9PT0gJ011bHRpTGluZVN0cmluZycgfHwgdHlwZSA9PT0gJ1BvbHlnb24nKSB7XG4gICAgICAgIHJpbmdzID0gW107XG4gICAgICAgIGZvciAoaSA9IDA7IGkgPCBjb29yZHMubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgICAgIHJpbmdzLnB1c2gocHJvamVjdChjb29yZHNbaV0sIHRvbGVyYW5jZSkpO1xuICAgICAgICB9XG4gICAgICAgIGZlYXR1cmVzLnB1c2goY3JlYXRlKHRhZ3MsIHR5cGUgPT09ICdQb2x5Z29uJyA/IDMgOiAyLCByaW5ncykpO1xuXG4gICAgfSBlbHNlIGlmICh0eXBlID09PSAnTXVsdGlQb2x5Z29uJykge1xuICAgICAgICByaW5ncyA9IFtdO1xuICAgICAgICBmb3IgKGkgPSAwOyBpIDwgY29vcmRzLmxlbmd0aDsgaSsrKSB7XG4gICAgICAgICAgICBmb3IgKGogPSAwOyBqIDwgY29vcmRzW2ldLmxlbmd0aDsgaisrKSB7XG4gICAgICAgICAgICAgICAgcmluZ3MucHVzaChwcm9qZWN0KGNvb3Jkc1tpXVtqXSwgdG9sZXJhbmNlKSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgZmVhdHVyZXMucHVzaChjcmVhdGUodGFncywgMywgcmluZ3MpKTtcblxuICAgIH0gZWxzZSBpZiAodHlwZSA9PT0gJ0dlb21ldHJ5Q29sbGVjdGlvbicpIHtcbiAgICAgICAgZm9yIChpID0gMDsgaSA8IGdlb20uZ2VvbWV0cmllcy5sZW5ndGg7IGkrKykge1xuICAgICAgICAgICAgY29udmVydEZlYXR1cmUoZmVhdHVyZXMsIHtcbiAgICAgICAgICAgICAgICBnZW9tZXRyeTogZ2VvbS5nZW9tZXRyaWVzW2ldLFxuICAgICAgICAgICAgICAgIHByb3BlcnRpZXM6IHRhZ3NcbiAgICAgICAgICAgIH0sIHRvbGVyYW5jZSk7XG4gICAgICAgIH1cblxuICAgIH0gZWxzZSB7XG4gICAgICAgIGNvbnNvbGUud2FybignVW5zdXBwb3J0ZWQgR2VvSlNPTiB0eXBlOiAnICsgZ2VvbS50eXBlKTtcbiAgICB9XG59XG5cbmZ1bmN0aW9uIGNyZWF0ZSh0YWdzLCB0eXBlLCBnZW9tZXRyeSkge1xuICAgIHZhciBmZWF0dXJlID0ge1xuICAgICAgICBnZW9tZXRyeTogZ2VvbWV0cnksXG4gICAgICAgIHR5cGU6IHR5cGUsXG4gICAgICAgIHRhZ3M6IHRhZ3MgfHwgbnVsbCxcbiAgICAgICAgbWluOiBbMSwgMV0sIC8vIGluaXRpYWwgYmJveCB2YWx1ZXM7XG4gICAgICAgIG1heDogWzAsIDBdICAvLyBub3RlIHRoYXQgYWxsIGNvb3JkcyBhcmUgaW4gWzAuLjFdIHJhbmdlXG4gICAgfTtcbiAgICBjYWxjQkJveChmZWF0dXJlKTtcbiAgICByZXR1cm4gZmVhdHVyZTtcbn1cblxuZnVuY3Rpb24gcHJvamVjdChsb25sYXRzLCB0b2xlcmFuY2UpIHtcbiAgICB2YXIgcHJvamVjdGVkID0gW107XG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBsb25sYXRzLmxlbmd0aDsgaSsrKSB7XG4gICAgICAgIHByb2plY3RlZC5wdXNoKHByb2plY3RQb2ludChsb25sYXRzW2ldKSk7XG4gICAgfVxuICAgIGlmICh0b2xlcmFuY2UpIHtcbiAgICAgICAgc2ltcGxpZnkocHJvamVjdGVkLCB0b2xlcmFuY2UpO1xuICAgICAgICBjYWxjU2l6ZShwcm9qZWN0ZWQpO1xuICAgIH1cbiAgICByZXR1cm4gcHJvamVjdGVkO1xufVxuXG5mdW5jdGlvbiBwcm9qZWN0UG9pbnQocCkge1xuICAgIHZhciBzaW4gPSBNYXRoLnNpbihwWzFdICogTWF0aC5QSSAvIDE4MCksXG4gICAgICAgIHggPSAocFswXSAvIDM2MCArIDAuNSksXG4gICAgICAgIHkgPSAoMC41IC0gMC4yNSAqIE1hdGgubG9nKCgxICsgc2luKSAvICgxIC0gc2luKSkgLyBNYXRoLlBJKTtcbiAgICByZXR1cm4gW3gsIHksIDBdO1xufVxuXG4vLyBjYWxjdWxhdGUgYXJlYSBhbmQgbGVuZ3RoIG9mIHRoZSBwb2x5XG5mdW5jdGlvbiBjYWxjU2l6ZShwb2ludHMpIHtcbiAgICB2YXIgYXJlYSA9IDAsXG4gICAgICAgIGRpc3QgPSAwO1xuXG4gICAgZm9yICh2YXIgaSA9IDAsIGEsIGI7IGkgPCBwb2ludHMubGVuZ3RoIC0gMTsgaSsrKSB7XG4gICAgICAgIGEgPSBiIHx8IHBvaW50c1tpXTtcbiAgICAgICAgYiA9IHBvaW50c1tpICsgMV07XG5cbiAgICAgICAgYXJlYSArPSBhWzBdICogYlsxXSAtIGJbMF0gKiBhWzFdO1xuXG4gICAgICAgIC8vIHVzZSBNYW5oYXR0YW4gZGlzdGFuY2UgaW5zdGVhZCBvZiBFdWNsaWRpYW4gb25lIHRvIGF2b2lkIGV4cGVuc2l2ZSBzcXVhcmUgcm9vdCBjb21wdXRhdGlvblxuICAgICAgICBkaXN0ICs9IE1hdGguYWJzKGJbMF0gLSBhWzBdKSArIE1hdGguYWJzKGJbMV0gLSBhWzFdKTtcbiAgICB9XG4gICAgcG9pbnRzLmFyZWEgPSBNYXRoLmFicyhhcmVhIC8gMik7XG4gICAgcG9pbnRzLmRpc3QgPSBkaXN0O1xufVxuXG4vLyBjYWxjdWxhdGUgdGhlIGZlYXR1cmUgYm91bmRpbmcgYm94IGZvciBmYXN0ZXIgY2xpcHBpbmcgbGF0ZXJcbmZ1bmN0aW9uIGNhbGNCQm94KGZlYXR1cmUpIHtcbiAgICB2YXIgZ2VvbWV0cnkgPSBmZWF0dXJlLmdlb21ldHJ5LFxuICAgICAgICBtaW4gPSBmZWF0dXJlLm1pbixcbiAgICAgICAgbWF4ID0gZmVhdHVyZS5tYXg7XG5cbiAgICBpZiAoZmVhdHVyZS50eXBlID09PSAxKSB7XG4gICAgICAgIGNhbGNSaW5nQkJveChtaW4sIG1heCwgZ2VvbWV0cnkpO1xuICAgIH0gZWxzZSB7XG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgZ2VvbWV0cnkubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgICAgIGNhbGNSaW5nQkJveChtaW4sIG1heCwgZ2VvbWV0cnlbaV0pO1xuICAgICAgICB9XG4gICAgfVxuICAgIHJldHVybiBmZWF0dXJlO1xufVxuXG5mdW5jdGlvbiBjYWxjUmluZ0JCb3gobWluLCBtYXgsIHBvaW50cykge1xuICAgIGZvciAodmFyIGkgPSAwLCBwOyBpIDwgcG9pbnRzLmxlbmd0aDsgaSsrKSB7XG4gICAgICAgIHAgPSBwb2ludHNbaV07XG4gICAgICAgIG1pblswXSA9IE1hdGgubWluKHBbMF0sIG1pblswXSk7XG4gICAgICAgIG1heFswXSA9IE1hdGgubWF4KHBbMF0sIG1heFswXSk7XG4gICAgICAgIG1pblsxXSA9IE1hdGgubWluKHBbMV0sIG1pblsxXSk7XG4gICAgICAgIG1heFsxXSA9IE1hdGgubWF4KHBbMV0sIG1heFsxXSk7XG4gICAgfVxufVxuIiwiJ3VzZSBzdHJpY3QnO1xuXG5tb2R1bGUuZXhwb3J0cyA9IGdlb2pzb252dDtcblxudmFyIGNvbnZlcnQgPSByZXF1aXJlKCcuL2NvbnZlcnQnKSwgLy8gR2VvSlNPTiBjb252ZXJzaW9uIGFuZCBwcmVwcm9jZXNzaW5nXG4gICAgY2xpcCA9IHJlcXVpcmUoJy4vY2xpcCcpLCAgICAgICAvLyBzdHJpcGUgY2xpcHBpbmcgYWxnb3JpdGhtXG4gICAgY3JlYXRlVGlsZSA9IHJlcXVpcmUoJy4vdGlsZScpOyAvLyBmaW5hbCBzaW1wbGlmaWVkIHRpbGUgZ2VuZXJhdGlvblxuXG5cbmZ1bmN0aW9uIGdlb2pzb252dChkYXRhLCBvcHRpb25zKSB7XG4gICAgcmV0dXJuIG5ldyBHZW9KU09OVlQoZGF0YSwgb3B0aW9ucyk7XG59XG5cbmZ1bmN0aW9uIEdlb0pTT05WVChkYXRhLCBvcHRpb25zKSB7XG4gICAgb3B0aW9ucyA9IHRoaXMub3B0aW9ucyA9IGV4dGVuZChPYmplY3QuY3JlYXRlKHRoaXMub3B0aW9ucyksIG9wdGlvbnMpO1xuXG4gICAgdmFyIGRlYnVnID0gb3B0aW9ucy5kZWJ1ZztcblxuICAgIGlmIChkZWJ1ZykgY29uc29sZS50aW1lKCdwcmVwcm9jZXNzIGRhdGEnKTtcblxuICAgIHZhciB6MiA9IDEgPDwgb3B0aW9ucy5iYXNlWm9vbSwgLy8gMl56XG4gICAgICAgIGZlYXR1cmVzID0gY29udmVydChkYXRhLCBvcHRpb25zLnRvbGVyYW5jZSAvICh6MiAqIG9wdGlvbnMuZXh0ZW50KSk7XG5cbiAgICB0aGlzLnRpbGVzID0ge307XG5cbiAgICBpZiAoZGVidWcpIHtcbiAgICAgICAgY29uc29sZS50aW1lRW5kKCdwcmVwcm9jZXNzIGRhdGEnKTtcbiAgICAgICAgY29uc29sZS50aW1lKCdnZW5lcmF0ZSB0aWxlcyB1cCB0byB6JyArIG9wdGlvbnMubWF4Wm9vbSk7XG4gICAgICAgIHRoaXMuc3RhdHMgPSBbXTtcbiAgICAgICAgdGhpcy50b3RhbCA9IDA7XG4gICAgfVxuXG4gICAgLy8gc3RhcnQgc2xpY2luZyBmcm9tIHRoZSB0b3AgdGlsZSBkb3duXG4gICAgdGhpcy5zcGxpdFRpbGUoZmVhdHVyZXMsIDAsIDAsIDApO1xuXG4gICAgaWYgKGRlYnVnKSB7XG4gICAgICAgIGNvbnNvbGUubG9nKCdmZWF0dXJlczogJWQsIHBvaW50czogJWQnLCB0aGlzLnRpbGVzWzBdLm51bUZlYXR1cmVzLCB0aGlzLnRpbGVzWzBdLm51bVBvaW50cyk7XG4gICAgICAgIGNvbnNvbGUudGltZUVuZCgnZ2VuZXJhdGUgdGlsZXMgdXAgdG8geicgKyBvcHRpb25zLm1heFpvb20pO1xuICAgICAgICBjb25zb2xlLmxvZygndGlsZXMgZ2VuZXJhdGVkOicsIHRoaXMudG90YWwsIHRoaXMuc3RhdHMpO1xuICAgIH1cbn1cblxuR2VvSlNPTlZULnByb3RvdHlwZS5vcHRpb25zID0ge1xuICAgIGJhc2Vab29tOiAxNCwgICAvLyBtYXggem9vbSB0byBwcmVzZXJ2ZSBkZXRhaWwgb25cbiAgICBtYXhab29tOiA0LCAgICAgLy8gem9vbSB0byBzbGljZSBkb3duIHRvIG9uIGZpcnN0IHBhc3NcbiAgICBtYXhQb2ludHM6IDEwMCwgLy8gc3RvcCBzbGljaW5nIGEgdGlsZSBiZWxvdyB0aGlzIG51bWJlciBvZiBwb2ludHNcbiAgICB0b2xlcmFuY2U6IDMsICAgLy8gc2ltcGxpZmljYXRpb24gdG9sZXJhbmNlIChoaWdoZXIgbWVhbnMgc2ltcGxlcilcbiAgICBleHRlbnQ6IDQwOTYsICAgLy8gdGlsZSBleHRlbnRcbiAgICBidWZmZXI6IDY0LCAgICAgLy8gdGlsZSBidWZmZXIgb24gZWFjaCBzaWRlXG4gICAgZGVidWc6IDAgICAgICAgIC8vIGxvZ2dpbmcgbGV2ZWwgKDAsIDEgb3IgMilcbn07XG5cbkdlb0pTT05WVC5wcm90b3R5cGUuc3BsaXRUaWxlID0gZnVuY3Rpb24gKGZlYXR1cmVzLCB6LCB4LCB5LCBjeiwgY3gsIGN5KSB7XG5cbiAgICB2YXIgc3RhY2sgPSBbZmVhdHVyZXMsIHosIHgsIHldLFxuICAgICAgICBvcHRpb25zID0gdGhpcy5vcHRpb25zLFxuICAgICAgICBkZWJ1ZyA9IG9wdGlvbnMuZGVidWcsXG4gICAgICAgIGV4dGVudCA9IG9wdGlvbnMuZXh0ZW50LFxuICAgICAgICBidWZmZXIgPSBvcHRpb25zLmJ1ZmZlcjtcblxuICAgIC8vIGF2b2lkIHJlY3Vyc2lvbiBieSB1c2luZyBhIHByb2Nlc3NpbmcgcXVldWVcbiAgICB3aGlsZSAoc3RhY2subGVuZ3RoKSB7XG4gICAgICAgIGZlYXR1cmVzID0gc3RhY2suc2hpZnQoKTtcbiAgICAgICAgeiA9IHN0YWNrLnNoaWZ0KCk7XG4gICAgICAgIHggPSBzdGFjay5zaGlmdCgpO1xuICAgICAgICB5ID0gc3RhY2suc2hpZnQoKTtcblxuICAgICAgICB2YXIgejIgPSAxIDw8IHosXG4gICAgICAgICAgICBpZCA9IHRvSUQoeiwgeCwgeSksXG4gICAgICAgICAgICB0aWxlID0gdGhpcy50aWxlc1tpZF0sXG4gICAgICAgICAgICB0aWxlVG9sZXJhbmNlID0geiA9PT0gb3B0aW9ucy5iYXNlWm9vbSA/IDAgOiBvcHRpb25zLnRvbGVyYW5jZSAvICh6MiAqIGV4dGVudCk7XG5cbiAgICAgICAgaWYgKCF0aWxlKSB7XG4gICAgICAgICAgICBpZiAoZGVidWcgPiAxKSBjb25zb2xlLnRpbWUoJ2NyZWF0aW9uJyk7XG5cbiAgICAgICAgICAgIHRpbGUgPSB0aGlzLnRpbGVzW2lkXSA9IGNyZWF0ZVRpbGUoZmVhdHVyZXMsIHoyLCB4LCB5LCB0aWxlVG9sZXJhbmNlLCBleHRlbnQsIHogPT09IG9wdGlvbnMuYmFzZVpvb20pO1xuXG4gICAgICAgICAgICBpZiAoZGVidWcpIHtcbiAgICAgICAgICAgICAgICBpZiAoZGVidWcgPiAxKSB7XG4gICAgICAgICAgICAgICAgICAgIGNvbnNvbGUubG9nKCd0aWxlIHolZC0lZC0lZCAoZmVhdHVyZXM6ICVkLCBwb2ludHM6ICVkLCBzaW1wbGlmaWVkOiAlZCknLFxuICAgICAgICAgICAgICAgICAgICAgICAgeiwgeCwgeSwgdGlsZS5udW1GZWF0dXJlcywgdGlsZS5udW1Qb2ludHMsIHRpbGUubnVtU2ltcGxpZmllZCk7XG4gICAgICAgICAgICAgICAgICAgIGNvbnNvbGUudGltZUVuZCgnY3JlYXRpb24nKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgdGhpcy5zdGF0c1t6XSA9ICh0aGlzLnN0YXRzW3pdIHx8IDApICsgMTtcbiAgICAgICAgICAgICAgICB0aGlzLnRvdGFsKys7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cblxuICAgICAgICBpZiAoIWN6ICYmICh6ID09PSBvcHRpb25zLm1heFpvb20gfHwgdGlsZS5udW1Qb2ludHMgPD0gb3B0aW9ucy5tYXhQb2ludHMgfHxcbiAgICAgICAgICAgICAgICBpc0NsaXBwZWRTcXVhcmUodGlsZS5mZWF0dXJlcywgZXh0ZW50LCBidWZmZXIpKSB8fCB6ID09PSBvcHRpb25zLmJhc2Vab29tIHx8IHogPT09IGN6KSB7XG4gICAgICAgICAgICB0aWxlLnNvdXJjZSA9IGZlYXR1cmVzO1xuICAgICAgICAgICAgY29udGludWU7IC8vIHN0b3AgdGlsaW5nXG4gICAgICAgIH1cblxuICAgICAgICBpZiAoY3opIHRpbGUuc291cmNlID0gZmVhdHVyZXM7XG4gICAgICAgIGVsc2UgdGlsZS5zb3VyY2UgPSBudWxsO1xuXG4gICAgICAgIGlmIChkZWJ1ZyA+IDEpIGNvbnNvbGUudGltZSgnY2xpcHBpbmcnKTtcblxuICAgICAgICAvLyB2YWx1ZXMgd2UnbGwgdXNlIGZvciBjbGlwcGluZ1xuICAgICAgICB2YXIgazEgPSAwLjUgKiBidWZmZXIgLyBleHRlbnQsXG4gICAgICAgICAgICBrMiA9IDAuNSAtIGsxLFxuICAgICAgICAgICAgazMgPSAwLjUgKyBrMSxcbiAgICAgICAgICAgIGs0ID0gMSArIGsxLFxuXG4gICAgICAgICAgICB0bCwgYmwsIHRyLCBiciwgbGVmdCwgcmlnaHQsXG4gICAgICAgICAgICBtLCBnb0xlZnQsIGdvVG9wO1xuXG4gICAgICAgIGlmIChjeikgeyAvLyBpZiB3ZSBoYXZlIGEgc3BlY2lmaWMgdGlsZSB0byBkcmlsbCBkb3duIHRvLCBjYWxjdWxhdGUgd2hlcmUgdG8gZ29cbiAgICAgICAgICAgIG0gPSAxIDw8IChjeiAtIHopO1xuICAgICAgICAgICAgZ29MZWZ0ID0gY3ggLyBtIC0geCA8IDAuNTtcbiAgICAgICAgICAgIGdvVG9wID0gY3kgLyBtIC0geSA8IDAuNTtcbiAgICAgICAgfVxuXG4gICAgICAgIHRsID0gYmwgPSB0ciA9IGJyID0gbGVmdCA9IHJpZ2h0ID0gbnVsbDtcblxuICAgICAgICBpZiAoIWN6IHx8ICBnb0xlZnQpIGxlZnQgID0gY2xpcChmZWF0dXJlcywgejIsIHggLSBrMSwgeCArIGszLCAwLCBpbnRlcnNlY3RYKTtcbiAgICAgICAgaWYgKCFjeiB8fCAhZ29MZWZ0KSByaWdodCA9IGNsaXAoZmVhdHVyZXMsIHoyLCB4ICsgazIsIHggKyBrNCwgMCwgaW50ZXJzZWN0WCk7XG5cbiAgICAgICAgaWYgKGxlZnQpIHtcbiAgICAgICAgICAgIGlmICghY3ogfHwgIGdvVG9wKSB0bCA9IGNsaXAobGVmdCwgejIsIHkgLSBrMSwgeSArIGszLCAxLCBpbnRlcnNlY3RZKTtcbiAgICAgICAgICAgIGlmICghY3ogfHwgIWdvVG9wKSBibCA9IGNsaXAobGVmdCwgejIsIHkgKyBrMiwgeSArIGs0LCAxLCBpbnRlcnNlY3RZKTtcbiAgICAgICAgfVxuXG4gICAgICAgIGlmIChyaWdodCkge1xuICAgICAgICAgICAgaWYgKCFjeiB8fCAgZ29Ub3ApIHRyID0gY2xpcChyaWdodCwgejIsIHkgLSBrMSwgeSArIGszLCAxLCBpbnRlcnNlY3RZKTtcbiAgICAgICAgICAgIGlmICghY3ogfHwgIWdvVG9wKSBiciA9IGNsaXAocmlnaHQsIHoyLCB5ICsgazIsIHkgKyBrNCwgMSwgaW50ZXJzZWN0WSk7XG4gICAgICAgIH1cblxuICAgICAgICBpZiAoZGVidWcgPiAxKSBjb25zb2xlLnRpbWVFbmQoJ2NsaXBwaW5nJyk7XG5cbiAgICAgICAgaWYgKHRsKSBzdGFjay5wdXNoKHRsLCB6ICsgMSwgeCAqIDIsICAgICB5ICogMik7XG4gICAgICAgIGlmIChibCkgc3RhY2sucHVzaChibCwgeiArIDEsIHggKiAyLCAgICAgeSAqIDIgKyAxKTtcbiAgICAgICAgaWYgKHRyKSBzdGFjay5wdXNoKHRyLCB6ICsgMSwgeCAqIDIgKyAxLCB5ICogMik7XG4gICAgICAgIGlmIChicikgc3RhY2sucHVzaChiciwgeiArIDEsIHggKiAyICsgMSwgeSAqIDIgKyAxKTtcbiAgICB9XG59O1xuXG5HZW9KU09OVlQucHJvdG90eXBlLmdldFRpbGUgPSBmdW5jdGlvbiAoeiwgeCwgeSkge1xuICAgIHZhciBpZCA9IHRvSUQoeiwgeCwgeSk7XG4gICAgaWYgKHRoaXMudGlsZXNbaWRdKSByZXR1cm4gdGhpcy50aWxlc1tpZF07XG5cbiAgICB2YXIgb3B0aW9ucyA9IHRoaXMub3B0aW9ucyxcbiAgICAgICAgZGVidWcgPSBvcHRpb25zLmRlYnVnO1xuXG4gICAgaWYgKGRlYnVnID4gMSkgY29uc29sZS5sb2coJ2RyaWxsaW5nIGRvd24gdG8geiVkLSVkLSVkJywgeiwgeCwgeSk7XG5cbiAgICB2YXIgejAgPSB6LFxuICAgICAgICB4MCA9IHgsXG4gICAgICAgIHkwID0geSxcbiAgICAgICAgcGFyZW50O1xuXG4gICAgd2hpbGUgKCFwYXJlbnQgJiYgejAgPiAwKSB7XG4gICAgICAgIHowLS07XG4gICAgICAgIHgwID0gTWF0aC5mbG9vcih4MCAvIDIpO1xuICAgICAgICB5MCA9IE1hdGguZmxvb3IoeTAgLyAyKTtcbiAgICAgICAgcGFyZW50ID0gdGhpcy50aWxlc1t0b0lEKHowLCB4MCwgeTApXTtcbiAgICB9XG5cbiAgICBpZiAoZGVidWcgPiAxKSBjb25zb2xlLmxvZygnZm91bmQgcGFyZW50IHRpbGUgeiVkLSVkLSVkJywgejAsIHgwLCB5MCk7XG5cbiAgICAvLyBpZiB3ZSBmb3VuZCBhIHBhcmVudCB0aWxlIGNvbnRhaW5pbmcgdGhlIG9yaWdpbmFsIGdlb21ldHJ5LCB3ZSBjYW4gZHJpbGwgZG93biBmcm9tIGl0XG4gICAgaWYgKHBhcmVudC5zb3VyY2UpIHtcbiAgICAgICAgaWYgKGlzQ2xpcHBlZFNxdWFyZShwYXJlbnQuZmVhdHVyZXMsIG9wdGlvbnMuZXh0ZW50LCBvcHRpb25zLmJ1ZmZlcikpIHJldHVybiBwYXJlbnQ7XG5cbiAgICAgICAgaWYgKGRlYnVnKSBjb25zb2xlLnRpbWUoJ2RyaWxsaW5nIGRvd24nKTtcbiAgICAgICAgdGhpcy5zcGxpdFRpbGUocGFyZW50LnNvdXJjZSwgejAsIHgwLCB5MCwgeiwgeCwgeSk7XG4gICAgICAgIGlmIChkZWJ1ZykgY29uc29sZS50aW1lRW5kKCdkcmlsbGluZyBkb3duJyk7XG4gICAgfVxuXG4gICAgcmV0dXJuIHRoaXMudGlsZXNbaWRdO1xufTtcblxuLy8gY2hlY2tzIHdoZXRoZXIgYSB0aWxlIGlzIGEgd2hvbGUtYXJlYSBmaWxsIGFmdGVyIGNsaXBwaW5nOyBpZiBpdCBpcywgdGhlcmUncyBubyBzZW5zZSBzbGljaW5nIGl0IGZ1cnRoZXJcbmZ1bmN0aW9uIGlzQ2xpcHBlZFNxdWFyZShmZWF0dXJlcywgZXh0ZW50LCBidWZmZXIpIHtcbiAgICBpZiAoZmVhdHVyZXMubGVuZ3RoICE9PSAxKSByZXR1cm4gZmFsc2U7XG5cbiAgICB2YXIgZmVhdHVyZSA9IGZlYXR1cmVzWzBdO1xuICAgIGlmIChmZWF0dXJlLnR5cGUgIT09IDMgfHwgZmVhdHVyZS5nZW9tZXRyeS5sZW5ndGggPiAxKSByZXR1cm4gZmFsc2U7XG5cbiAgICBmb3IgKHZhciBpID0gMDsgaSA8IGZlYXR1cmUuZ2VvbWV0cnlbMF0ubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgdmFyIHAgPSBmZWF0dXJlLmdlb21ldHJ5WzBdW2ldO1xuICAgICAgICBpZiAoKHBbMF0gIT09IC1idWZmZXIgJiYgcFswXSAhPT0gZXh0ZW50ICsgYnVmZmVyKSB8fFxuICAgICAgICAgICAgKHBbMV0gIT09IC1idWZmZXIgJiYgcFsxXSAhPT0gZXh0ZW50ICsgYnVmZmVyKSkgcmV0dXJuIGZhbHNlO1xuICAgIH1cbiAgICByZXR1cm4gdHJ1ZTtcbn1cblxuZnVuY3Rpb24gdG9JRCh6LCB4LCB5KSB7XG4gICAgcmV0dXJuICgoKDEgPDwgeikgKiB5ICsgeCkgKiAzMikgKyB6O1xufVxuXG5mdW5jdGlvbiBpbnRlcnNlY3RYKGEsIGIsIHgpIHtcbiAgICByZXR1cm4gW3gsICh4IC0gYVswXSkgKiAoYlsxXSAtIGFbMV0pIC8gKGJbMF0gLSBhWzBdKSArIGFbMV0sIDFdO1xufVxuZnVuY3Rpb24gaW50ZXJzZWN0WShhLCBiLCB5KSB7XG4gICAgcmV0dXJuIFsoeSAtIGFbMV0pICogKGJbMF0gLSBhWzBdKSAvIChiWzFdIC0gYVsxXSkgKyBhWzBdLCB5LCAxXTtcbn1cblxuZnVuY3Rpb24gZXh0ZW5kKGRlc3QsIHNyYykge1xuICAgIGZvciAodmFyIGkgaW4gc3JjKSB7XG4gICAgICAgIGRlc3RbaV0gPSBzcmNbaV07XG4gICAgfVxuICAgIHJldHVybiBkZXN0O1xufVxuIiwiJ3VzZSBzdHJpY3QnO1xuXG5tb2R1bGUuZXhwb3J0cyA9IHNpbXBsaWZ5O1xuXG4vLyBjYWxjdWxhdGUgc2ltcGxpZmljYXRpb24gZGF0YSB1c2luZyBvcHRpbWl6ZWQgRG91Z2xhcy1QZXVja2VyIGFsZ29yaXRobVxuXG5mdW5jdGlvbiBzaW1wbGlmeShwb2ludHMsIHRvbGVyYW5jZSkge1xuXG4gICAgdmFyIHNxVG9sZXJhbmNlID0gdG9sZXJhbmNlICogdG9sZXJhbmNlLFxuICAgICAgICBsZW4gPSBwb2ludHMubGVuZ3RoLFxuICAgICAgICBmaXJzdCA9IDAsXG4gICAgICAgIGxhc3QgPSBsZW4gLSAxLFxuICAgICAgICBzdGFjayA9IFtdLFxuICAgICAgICBpLCBtYXhTcURpc3QsIHNxRGlzdCwgaW5kZXg7XG5cbiAgICAvLyBhbHdheXMgcmV0YWluIHRoZSBlbmRwb2ludHMgKDEgaXMgdGhlIG1heCB2YWx1ZSlcbiAgICBwb2ludHNbZmlyc3RdWzJdID0gMTtcbiAgICBwb2ludHNbbGFzdF1bMl0gPSAxO1xuXG4gICAgLy8gYXZvaWQgcmVjdXJzaW9uIGJ5IHVzaW5nIGEgc3RhY2tcbiAgICB3aGlsZSAobGFzdCkge1xuXG4gICAgICAgIG1heFNxRGlzdCA9IDA7XG5cbiAgICAgICAgZm9yIChpID0gZmlyc3QgKyAxOyBpIDwgbGFzdDsgaSsrKSB7XG4gICAgICAgICAgICBzcURpc3QgPSBnZXRTcVNlZ0Rpc3QocG9pbnRzW2ldLCBwb2ludHNbZmlyc3RdLCBwb2ludHNbbGFzdF0pO1xuXG4gICAgICAgICAgICBpZiAoc3FEaXN0ID4gbWF4U3FEaXN0KSB7XG4gICAgICAgICAgICAgICAgaW5kZXggPSBpO1xuICAgICAgICAgICAgICAgIG1heFNxRGlzdCA9IHNxRGlzdDtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuXG4gICAgICAgIGlmIChtYXhTcURpc3QgPiBzcVRvbGVyYW5jZSkge1xuICAgICAgICAgICAgcG9pbnRzW2luZGV4XVsyXSA9IG1heFNxRGlzdDsgLy8gc2F2ZSB0aGUgcG9pbnQgaW1wb3J0YW5jZSBpbiBzcXVhcmVkIHBpeGVscyBhcyBhIHogY29vcmRpbmF0ZVxuICAgICAgICAgICAgc3RhY2sucHVzaChmaXJzdCwgaW5kZXgsIGluZGV4LCBsYXN0KTtcbiAgICAgICAgfVxuXG4gICAgICAgIGxhc3QgPSBzdGFjay5wb3AoKTtcbiAgICAgICAgZmlyc3QgPSBzdGFjay5wb3AoKTtcbiAgICB9XG59XG5cbi8vIHNxdWFyZSBkaXN0YW5jZSBmcm9tIGEgcG9pbnQgdG8gYSBzZWdtZW50XG5mdW5jdGlvbiBnZXRTcVNlZ0Rpc3QocCwgYSwgYikge1xuXG4gICAgdmFyIHggPSBhWzBdLCB5ID0gYVsxXSxcbiAgICAgICAgYnggPSBiWzBdLCBieSA9IGJbMV0sXG4gICAgICAgIHB4ID0gcFswXSwgcHkgPSBwWzFdLFxuICAgICAgICBkeCA9IGJ4IC0geCxcbiAgICAgICAgZHkgPSBieSAtIHk7XG5cbiAgICBpZiAoZHggIT09IDAgfHwgZHkgIT09IDApIHtcblxuICAgICAgICB2YXIgdCA9ICgocHggLSB4KSAqIGR4ICsgKHB5IC0geSkgKiBkeSkgLyAoZHggKiBkeCArIGR5ICogZHkpO1xuXG4gICAgICAgIGlmICh0ID4gMSkge1xuICAgICAgICAgICAgeCA9IGJ4O1xuICAgICAgICAgICAgeSA9IGJ5O1xuXG4gICAgICAgIH0gZWxzZSBpZiAodCA+IDApIHtcbiAgICAgICAgICAgIHggKz0gZHggKiB0O1xuICAgICAgICAgICAgeSArPSBkeSAqIHQ7XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICBkeCA9IHB4IC0geDtcbiAgICBkeSA9IHB5IC0geTtcblxuICAgIHJldHVybiBkeCAqIGR4ICsgZHkgKiBkeTtcbn1cbiIsIid1c2Ugc3RyaWN0JztcblxubW9kdWxlLmV4cG9ydHMgPSBjcmVhdGVUaWxlO1xuXG5mdW5jdGlvbiBjcmVhdGVUaWxlKGZlYXR1cmVzLCB6MiwgdHgsIHR5LCB0b2xlcmFuY2UsIGV4dGVudCwgbm9TaW1wbGlmeSkge1xuICAgIHZhciB0aWxlID0ge1xuICAgICAgICBmZWF0dXJlczogW10sXG4gICAgICAgIG51bVBvaW50czogMCxcbiAgICAgICAgbnVtU2ltcGxpZmllZDogMCxcbiAgICAgICAgbnVtRmVhdHVyZXM6IDAsXG4gICAgICAgIHNvdXJjZTogbnVsbFxuICAgIH07XG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBmZWF0dXJlcy5sZW5ndGg7IGkrKykge1xuICAgICAgICB0aWxlLm51bUZlYXR1cmVzKys7XG4gICAgICAgIGFkZEZlYXR1cmUodGlsZSwgZmVhdHVyZXNbaV0sIHoyLCB0eCwgdHksIHRvbGVyYW5jZSwgZXh0ZW50LCBub1NpbXBsaWZ5KTtcbiAgICB9XG4gICAgcmV0dXJuIHRpbGU7XG59XG5cbmZ1bmN0aW9uIGFkZEZlYXR1cmUodGlsZSwgZmVhdHVyZSwgejIsIHR4LCB0eSwgdG9sZXJhbmNlLCBleHRlbnQsIG5vU2ltcGxpZnkpIHtcblxuICAgIHZhciBnZW9tID0gZmVhdHVyZS5nZW9tZXRyeSxcbiAgICAgICAgdHlwZSA9IGZlYXR1cmUudHlwZSxcbiAgICAgICAgdHJhbnNmb3JtZWQgPSBbXSxcbiAgICAgICAgc3FUb2xlcmFuY2UgPSB0b2xlcmFuY2UgKiB0b2xlcmFuY2UsXG4gICAgICAgIGksIGosIHJpbmcsIHA7XG5cbiAgICBpZiAodHlwZSA9PT0gMSkge1xuICAgICAgICBmb3IgKGkgPSAwOyBpIDwgZ2VvbS5sZW5ndGg7IGkrKykge1xuICAgICAgICAgICAgdHJhbnNmb3JtZWQucHVzaCh0cmFuc2Zvcm1Qb2ludChnZW9tW2ldLCB6MiwgdHgsIHR5LCBleHRlbnQpKTtcbiAgICAgICAgICAgIHRpbGUubnVtUG9pbnRzKys7XG4gICAgICAgICAgICB0aWxlLm51bVNpbXBsaWZpZWQrKztcbiAgICAgICAgfVxuXG4gICAgfSBlbHNlIHtcblxuICAgICAgICAvLyBzaW1wbGlmeSBhbmQgdHJhbnNmb3JtIHByb2plY3RlZCBjb29yZGluYXRlcyBmb3IgdGlsZSBnZW9tZXRyeVxuICAgICAgICBmb3IgKGkgPSAwOyBpIDwgZ2VvbS5sZW5ndGg7IGkrKykge1xuICAgICAgICAgICAgcmluZyA9IGdlb21baV07XG5cbiAgICAgICAgICAgIC8vIGZpbHRlciBvdXQgdGlueSBwb2x5bGluZXMgJiBwb2x5Z29uc1xuICAgICAgICAgICAgaWYgKCFub1NpbXBsaWZ5ICYmICgodHlwZSA9PT0gMiAmJiByaW5nLmRpc3QgPCB0b2xlcmFuY2UpIHx8XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICh0eXBlID09PSAzICYmIHJpbmcuYXJlYSA8IHNxVG9sZXJhbmNlKSkpIHtcbiAgICAgICAgICAgICAgICB0aWxlLm51bVBvaW50cyArPSByaW5nLmxlbmd0aDtcbiAgICAgICAgICAgICAgICBjb250aW51ZTtcbiAgICAgICAgICAgIH1cblxuICAgICAgICAgICAgdmFyIHRyYW5zZm9ybWVkUmluZyA9IFtdO1xuXG4gICAgICAgICAgICBmb3IgKGogPSAwOyBqIDwgcmluZy5sZW5ndGg7IGorKykge1xuICAgICAgICAgICAgICAgIHAgPSByaW5nW2pdO1xuICAgICAgICAgICAgICAgIC8vIGtlZXAgcG9pbnRzIHdpdGggaW1wb3J0YW5jZSA+IHRvbGVyYW5jZVxuICAgICAgICAgICAgICAgIGlmIChub1NpbXBsaWZ5IHx8IHBbMl0gPiBzcVRvbGVyYW5jZSkge1xuICAgICAgICAgICAgICAgICAgICB0cmFuc2Zvcm1lZFJpbmcucHVzaCh0cmFuc2Zvcm1Qb2ludChwLCB6MiwgdHgsIHR5LCBleHRlbnQpKTtcbiAgICAgICAgICAgICAgICAgICAgdGlsZS5udW1TaW1wbGlmaWVkKys7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIHRpbGUubnVtUG9pbnRzKys7XG4gICAgICAgICAgICB9XG5cbiAgICAgICAgICAgIHRyYW5zZm9ybWVkLnB1c2godHJhbnNmb3JtZWRSaW5nKTtcbiAgICAgICAgfVxuICAgIH1cblxuICAgIGlmICh0cmFuc2Zvcm1lZC5sZW5ndGgpIHtcbiAgICAgICAgdGlsZS5mZWF0dXJlcy5wdXNoKHtcbiAgICAgICAgICAgIGdlb21ldHJ5OiB0cmFuc2Zvcm1lZCxcbiAgICAgICAgICAgIHR5cGU6IHR5cGUsXG4gICAgICAgICAgICB0YWdzOiBmZWF0dXJlLnRhZ3MgfHwgbnVsbFxuICAgICAgICB9KTtcbiAgICB9XG59XG5cbmZ1bmN0aW9uIHRyYW5zZm9ybVBvaW50KHAsIHoyLCB0eCwgdHksIGV4dGVudCkge1xuICAgIHZhciB4ID0gTWF0aC5yb3VuZChleHRlbnQgKiAocFswXSAqIHoyIC0gdHgpKSxcbiAgICAgICAgeSA9IE1hdGgucm91bmQoZXh0ZW50ICogKHBbMV0gKiB6MiAtIHR5KSk7XG4gICAgcmV0dXJuIFt4LCB5XTtcbn1cbiJdfQ==
