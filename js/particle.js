/**
 * Particle model
 * @constructor
 */
function Particle(x,y, vx, vy) {
    // real position
    this._x = x;
    this._y = y;

    // past position
    this._x0 = x;
    this._y0 = y;

    // velocity
    this._vx = vx || 0;
    this._vy = vy || 0;
}

Particle.prototype.distanceSquaredTo = function(p2) {
    var dx = this._x - p2._x;
    var dy = this._y - p2._y;
    return dx*dx + dy*dy;
};

Particle.prototype.subtractPositionVectorOf = function(p2) {
    return {x: this._x - p2._x, y: this._y - p2._y};
};