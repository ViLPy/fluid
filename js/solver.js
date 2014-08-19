'use strict';

/**
 * FluidSolver based on 'Particle-based Viscoelastic Fluid Simulation' paper
 * by Simon Clavet, Philippe Beaudoin, and Pierre Poulin
 *
 * This is simplified version without springs between particles
 *
 * @param width - container width
 * @param height - container height
 * @constructor
 */
function FluidSolver(width, height) {
    this._particles = [];

    this._width = width;
    this._height = height;

    this._h = 3.4; // interaction radius
    this._rho0 = 15; // rest density
    this._k = 0.5; // stiffness parameter
    this._kNear = 5; // near stiffness parameter
    this._sigma = 0; // linear viscosity
    this._beta = 0.3; // quadratic viscosity

    // Grid for neighbor clustering
    this._NGrid = [];
    this._NGridWidth = Math.ceil(width / this._h);
    this._NGridHeight = Math.ceil(height / this._h);

    // Neighbor map, contains all neighbors for i-th particle
    this._Ni = {};


    // NGrid initial clustering
    var i;
    for (i = 0; i < (this._NGridWidth * this._NGridHeight); i++) {
        this._NGrid.push([]);
    }
}

FluidSolver.prototype.addParticle = function (particle) {
    this.updateParticleIndex(particle);
    this._particles.push(particle);
    this._NGrid[particle._NIndex].push(this._particles.length - 1);
};

FluidSolver.prototype.solve = function (dt) {
    var i, particle,
        particlesCount = this._particles.length,
        g = -9.8;

    // apply gravity
    for (i = 0; i < particlesCount; i++) {
        particle = this._particles[i];
        particle._vy += g * dt;
    }
    this.findNeighbors();

    this.applyViscosity(dt);

    for (i = 0; i < particlesCount; i++) {
        // save previous position
        particle = this._particles[i];
        particle._x0 = particle._x;
        particle._y0 = particle._y;

        // advance to predicted position
        particle._x += dt * particle._vx;
        particle._y += dt * particle._vy;
    }

    this.doubleDensityRelaxation(dt);

    for (i = 0; i < particlesCount; i++) {
        particle = this._particles[i];

        //use previous position to compute next velocity
        particle._vx = (particle._x - particle._x0) / dt;
        particle._vy = (particle._y - particle._y0) / dt;

        // collide with walls and update velocity after collision
        var mu1 = 0, mu2 = 0.9;
        if (particle._x < 0) {
            particle._x = 0;
            particle._vx *= mu1;
            particle._vy *= mu2;
        } else if (particle._x > this._width) {
            particle._x = this._width;
            particle._vx *= mu1;
            particle._vy *= mu2;
        }

        if (particle._y < 0) {
            particle._y = 0;
            particle._vx *= mu2;
            particle._vy *= mu1;
        } else if (particle._y > this._height) {
            particle._y = this._height;
            particle._vx *= mu2;
            particle._vy *= mu1;
        }

        // update NGrid indices if needed
        var indexPrev = particle._NIndex;

        // update index to get new value
        this.updateParticleIndex(particle);
        var indexNew = particle._NIndex;

        if (indexPrev != indexNew) {
            this._NGrid[indexPrev].splice(this._NGrid[indexPrev].indexOf(i), 1);
            this._NGrid[indexNew].push(i);
        }
    }
};

FluidSolver.prototype.doubleDensityRelaxation = function (dt) {
    var j, neighborIdx , neighbor, rij, rijLen, q, coeff,
        particlesCount = this._particles.length;

    for (var i = 0; i < particlesCount; i++) {
        var particle = this._particles[i];
        var rho = 0,
            rhoNear = 0;

        if (this._Ni[i]) {
            for (j = 0; j < this._Ni[i].length; j++) {
                neighborIdx = this._Ni[i][j];
                neighbor = this._particles[neighborIdx];

                rij = neighbor.subtractPositionVectorOf(particle);
                rijLen = vectorLength(rij);
                q = rijLen / this._h;

                if (q < 1) {
                    coeff = 1 - q;
                    rho += coeff * coeff;
                    rhoNear += coeff * coeff * coeff;
                }
            }
        }

        var P = this._k * (rho - this._rho0),
            PNear = this._kNear * rhoNear;

        var dx = {x: 0, y: 0};

        if (this._Ni[i]) {
            for (j = 0; j < this._Ni[i].length; j++) {
                neighborIdx = this._Ni[i][j];
                neighbor = this._particles[neighborIdx];

                rij = neighbor.subtractPositionVectorOf(particle);
                rijLen = vectorLength(rij);
                q = rijLen / this._h;

                if (q < 1) {
                    coeff = 1 - q;
                    var DCommon = dt * dt * (P * coeff + PNear * coeff * coeff);
                    var D = {x: rij.x * DCommon, y: rij.y * DCommon};

                    neighbor._x += D.x / 2;
                    neighbor._y += D.y / 2;

                    dx.x -= D.x / 2;
                    dx.y -= D.y / 2;
                }
            }
        }

        particle._x += dx.x;
        particle._y += dx.y;
    }
};

FluidSolver.prototype.applyViscosity = function (dt) {
    var keys = Object.keys(this._Ni),
        kTotal = keys.length;

    for (var k = 0; k < kTotal; k++) {
        var i = keys[k];
        if (this._Ni[i]) {
            var particle = this._particles[i];
            for (var j = 0; j < this._Ni[i].length; j++) {
                var neighborIdx = this._Ni[i][j];
                var neighbor = this._particles[neighborIdx];

                var rij = neighbor.subtractPositionVectorOf(particle);
                var rijLen = vectorLength(rij);
                var q = rijLen / this._h;

                if (q < 1) {
                    var dv = {x: particle._vx - neighbor._vx, y: particle._vy - neighbor._vy};
                    var rijNorm = unitVector(rij, rijLen);
                    var u = dotProductVectors(dv, rijNorm);

                    if (u > 0) {
                        var common = dt * (1 - q) * (this._sigma * u + this._beta * u * u);
                        var I = {x: rijNorm.x * common, y: rijNorm.y * common};
                        particle._vx -= I.x / 2;
                        particle._vy -= I.y / 2;

                        neighbor._vx += I.x / 2;
                        neighbor._vy += I.y / 2;
                    }
                }
            }
        }
    }
};

// Neighbor calculation
FluidSolver.prototype.updateParticleIndex = function (particle) {
    particle._NX = ~~(particle._x / this._h);
    particle._NY = ~~(particle._y / this._h);
    particle._NIndex = particle._NX + this._NGridWidth * particle._NY
};

FluidSolver.prototype.getNeighbors = function (particle) {
    var gridWidth = this._NGridWidth,
        gridHeight = this._NGridHeight;

    var NGrid = this._NGrid;

    function getByIndex(i, j) {
        if (i >= 0 && i < gridWidth
            && j >= 0 && j < gridHeight) {
            return NGrid[i + gridWidth * j];
        } else {
            return [];
        }
    }

    var i = particle._NX,
        j = particle._NY;

    return this._NGrid[particle._NIndex].concat(
        getByIndex(i - 1, j - 1), getByIndex(i - 1, j), getByIndex(i - 1, j + 1),
        getByIndex(i, j - 1), getByIndex(i, j + 1),
        getByIndex(i + 1, j - 1), getByIndex(i + 1, j), getByIndex(i + 1, j + 1));
};

FluidSolver.prototype.findNeighbors = function () {
    var particlesCount = this._particles.length,
        distanceToNeighbourSq = this._h * this._h;

    this._Ni = {};

    for (var i = 0; i < particlesCount; i++) {
        var particle = this._particles[i];
        var neighborIndexes = this.getNeighbors(particle);
        var neighborLength = neighborIndexes.length;
        for (var j = 0; j < neighborLength; j++) {
            var neighborIndex = neighborIndexes[j];
            if (i <= neighborIndex) continue;
            var particleToCheck = this._particles[neighborIndex];
            var distanceSq = particle.distanceSquaredTo(particleToCheck);

            if (distanceSq <= distanceToNeighbourSq) {
                if (!this._Ni[i]) this._Ni[i] = [];
                if (!this._Ni[neighborIndex]) this._Ni[neighborIndex] = [];

                this._Ni[i].push(neighborIndex);
                this._Ni[neighborIndex].push(i);
            }
        }

    }
};