/* Misc. vector utils */
function vectorLengthSquared(r) {
    return r.x * r.x + r.y * r.y;
}

function vectorLength(r) {
    return Math.sqrt(vectorLengthSquared(r));
}

function dotProductVectors(p1, p2) {
    return p1.x * p2.x + p1.y * p2.y;
}

function unitVector(v, len) {
    return {x: v.x / len, y: v.y / len};
}