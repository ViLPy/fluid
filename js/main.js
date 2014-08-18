(function() {
    var canvas = document.getElementById('myCanvas'),
        ctx = canvas.getContext('2d');

    // Simulation step and overall time
    var dt = 0.1,
        time = 0;

    var width = parseInt(canvas.getAttribute('width')),
        height = parseInt(canvas.getAttribute('height'));

    var solver, timeOut, isPaused, lastTime;

    function init(demoCase) {
        var i,j;

        clearTimeout(timeOut);
        isPaused = false;
        time = 0;
        ctx.fillStyle = '#0000ff';

        solver = new FluidSolver(width, height, []);

        switch (demoCase) {
            case 'demo1':
                for (i = 30; i< 70; i++) {
                    for (j = 30; j<120; j++) {
                        solver.addParticle(new Particle(i,j, 0, 0));
                    }
                }
                break;

            case 'demo2':
                for (i = 30; i< 70; i++) {
                    for (j = 30; j<120; j++) {
                        solver.addParticle(new Particle(i,j, 20, 0));
                    }
                }
                break;

            case 'raindrop':
                for (i = 90; i< 110; i++) {
                    for (j = 270; j<300; j++) {
                        solver.addParticle(new Particle(i,j, 0, 0));
                    }
                }
                break;

            case 'flatsplash':
                for (i = 0; i< 200; i++) {
                    for (j = 280; j<300; j++) {
                        solver.addParticle(new Particle(i,j, 0, 0));
                    }
                }
                break;
        }

        lastTime = Date.now();
    }

    function loop() {
        render();
        var currentTime = Date.now();
        var deltaT = currentTime - lastTime;
        lastTime = currentTime;
        if (!isPaused) timeOut = setTimeout(loop, 1000*dt - deltaT);
    }

    function render() {
        var particlesCount = solver._particles.length;

        time += dt;

        solver.solve(dt);
        ctx.clearRect(0, 0, width, height);
        for (var i = 0; i < particlesCount; i++) {
            var p = solver._particles[i];
            renderParticle(p._x, height - p._y);
        }

        document.getElementById('time').innerHTML = (~~(time * 100)/100).toString();
    }

    function renderParticle(x, y) {
        ctx.beginPath();
        ctx.rect(x, y, 1, 1);
        ctx.fill();
    }


    function getDemoCase() {
        return document.getElementById('demoCase').value;
    }


    // Init controls
    document.getElementById('playSim').onclick = function() {
        if (isPaused) {
            isPaused = false;
            loop();
        } else {
            init(getDemoCase());
            loop();
        }
    };

    document.getElementById('pauseSim').onclick = function() {
        isPaused = true;
    };

    document.getElementById('stepForward').onclick = function() {
        render();
    };

    document.getElementById('stopSim').onclick = function() {
        init(getDemoCase());
    };
})();