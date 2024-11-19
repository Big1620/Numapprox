let myChart = null;

function evaluateFunction(func, x) {
    const scope = { x: x };
    return math.evaluate(func, scope);
}

// Improved Bisection Method
function bisectionMethod(f, a, b, tol, maxIter) {
    let iterations = [];
    let fa = evaluateFunction(f, a);
    let fb = evaluateFunction(f, b);

    if (fa * fb > 0) {
        throw new Error("La función debe cambiar de signo en los extremos del intervalo [a, b].");
    }

    let i = 0;
    let c_prev = null;
    let error = tol + 1;

    while (error > tol && i < maxIter) {
        let c = (a + b) / 2;
        let fc = evaluateFunction(f, c);

        if (c_prev !== null) {
            error = Math.abs(c - c_prev);
        }

        iterations.push({
            iteration: i + 1,
            value: c,
            fc: fc,
            error: error
        });

        if (fc === 0) break;
        
        if (fa * fc < 0) {
            b = c;
        } else {
            a = c;
            fa = fc;
        }

        c_prev = c;
        i++;
    }

    return iterations;
}

// Improved Newton-Raphson Method
function newtonMethod(f, x0, tol, maxIter) {
    let iterations = [];
    let i = 0;
    let x = x0;
    let error = tol + 1;

    const calculateDerivative = (func, x, h = 1e-5) => {
        const fx = evaluateFunction(func, x);
        const fxh = evaluateFunction(func, x + h);
        return (fxh - fx) / h;
    };

    while (error > tol && i < maxIter) {
        let fx = evaluateFunction(f, x);
        let derivative = calculateDerivative(f, x);

        if (Math.abs(derivative) < 1e-10) {
            throw new Error("Derivada cercana a cero. No se puede continuar.");
        }

        let x_prev = x;
        x = x - fx / derivative;

        error = Math.abs(x - x_prev);

        iterations.push({
            iteration: i + 1,
            value: x,
            fc: fx,
            error: error
        });

        if (Math.abs(fx) < tol) break;
        i++;
    }

    return iterations;
}

// Improved Secant Method
function secantMethod(f, x0, x1, tol, maxIter) {
    let iterations = [];
    let i = 0;
    let error = tol + 1;

    while (error > tol && i < maxIter) {
        let fx0 = evaluateFunction(f, x0);
        let fx1 = evaluateFunction(f, x1);

        if (Math.abs(fx1 - fx0) < 1e-10) {
            throw new Error("Problemas de convergencia en el método de la secante");
        }

        let x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
        let fx2 = evaluateFunction(f, x2);

        error = Math.abs(x2 - x1);

        iterations.push({
            iteration: i + 1,
            value: x2,
            fc: fx2,
            error: error
        });

        x0 = x1;
        x1 = x2;

        if (Math.abs(fx2) < tol) break;
        i++;
    }

    return iterations;
}
// Interpolación Lineal
function linearInterpolation(points, x) {
    // Ordenar puntos por x para asegurar la interpolación correcta
    points.sort((a, b) => a[0] - b[0]);

    // Encontrar los dos puntos más cercanos
    for (let i = 0; i < points.length - 1; i++) {
        const [x0, y0] = points[i];
        const [x1, y1] = points[i + 1];

        // Verificar si x está dentro del rango
        if (x >= x0 && x <= x1) {
            // Interpolación lineal
            return y0 + ((y1 - y0) / (x1 - x0)) * (x - x0);
        }
    }

    // Manejar casos fuera de rango
    if (x < points[0][0]) {
        // Extrapolar al inicio
        const [x0, y0] = points[0];
        const [x1, y1] = points[1];
        return y0 - ((y1 - y0) / (x1 - x0)) * (x0 - x);
    }

    if (x > points[points.length - 1][0]) {
        // Extrapolar al final
        const [x0, y0] = points[points.length - 2];
        const [x1, y1] = points[points.length - 1];
        return y1 + ((y1 - y0) / (x1 - x0)) * (x - x1);
    }

    throw new Error(`No se pudo interpolar para x = ${x}`);
}

// Interpolación Polinómica - Mejorada (Método de Lagrange)
function polynomialInterpolation(points, x) {
    let result = 0;

    for (let i = 0; i < points.length; i++) {
        let term = points[i][1];
        for (let j = 0; j < points.length; j++) {
            if (i !== j) {
                // Cálculo del término de Lagrange
                term *= (x - points[j][0]) / (points[i][0] - points[j][0]);
            }
        }
        result += term;
    }

    return result;
}

// Regla del Trapecio
function trapezoidRule(f, a, b, n) {
    // Validar número de intervalos
    if (n <= 0) throw new Error("El número de intervalos debe ser positivo");

    const h = (b - a) / n;
    let sum = 0.5 * (evaluateFunction(f, a) + evaluateFunction(f, b));

    for (let i = 1; i < n; i++) {
        const x = a + i * h;
        sum += evaluateFunction(f, x);
    }

    return sum * h;
}

// Regla de Simpson - Mejorada
function simpsonRule(f, a, b, n) {
    // Asegurar número par de intervalos
    if (n % 2 !== 0) n++;

    const h = (b - a) / n;
    let sum = evaluateFunction(f, a) + evaluateFunction(f, b);

    // Términos pares
    for (let i = 1; i < n; i += 2) {
        sum += 4 * evaluateFunction(f, a + i * h);
    }

    // Términos impares
    for (let i = 2; i < n - 1; i += 2) {
        sum += 2 * evaluateFunction(f, a + i * h);
    }

    return (sum * h) / 3;
}

// Gauss-Legendre - Mejorada
function gaussLegendre(f, a, b) {
    // Puntos y pesos para cuadratura de Gauss-Legendre
    const roots = [
        -0.577350269189626,  // -1/√3
        0.577350269189626    // 1/√3
    ];
    const weights = [1, 1];

    const midpoint = (a + b) / 2;
    const halfRange = (b - a) / 2;

    return halfRange * roots.reduce((sum, root, i) => {
        const x = midpoint + halfRange * root;
        return sum + weights[i] * evaluateFunction(f, x);
    }, 0);
}
// Actualización de Gráfica
function updateChart(f, a, b, iterations, type = "function", points = []) {
    const ctx = document.getElementById('chart').getContext('2d');
    if (myChart) myChart.destroy();
    
    let datasets = [];

    // Manejar gráficas de interpolación
    if (type === "interpolation") {
        datasets.push({
            label: "Puntos de Interpolación",
            data: points.map(([x, y]) => ({ x, y })),
            backgroundColor: "rgb(52, 152, 219)",
            borderColor: "rgb(52, 152, 219)",
            type: 'scatter',
            pointRadius: 7
        });
    } 
    // Manejar gráficas de función y métodos numéricos
    else if (type === "function" || type === "integration") {
        // Graficar función base
        if (f) {
            const step = (b - a) / 200;
            const xValues = Array.from({ length: 201 }, (_, i) => a + i * step);
            const yValues = xValues.map(x => evaluateFunction(f, x));

            datasets.push({
                label: "f(x)",
                data: xValues.map((x, i) => ({ x, y: yValues[i] })),
                borderColor: "rgb(52, 152, 219)",
                fill: false
            });
        }

        // Graficar puntos de iteración
        if (iterations && iterations.length) {
            datasets.push({
                label: "Iteraciones",
                data: iterations.map(iter => ({
                    x: iter.value,
                    y: evaluateFunction(f, iter.value)
                })),
                borderColor: "rgb(231, 76, 60)",
                backgroundColor: "rgb(231, 76, 60)",
                type: 'scatter',
                pointRadius: 7
            });
        }
    }

    // Configuración final de la gráfica
    myChart = new Chart(ctx, {
        type: "line",
        data: { datasets },
        options: { 
            responsive: true,
            scales: {
                x: { 
                    type: 'linear', 
                    position: 'bottom',
                    title: {
                        display: true,
                        text: 'x'
                    }
                },
                y: { 
                    type: 'linear',
                    title: {
                        display: true,
                        text: 'f(x)'
                    }
                }
            },
            plugins: {
                title: {
                    display: true,
                    text: 'Visualización de Método Numérico'
                }
            }
        }
    });
}

// Función Principal
function calculate() {
    try {
        const method = document.getElementById('methodSelect').value;

        if (["bisection", "newton", "secant"].includes(method)) {
            const f = document.getElementById('function').value.trim();
            if (!f) throw new Error("Debe ingresar una función válida.");

            const intervalInput = document.getElementById('interval').value.trim();
            const [a, b] = intervalInput.split(',').map(Number);

            if (isNaN(a) || isNaN(b)) throw new Error("Debe ingresar un intervalo válido [a, b].");

            let iterations;
            if (method === "bisection") {
                const tol = parseFloat(document.getElementById('tolerance').value);
                const maxIter = parseInt(document.getElementById('maxIterations').value);
                iterations = bisectionMethod(f, a, b, tol, maxIter);
            } else if (method === "newton") {
                const x0 = parseFloat(document.getElementById('initialValue').value);
                const tol = parseFloat(document.getElementById('tolerance').value);
                iterations = newtonMethod(f, x0, tol, 100);
            } else if (method === "secant") {
                const tol = parseFloat(document.getElementById('tolerance').value);
                const maxIter = parseInt(document.getElementById('maxIterations').value);
                iterations = secantMethod(f, a, b, tol, maxIter);
            }

            updateChart(f, a, b, iterations, "function");
            updateTable(f, iterations);

        } else if (["linear-interpolation", "polynomial-interpolation"].includes(method)) {
            const pointsInput = document.getElementById('points').value.trim();
            const x = parseFloat(document.getElementById('xValue').value);

            if (!pointsInput || isNaN(x)) {
                throw new Error("Debe ingresar puntos válidos y un valor de x.");
            }

            const points = pointsInput
                .split('\n')
                .map(line => line.split(',').map(Number))
                .filter(pair => pair.length === 2);

            if (points.length < 2) {
                throw new Error("Se necesitan al menos dos puntos para interpolar.");
            }

            let result;
            if (method === "linear-interpolation") {
                result = linearInterpolation(points, x);
            } else if (method === "polynomial-interpolation") {
                result = polynomialInterpolation(points, x);
            }

            updateChart(null, null, null, null, "interpolation", points);
            document.getElementById('statusBar').textContent = `Resultado de interpolación: f(${x}) = ${result}`;

        } else if (["trapezoid", "simpson", "gauss-legendre"].includes(method)) {
            const f = document.getElementById('function').value.trim();
            const intervalInput = document.getElementById('interval').value.trim();
            const [a, b] = intervalInput.split(',').map(Number);
            const n = parseInt(document.getElementById('intervalCount').value);

            if (isNaN(a) || isNaN(b) || a >= b) throw new Error("Intervalo [a, b] inválido.");
            if (isNaN(n) || n <= 0) throw new Error("El número de intervalos debe ser un entero positivo.");

            let result;
            if (method === "trapezoid") result = trapezoidRule(f, a, b, n);
            else if (method === "simpson") result = simpsonRule(f, a, b, n);
            else if (method === "gauss-legendre") result = gaussLegendre(f, a, b);

            updateChart(f, a, b, null, "integration");
            document.getElementById('statusBar').textContent = `Resultado de integración: ${result}`;
        }
    } catch (error) {
        document.getElementById('statusBar').textContent = `Error: ${error.message}`;
    }
}

function updateTable(f, iterations) {
    const tbody = document.getElementById('resultsBody');
    tbody.innerHTML = ''; // Limpia la tabla antes de agregar los nuevos resultados

    iterations.forEach((iter, index) => {
        const row = tbody.insertRow();
        const functionValue = evaluateFunction(f, iter.value);

        row.insertCell(0).textContent = iter.iteration;
        row.insertCell(1).textContent = iter.value.toFixed(6);
        row.insertCell(2).textContent = functionValue.toFixed(6);
        row.insertCell(3).textContent = iter.error.toFixed(6);
        row.insertCell(4).textContent =
            index > 0 ? Math.abs(iterations[index].value - iterations[index - 1].value).toFixed(6) : '-';
    });
}

function updateMethodInputs() {
    const method = document.getElementById('methodSelect').value;
    const dynamicParams = document.getElementById('dynamic-parameters');
    dynamicParams.innerHTML = ''; // Limpia los parámetros anteriores

    if (["bisection", "newton", "secant"].includes(method)) {
        dynamicParams.innerHTML = `
            <label>Función f(x):</label>
            <input type="text" id="function" placeholder="Ejemplo: x^2 - 2">
            <label>Intervalo [a, b]:</label>
            <input type="text" id="interval" placeholder="Ejemplo: 1, 2">
            ${method !== "newton" ? `
            <label>Tolerancia:</label>
            <input type="number" id="tolerance" value="0.0001" step="any">
            <label>Iteraciones máximas:</label>
            <input type="number" id="maxIterations" value="100">` : ''}
            ${method === "newton" ? `
            <label>Valor inicial x₀:</label>
            <input type="number" id="initialValue" placeholder="Ejemplo: 1.5">` : ''}
        `;
    } else if (["linear-interpolation", "polynomial-interpolation"].includes(method)) {
        dynamicParams.innerHTML = `
            <label>Puntos (x, y):</label>
            <textarea id="points" rows="4" placeholder="Formato: 1,2\n2,4\n3,6"></textarea>
            <label>Valor x a interpolar:</label>
            <input type="number" id="xValue" step="any">
        `;
    } else if (["trapezoid", "simpson", "gauss-legendre"].includes(method)) {
        dynamicParams.innerHTML = `
            <label>Función f(x):</label>
            <input type="text" id="function" placeholder="Ejemplo: x^2">
            <label>Intervalo [a, b]:</label>
            <input type="text" id="interval" placeholder="Ejemplo: 0, 2">
            <label>Número de intervalos (n):</label>
            <input type="number" id="intervalCount" value="10" step="1">
        `;
    }
}

function clearAll() {
    // Limpiar campos de entrada
    const dynamicParams = document.getElementById('dynamic-parameters');
    const inputs = dynamicParams.querySelectorAll('input, textarea');
    inputs.forEach(input => input.value = '');

    // Reiniciar la tabla de resultados
    const tbody = document.getElementById('resultsBody');
    tbody.innerHTML = '';

    // Reiniciar la barra de estado
    const statusBar = document.getElementById('statusBar');
    statusBar.textContent = 'Estado: Listo para calcular | Último cálculo: --:--:-- | Error: 0.00%';

    // Destruir la gráfica actual
    if (myChart) {
        myChart.destroy();
        myChart = null;
    }

    // Restablecer el selector de método a su valor predeterminado
    const methodSelect = document.getElementById('methodSelect');
    methodSelect.selectedIndex = 0;
    
    // Actualizar los parámetros de entrada según el método predeterminado
    updateMethodInputs();
}
