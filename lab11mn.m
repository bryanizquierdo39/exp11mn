clear; clc; close all;

mu = 1.5; % coeficiente de amortiguamiento no lineal
t0 = 0; % Tiempo inicial
tf = 20; % Tiempo final
h = 0.05; % Paso
t = t0:h:tf; % Vector de tiempo
N = length(t); % Número de pasos

x0 = 2.0; % x(0) = desplazamiento inicial
v0 = 0; % v(0) = x'(0) = velocidad inicial
Y = zeros(2, N); % Y(1,:) = x(t), Y(2,:) = v(t)
Y(:,1) = [x0; v0];

% función del sistema van der Pol
% Y = [x; v]
f = @(t, Y) [Y(2);
    mu*(1 - Y(1)^2)*Y(2) - Y(1)];

% metodo RK4
for i = 1:N-1
    k1 = f(t(i), Y(:,i));
    k2 = f(t(i) + h/2, Y(:,i) + h/2*k1);
    k3 = f(t(i) + h/2, Y(:,i) + h/2*k2);
    k4 = f(t(i) + h, Y(:,i) + h*k3);

    Y(:,i+1) = Y(:,i) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
end

x = Y(1,:); % Desplazamiento
v = Y(2,:); % Velocidad

% Gráfica 1: Evolución temporal x(t) y v(t)
figure(1);
plot(t, x, 'b-', 'LineWidth', 1.5); hold on;
plot(t, v, 'r--', 'LineWidth', 1.5);
xlabel('Tiempo t [s]');
ylabel('Amplitud');
title('Evolución temporal - Ecuación de van der Pol');
legend('Desplazamiento x(t)', 'Velocidad v(t)', 'Location', 'best');
grid on;

% Gráfica 2: Diagrama de fases x vs v (ciclo límite)
figure(2);
plot(x, v, 'g-', 'LineWidth', 1.5);
xlabel('Desplazamiento x');
ylabel('Velocidad v = x''');
title('Diagrama de Fases - Ciclo Límite van der Pol \mu = 1.5');
grid on;
axis equal;

% primeras 10 iteraciones y estados finales
fprintf('\n=== TABLA DE DATOS RK4 - van der Pol ===\n');
fprintf('h = %.2f s, mu = %.1f, t = [0, %d] s\n', h, mu, tf);
fprintf('Condiciones iniciales: x(0) = %.1f, v(0) = %.1f\n', x0, v0);

fprintf('%-8s %-15s %-15s\n', 'Iter', 't [s]', 'x(t)', 'v(t)');
fprintf('------------------------------------------------\n');

% Primeras 10 iteraciones
for i = 1:min(10, N)
    fprintf('%-8d %-15.2f %-15.6f %-15.6f\n', i-1, t(i), x(i), v(i));
end

fprintf('...\n');

% estados finales, últimas 5
for i = max(1, N-4):N
    fprintf('%-8d %-15.2f %-15.6f %-15.6f\n', i-1, t(i), x(i), v(i));
end
fprintf('------------------------------------------------\n');
