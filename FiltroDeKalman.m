%Definición del modelo
syms x1 x2 u T

epsilon = 0.1;

%Ecuaciones
ode1 = x2;
ode2 = -x1 + epsilon * (1-x1^2)*x2;
y = x1;
odes = [ode1; ode2]; vars = [x1, x2];

%Linealización
A = jacobian(odes,vars); B = jacobian(odes, u); C = jacobian(y,vars); 
sizeB = size(B); sizeC = size(C);D = zeros(sizeC(1),sizeB(2));



A = eye(2) + A*T;
%A(T,x1,x2): función de la matriz del modelo linealizado y discretizado
A = matlabFunction(A);
Ts = 0.1; % tiempo de muestreo T.
tiempo = 0:Ts:50; % periodo de tiempo a usar.
x0 = [1 1]; % Condición inicial

%%
% Entradas
u = sin(5*tiempo) + randn(1,length(tiempo));
%u = ones(1, length(tiempo));

%Visualizamos la entrada
plot(tiempo, u)
title('Entrada');
xlabel('k');
ylabel('u (k)');
grid on

%%
[t, x] = ode45(@(t,x) ode(t, x, u, tiempo, odes), tiempo, x0);
y = x(:,1);
%%
%Condiciones iniciales
%M : Matriz de covarianzas del ruido del estado
%R : Matriz de covarianzas del ruido de la entrada
%P0 : estimación inicial de nuestro pronóstico
%x0 : valores iniciales de los estados (no importa)
%Ts : tiempo de muestreo
M = 0.01*eye(2); R = 0.5; P0 = M; x0 = [1 1]'; Ts = 0.1;

%Cosas que conocemos
%PTraza : error que se esta cometiendo
%KNorm : norma de la K de Kalman
n = length(y); xAprox = zeros(n,2); PTrace = zeros(n,1); KNorm = zeros(n,1);
xK = x0; PK = P0;

%%


for k = 1:n
    Anew = A(Ts, xK(1), xK(2));
    % Predicción:
    xPrediction = [xK(1) + Ts*(xK(1)*xK(2)-xK(1));
                   xK(2) + Ts*(xK(1)^2 - xK(2) + u(k))];
    PPrediction = Anew*PK*Anew' + M;
    % Corrección:
    K = PPrediction*C'*inv(C*PPrediction*C' + R);
    xK = xPrediction + K*(y(k) - C*xPrediction);
    PK = (eye(2) - K*C)*PPrediction;
    % Estimación:
    xAprox(k,:) = xK';
    % Medidas:
    PTrace(k) = trace(PK);
    KNorm(k) = norm(K);
end

%%
plot(t, x(:,2), t, x_est(:,2))
legend({'Modelo original', 'Estimación Kalman Extendido'}, 'Location', 'best');
xlabel('k');
ylabel('x2 (k)');
title('Comparación modelo original con estimación kalman extendido');
grid on

