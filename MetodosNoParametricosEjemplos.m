%función de transferencia
G1 = tf(2,[1 4])
%respuesta temporal a una entrada de escalón unitario
tiempofinal = 5;
clf
figure(1)
step(G1,tiempofinal)
[y,t] = step(G1,tiempofinal);
inicio = find(y > 0, 1);%identificamos el dato a partir del cual el modelo responde a la entrada
yvalores = y(inicio:end);
tvalores = t(inicio:end);
k = yvalores(end)-y(1) %valor de k según la fórmula
k = 0.5;
[p,s] = polyfit(tvalores,log(1-yvalores/(k)),1);%regresión lineal al logaritmo de y
T = -1/p(1);
tau = -p(2)/p(1);%valores de parámetros  T y tau según las fórmulas
tfEstimada = tf(k,[T 1],'InputDelay',abs(tau))%función de transferencia resultante
figure(1)
hold on
step(tfEstimada)
xline(T+tau, '--')
yline(k*0.632,'--')
legend({'Original','Regresión'},'Location','best')

%%
%función de transferencia
G3 = zpk([-3.9],[-4, -10, -(0.5 + 1.32287*i),-( 0.5-1.32287*i)],5,'InputDelay',1)
%respuesta temporal a una entrada de escalón unitario
tiempofinal = 20;
figure(3)
clf
step(G3,tiempofinal)
grid on
legend({'G(s)'})
[y,t] = step(G3);%recolección de datos por recomendación de MatLab
inicio = find(y > 0, 1);%dato a partir del cual el sistema reacciona
tau = t(inicio)%tiempo a partir del cual el sistema reacciona
deltaY = y(end) - y(1)%delta de y por definición
Mp = max(y) - y(end);%sobreimpulso máximo pode definición
Tp = t(y == max(y)) - tau;%tiempo de pico, se le debe restar el retardo
zeta = 1/sqrt(1+(pi/log(Mp/deltaY))^2)%valor de z (diapositivas)
%zeta = -log(Mp)/sqrt(pi^2+log(Mp)^2)%valor de z (Soderstrom)
omega0 = pi/(Tp*sqrt(1-zeta^2))%valor de omega0 por fórmulas
tfEstimada=tf(deltaY*omega0^2,[1 2*zeta*omega0 omega0^2],'InputDelay',tau)%función de transferencia estimada
figure(3)
clf
step(G3)
hold on
step(tfEstimada)

legend({'Original','Estimado'},'Location','best')