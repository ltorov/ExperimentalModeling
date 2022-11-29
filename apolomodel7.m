%%
%1. Define equations for the model and load timeseries data

syms  Ms(t) Me(t) Mi(t) Hs(t) He(t) Hi(t) Hr(t) Hit(t)

syms lambda beta_m mu_m theta_m mu_h beta_h theta_h...
    gamma_h

%pulse(t)=1;
H=Hs+He+Hi+Hr;
M=Ms+Me+Mi;
ode1 = diff(Ms) == lambda - beta_m*Hi*Ms/H - (mu_m)*Ms;
ode2 = diff(Me) == beta_m*Hi*Ms/H - (theta_m+mu_m)*Me;
ode3 = diff(Mi) == theta_m*Me - mu_m*Mi;
ode4 = diff(Hs) == -beta_h*Mi*Hs/M + (He+Hi+Hr)*mu_h;
ode5 = diff(He) == beta_h*Mi*Hs/M - (theta_h+mu_h)*He;
ode6 = diff(Hi) == theta_h*He - (gamma_h+mu_h)*Hi;
ode7 = diff(Hr) == gamma_h*Hi - mu_h*Hr;
ode8 = diff(Hit) == theta_h*He;
odes=[ode1; ode2 ;ode3; ode4; ode5 ;ode6; ode7; ode8];
vars=[Hit Hi Me Hr Hs He Ms Mi];
opts = odeset('NonNegative',1:8);

%2. Data preparation

load Range7.mat
[T,~]=gsua_dpmat(odes,vars,[0 150],'7m','output',1,'opt',opts,'Range',Range);
load('DataBello_full.mat');
ydata=DataBello.cases(153:303)';
xdata=linspace(0,length(ydata)-1,length(ydata));
ydata2=ydata;
for i=1:length(ydata)
ydata2(i)=sum(ydata(1:i));
end

%%
%3.Analisis de sensibilidad

%3.1. Data prep

[T,~]=gsua_dpmat(odes,vars,[0 150],'7m','output',1,'opt',opts,'Range',Range);
T.Properties.CustomProperties.output = 1;
M = gsua_dmatrix(T,100);

ynom = gsua_eval(T.Nominal,T);

plot(cumsum(ydata),'b')
title('Original Acumulated Human Infections (Hit)')
xlabel('Weeks')
ylabel('Cases')
savefig('iteration1/figures/NominalValues.fig')

%%


%3.2. Sensibility analysis

Tsa = gsua_sa(M,T,'parallel', false, 'SensMethod', 'Xiao', 'ynom', ynom);

gsua_plot('Bar',Tsa,Tsa.STi)
savefig('iteration1/figures/SensibilityAnalisis.fig')

%3.3. Confiabilidad

c = sum(Tsa.Si)/sum(abs(Tsa.Si))

%%
%3.4 Fixing parameters

%According to the sensitivity analysis, we will fix the 4 parameters with
%the lower sensitivity indexes (with the nominal values).

%These are He0 and Me0, añadir Ms0

vars=[Hit Hi He Ms Me Hr Hs Mi];
HeO = T.Nominal('He0'); MsO = T.Nominal('Ms0'); MeO = T.Nominal('Me0');

RangeTemp = Range;
Range(3,:) = HeO; 
Range(4,:) = MsO;
Range(5,:) = MeO;
Range(6:8,:) = [RangeTemp(4,:); RangeTemp(5,:); RangeTemp(8,:)];

[T,~]=gsua_dpmat(odes,vars,[0 150],'7m','output',1,'opt',opts,'Range',Range);
T.Properties.CustomProperties.output = 1;
M = gsua_dmatrix(T,100);

%%
%3.5. Parameter estimation

%(takes a lot of time)

opt = optimoptions('lsqcurvefit','UseParallel',false,'Display','iter');
[T7,res] = gsua_pe(T,xdata,ydata2,'N',100,'opt',opt); %mover N ( 100 buena, 50 aceptable, 10 :( )
save('iteration1/values/Results7.mat','T7','res','xdata','ydata2')

%%
%3.6 Identifiability analysis 
%Cogemos el 30 % para solo tomar una de las familias
th = sum(res<res(1)*1.3)
y7 = gsua_eval(T7.Estlsqc(:,1:th),T7,xdata,ydata2);
savefig('iteration1/figures/curves.fig')

%%
bestest = y7(1,:);
trend = [bestest(1),bestest(2:end)-bestest(1:end-1)];
plot(trend,'b')
hold on
plot(diff(ydata2),'r')
title('Estimated vs Real Weekly Infections')
xlabel('Weeks')
ylabel('Cases')
legend({'Estimated','Real'})
savefig('iteration1/figures/EstimatedvsReal.fig')

%%
T7.Nominal = T7.Estlsqc(:,1);
%res : funciones de costo
th = sum(res<res(1)*1.3) % threshold, cambiar el 1.01
%th = 476 %el parametro que habia dado el profe, a nosotros nos da 1????

T7_ia = gsua_ia(T7,T7.Estlsqc(:,1:th), true, false, true); 
savefig('iteration1/figures/CorrelationsDiag.fig')
save('iteration1/values/T.mat','T')
%% 
% Segunda iteracion de los calculos


%En teoria deberiamos aumentar los rangos pero eso lleva a un ciclo
%infinito de aumentar rangos, entonces no lo haremos

%5.1 Aumentamos los rangos de thetaH
T.Range('theta_h',1) = 0;
T.Range('lambda',1) = 800;
T.Range('gamma_h',1) = 0;


%Realizamos analisis otra vez
Tsa = gsua_sa(M,T,'parallel', false, 'SensMethod', 'Xiao', 'ynom', ynom);

gsua_plot('Bar',Tsa,Tsa.STi)
savefig('iteration2/figures/SensibilityAnalisis.fig')


% Confiabilidad

c = sum(Tsa.Si)/sum(abs(Tsa.Si))

%%
%5.2 estimation
opt = optimoptions('lsqcurvefit','UseParallel',true,'Display','iter');
[T7_1,res_1] = gsua_pe(T,xdata,ydata2,'N',100,'opt',opt); 
save('iteration2/values/Results7.mat','T7_1','res_1','xdata','ydata2')

%%
%5.2 Identifiability analysis 
%Como hay solo una familia se cogen todas las estimaciones
th_1 = sum(res_1<res_1(1)*1.15)
y7_1 = gsua_eval(T7_1.Estlsqc(:,1:end),T7_1,xdata,ydata2);
savefig('iteration2/figures/curves.fig')

%%
bestest_1 = y7_1(1,:);
trend = [bestest_1(1),bestest_1(2:end)-bestest_1(1:end-1)];
plot(trend,'b')
hold on
plot(diff(ydata2),'r')
title('Estimated vs Real Weekly Infections')
xlabel('Weeks')
ylabel('Cases')
legend({'Estimated','Real'})
savefig('iteration2/figures/EstimatedvsReal.fig')

%%
T7_1.Nominal = T7_1.Estlsqc(:,1);
%res : funciones de costo

T7_1ia = gsua_ia(T7_1,T7_1.Estlsqc(:,1:end), false, true); 
savefig('iteration2/figures/CorrelationsDiag.fig')

save('iteration2/values/T.mat','T')

%%
% Tercera Iteracion
%Fijamos Hs_0


%vars=[Hit Hi He Ms Me Hr Hs Mi];
%vars=[Hit Hi He Ms Me Hs Hr Mi];
Betah = T7.Nominal('beta_h'); 

Range(9,:) = Betah;

[T,~]=gsua_dpmat(odes,vars,[0 150],'7m','output',1,'opt',opts,'Range',Range);
T.Properties.CustomProperties.output = 1;
M = gsua_dmatrix(T,100);

%%
%5.2 estimation


opt1 = optimoptions('lsqcurvefit','UseParallel',false,'Display','iter');
[T7_1,res_1] = gsua_pe(T,xdata,ydata2,'N',100,'opt',opt1); %mover N ( 100 buena, 50 aceptable, 10 :( )
save('iteration2/values/Results7.mat','T7_1','res_1','xdata','ydata2')


%%
%5.2 Identifiability analysis 
th_2 = sum(res_2<res_2(1)*1.02)
y7_2 = gsua_eval(T7_2.Estlsqc(:,1:th_2),T7_2,xdata,ydata2);
savefig('iteration3/figures/curves.fig')

%%
bestest_2 = y7_2(1,:);
trend = [bestest_2(1),bestest_2(2:end)-bestest_2(1:end-1)];
plot(trend,'b')
hold on
plot(diff(ydata2),'r')
title('Estimated vs Real Weekly Infections')
xlabel('Weeks')
ylabel('Cases')
legend({'Estimated','Real'})
savefig('iteration3/figures1/EstimatedvsReal.fig')

%%
T7_2.Nominal = T7_2.Estlsqc(:,1);
%res : funciones de costo
th_2 = sum(res_2<res_2(1)*1.02) % threshold, cambiar el 1.01
%th = 476 %el parametro que habia dado el profe, a nosotros nos da 1????

T7_2 = gsua_ia(T7_2,T7_2.Estlsqc(:,1:th_2), false, true); 
savefig('iteration3/figures/Correlations.fig')
save('iteration3/values/T.mat','T')
%%
%4.4 Uncertainty analysis

%analisis de incertidumbre con T7_2. veamos que las curvas hagan una banda 
%al lado de los datos que estamos estimando. Si la banda es chevere,
%terminamos, si no, fijar algun parametro. 

Ua = gsua_ua(M, T7, 'parallel', false, 'ynom',ydata2);
savefig('iteration3/figures/Montecarlo.fig')

%%
% Cuarta Iteracion
%Fijamos Beta_h


Range(9,:) = T7_3.Nominal('beta_h');
[T,~]=gsua_dpmat(odes,vars,[0 150],'7m','output',1,'opt',opts,'Range',Range);
T.Properties.CustomProperties.output = 1;
M = gsua_dmatrix(T,100);
%%
%Realizamos analisis otra vez
Tsa = gsua_sa(M,T,'parallel', false, 'SensMethod', 'Xiao', 'ynom', ynom);

gsua_plot('Bar',Tsa,Tsa.STi)
savefig('iteration4/figures/SensibilityAnalisis.fig')


% Confiabilidad

c = sum(Tsa.Si)/sum(abs(Tsa.Si))

%%
%5.2 estimation
opt = optimoptions('lsqcurvefit','UseParallel',false,'Display','iter');
[T7_3,res_3] = gsua_pe(T,xdata,ydata2,'N',100,'opt',opt); 
save('iteration4/values/Results7.mat','T7_3','res_3','xdata','ydata2')

%%
%5.2 Identifiability analysis 
th_3 = sum(res_3<res_3(1)*1.02)
y7_3 = gsua_eval(T7_3.Estlsqc(:,1:th_3),T7_3,xdata,ydata2);
savefig('iteration4/figures/curves.fig')

%%
bestest_3 = y7_3(1,:);
trend = [bestest_3(1),bestest_3(2:end)-bestest_3(1:end-1)];
plot(trend,'b')
hold on
plot(diff(ydata2),'r')
title('Estimated vs Real Weekly Infections')
xlabel('Weeks')
ylabel('Cases')
legend({'Estimated','Real'})
savefig('iteration4/figures/EstimatedvsReal.fig')

%%
T7_3.Nominal = T7_3.Estlsqc(:,1);
%res : funciones de costo
th_3 = sum(res_3<res_3(1)*1.02) % threshold, cambiar el 1.01
%th = 476 %el parametro que habia dado el profe, a nosotros nos da 1????

T7_3 = gsua_ia(T7_3,T7_3.Estlsqc(:,1:th_3), false, true); 
savefig('iteration4/figures/Correlations.fig')

%%
%Uncertainty analysis

%analisis de incertidumbre con T7_3. veamos que las curvas hagan una banda 
%al lado de los datos que estamos estimando. Si la banda es chevere,
%terminamos, si no, fijar algun parametro. 

Ua = gsua_ua(M, T7_3, 'parallel', false, 'ynom',ydata2);
savefig('iteration4/figures/Montecarlo.fig')

%%
% Quinta Iteracion
%Fijamos Lambda


Range(9,:) = T7_3.Nominal('beta_h');
[T,~]=gsua_dpmat(odes,vars,[0 150],'7m','output',1,'opt',opts,'Range',Range);
T.Properties.CustomProperties.output = 1;
M = gsua_dmatrix(T,100);
%%
%Realizamos analisis otra vez
Tsa = gsua_sa(M,T,'parallel', false, 'SensMethod', 'Xiao', 'ynom', ynom);

gsua_plot('Bar',Tsa,Tsa.STi)
savefig('iteration4/figures/SensibilityAnalisis.fig')


% Confiabilidad

c = sum(Tsa.Si)/sum(abs(Tsa.Si))

%%
%5.2 estimation
opt = optimoptions('lsqcurvefit','UseParallel',false,'Display','iter');
[T7_3,res_3] = gsua_pe(T,xdata,ydata2,'N',100,'opt',opt); 
save('iteration4/values/Results7.mat','T7_3','res_3','xdata','ydata2')

%%
%5.2 Identifiability analysis 
th_3 = sum(res_3<res_3(1)*1.02)
y7_3 = gsua_eval(T7_3.Estlsqc(:,1:th_3),T7_3,xdata,ydata2);
savefig('iteration4/figures/curves.fig')

%%
bestest_3 = y7_3(1,:);
trend = [bestest_3(1),bestest_3(2:end)-bestest_3(1:end-1)];
plot(trend,'b')
hold on
plot(diff(ydata2),'r')
title('Estimated vs Real Weekly Infections')
xlabel('Weeks')
ylabel('Cases')
legend({'Estimated','Real'})
savefig('iteration4/figures/EstimatedvsReal.fig')

%%
T7_3.Nominal = T7_3.Estlsqc(:,1);
%res : funciones de costo
th_3 = sum(res_3<res_3(1)*1.02) % threshold, cambiar el 1.01
%th = 476 %el parametro que habia dado el profe, a nosotros nos da 1????

T7_3 = gsua_ia(T7_3,T7_3.Estlsqc(:,1:th_3), false, true); 
savefig('iteration4/figures/Correlations.fig')

%%
%Uncertainty analysis

%analisis de incertidumbre con T7_3. veamos que las curvas hagan una banda 
%al lado de los datos que estamos estimando. Si la banda es chevere,
%terminamos, si no, fijar algun parametro. 

Ua = gsua_ua(M, T7_3, 'parallel', false, 'ynom',ydata2);
savefig('iteration4/figures/Montecarlo.fig')

%%
T7_3.Properties.CustomProperties.output = 1:4;
gsua_eval(T7_3.Estlsqc(:,1:th_3),T7,xdata,ydata2)