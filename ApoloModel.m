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

[T, M, ynom,Tsa] = Sensitivity(odes,vars,Range,opts,ydata)

%%
%4. Fixing parameters

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
%5. Parameter estimation



opt = optimoptions('lsqcurvefit','UseParallel',false,'Display','iter');
[T7,res] = gsua_pe(T,xdata,ydata2,'N',100,'opt',opt); %mover N ( 100 buena, 50 aceptable, 10 :( )
save('iteration1/values/Results7.mat','T7','res','xdata','ydata2')