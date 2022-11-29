A = [2 0; -3 3];
B = [1;3];
C = [9 -3];
D = 0;
ss(A,B,C,D)

observabilidad =[C;C*A]

cond(observabilidad)

rank(observabilidad)

inv(observabilidad)


%%
A = [2 0; -3 3];
B = [1;3];
C = [9 -3];
D = 0;
ss(A,B,C,D)

observabilidad =[C;C*A]

cond(observabilidad)

rank(observabilidad)

inv(observabilidad)

%%

A = [0 1; 3 0];
B = [1; 3];
C = [1 0];
D = 0;

observabilidad =[C;C*A]

rank(observabilidad)
cond(observabilidad)
inv(observabilidad)

%%

A = [0 1; 3 0];
B = [1; 3];
C = [0 1];
D = 0;

observabilidad =[C;C*A]

rank(observabilidad)
cond(observabilidad)
inv(observabilidad)
%%

A = [2 0; 0 1];
B = [0; 1];
C = [1 0];
D = 0;

observabilidad =[C;C*A]

rank(observabilidad)
cond(observabilidad)
inv(observabilidad)
%%

A = [2 0; 0 1];
B = [0; 1];
C = [1 1];
D = 0;

observabilidad =[C;C*A]

rank(observabilidad)
cond(observabilidad)
inv(observabilidad)
%%

A = [2 0; 0.02 1];
B = [0; 1];
C = [1 0];
D = 0;

observabilidad =[C;C*A]

rank(observabilidad)
cond(observabilidad)
inv(observabilidad)
%%
syms s
den = s^2 + 3*s + 2;
factor(den);
G = zpk([-3],[-1 -2],1)

ss(G)

%%



G = zpk([-2 -3 -4],[-1 -2 -3 -4],1)

ss(G)

%%
syms x1 x2
Ahat = [0 1; -(1+0.2*x1*x2) 0.1*(1-x1^2)] == 0;

f1 = -(1+0.2*x1*x2) ==0;
As = solve(Ahat)
