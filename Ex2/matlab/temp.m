clc; clear;


syms Q1 Q2 Q3 gamma

rho = Q1;
u = Q2/Q1;
p = (gamma-1)*(Q3-1/2*Q2^2/Q1);
e = p/(gamma-1)+1/2*rho*u^2;

E1 = rho*u;
E2 = p+rho*u^2;
E3 = u*(e+p);

dE_dQ = [diff(E1,Q1), diff(E1,Q2), diff(E1,Q3);
         diff(E2,Q1), diff(E2,Q2), diff(E2,Q3);
         diff(E3,Q1), diff(E3,Q2), diff(E3,Q3)];

syms rho u a
Q1 = rho;
Q2 = rho*u;
e = (a^2*rho/gamma)/(gamma-1)+1/2*rho*u^2;
Q3 = e;

dE_dQ = simplify(subs(dE_dQ));

[T_R,Lambda_R] = eig(dE_dQ);
Lambda_R = simplify(Lambda_R)
T_R = simplify(T_R)
T_R*Lambda_R/T_R

%%
clc; clear;


syms Q1 Q2 Q3 gamma R A B C D u_x T_x mu kappa

rho = Q1;
u = Q2/Q1;
p = (gamma-1)*(Q3-1/2*Q2^2/Q1);
e = p/(gamma-1)+1/2*rho*u^2;
T = p/(rho*R);
% mu = A*T^(2/3)/(T+B);
% kappa = C*T^(2/3)/(T+D);

v1 = 0;
v2 = 4/3*mu*u_x;
v3 = 4/3*mu*u*u_x+kappa*T_x;

dv_dQ = [diff(v1,Q1), diff(v1,Q2), diff(v1,Q3);
         diff(v2,Q1), diff(v2,Q2), diff(v2,Q3);
         diff(v3,Q1), diff(v3,Q2), diff(v3,Q3)];

syms rho u a
Q1 = rho;
Q2 = rho*u;
e = (a^2*rho/gamma)/(gamma-1)+1/2*rho*u^2;
Q3 = e;

dv_dQ = simplify(subs(dv_dQ))
[T_P,Lambda_P] = eig(dv_dQ);
Lambda_P = simplify(Lambda_P)
T_P = simplify(T_P)
