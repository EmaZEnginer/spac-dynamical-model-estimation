% Definimos variables y parametros
% Definimos variables y parametros
clear
syms x_eq w_eq v_eq u_eq              %variables de estado y control
%syms kappa g m delta r beta k gamma rho
kappa=0.99373; g=0.10390; m=1.07359;
delta=1.02885; r=1.22763; beta=1.3547;
k=0.01003; gamma=0.00001; rho=0.01146;

%Definir las tres ecuaciones (Ecuación == 0)
%ecuaciones (x',w' y v')
f=[
    (kappa*x_eq*w_eq)/(1+g*x_eq)-m*x_eq;
    (delta*w_eq*v_eq)/(1+r*w_eq) - beta*w_eq - (k*x_eq*w_eq)/(1+g*x_eq);
    -gamma*v_eq - rho*(w_eq*v_eq)/(1+r*w_eq) + u_eq;
    ];

eqns=[f(1)==0, f(2)==0, f(3)==0];

sln=vpasolve(eqns,[x_eq w_eq v_eq]);
sln_m=[sln.x_eq sln.w_eq sln.v_eq];
digits(7);
sln_reducida=vpa(sln_m);

sln1=sln_reducida(1,:);   %solucion valida en funcion de u
sln1=subs(sln1,u_eq,2//.5);

%linealizacion de matriz A
x=[x_eq;w_eq;v_eq];
A_lin=jacobian(f,x);
A_lin_eval=subs(A_lin,{x_eq,w_eq,v_eq},sln1); %queda en funcion de U

C=[1 0 0]; %matriz de salida para los 3 estados

%observabilidad (matriz de observabilidad)
obs=[C;C*A_lin_eval;C*(A_lin_eval)^2];
rango=rank(obs) %se obtiene que el rango es completo

%valores para los cuales no es observable
det_obs=det(obs)
u_no_obs=vpasolve(det_obs==0,u_eq)




