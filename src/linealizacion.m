clear;
close all;
%activacion de graficas
graficas=true;        
%espacio para simular el modelo del arituclo 2
%parametros
kappa=0.99373; g=0.10390; m=1.07359;
delta=1.02885; r=1.22763; beta=1.3547;
k=0.01003; gamma=0.00001; rho=0.01146;

%% Linealizacion y Parametros Base
params=[kappa; g; m; delta; r; beta; k; gamma; rho];
u_eq=4.0;        %entrada de equilibrio (usada para linealizar)
perturbacion=2.0; %magnitud del cambio
tiempo_step=100;  %tiempo en que ocurre el cambio
tspan=0:0.5:1000; %vector de tiempo de simulación

%calculamos punto de equilibrio (usando u_eq)
p_eq_sym=puntos_eq(params, u_eq);
p_eq=double(p_eq_sym(1,:)); %usamos la primera solucion real positiva
x0=p_eq'; %condiciones iniciales en el equilibrio (columna)

%la matriz B es calculada en u_eq (0.08).
[A_lin, B_lin] = funcion_linealizacion(params, p_eq, u_eq);
valores_propios=eig(A_lin);
tiempo_establecimiento=4/max(valores_propios);
disp(['Matriz lineal A con entrada u=' num2str(u_eq)])
disp(A_lin)
disp('Valores propios: ')
disp(valores_propios)

%% MODELO NO LINEAL CON ENTRADA VARIABLE
% x:tamano de masa en la planta, w:agua dentro de la planta, v:agua de
% entrada (riego)

%modificamos la función f para que reciba el tiempo y calcule u(t)
f=@(t,x)[
    kappa*x(1)*x(2)/(1+g*x(1)) - (m*x(1));
    (delta*x(2)*x(3))/(1+r*x(2)) - (beta*x(2)) - (k*x(1)*x(2))/(1+g*x(1));
    (-gamma*x(3)) - (rho*x(2)*x(3)/(1+r*x(2))) + (u_eq+perturbacion*(t>=tiempo_step)) % u(t)
    ];

[t,x]=ode45(f, tspan, x0);

%% MODELO LINEAL CON ENTRADA VARIABLE
%desviacion inicial (es cero si empezamos en el equilibrio)
dx0=zeros(3, 1); 

%definimos la función f_lin para que calcule la desviaciòn de la entrada du(t)
f_lin=@(t,dx) A_lin*dx + B_lin*(perturbacion*(t>=tiempo_step));
%nota: du(t)=u(t)-u_eq. Como u(t)=u_eq+perturbacion*(...), entonces du(t)=perturbacion*(...)

[t,dx]=ode45(f_lin, tspan, dx0);

%recuperamos el estado lineal (sumando el equilibrio inicial)
x_lin=dx+p_eq; 

%% Graficación y Comparación
if graficas
    figure(1)
    plot(t, x);
    hold on
    plot(t, x_lin, '--', 'LineWidth', 1.5);
    title(['Sistema lineal (--) vs no lineal (-) con perturbación de u: ', ...
        num2str(u_eq), ' \rightarrow ', num2str(u_eq + perturbacion)]);
    xlabel('Tiempo (días)');
    ylabel('Estados');
    legend('x (NL)', 'w (NL)', 'v (NL)', 'x (L, --)', 'w (L, --)', 'v (L, --)','Perturbacion');
    xline(tiempo_step, 'r--', 'Perturbacion');
    grid on;
end

%calculamos la diferencia entre el lineal y el no lineal
size(x_lin);
size(x);
error=abs(x-x_lin);

if graficas
    figure(2)
    plot(t, x_lin)
    legend('Δx_1: Biomasa','Δx_2: Agua en planta','Δx_3: Agua en suelo')
    xlabel('Tiempo (días)')
    ylabel('Desviación respecto al equilibrio')
    title('Modelo lineal (con C.I equilibrio)')
    grid on
end

if graficas
    figure(3)
    plot(t,error)
    title('Diferencia entre sistema linealizado y no linealizado (con C.I equilibrio)')
    grid on
    legend('Error variable x', 'Error variable w', 'Error variable v')

    figure(4)
    grid on
    hold on
    plot(t,x)
    plot(t,x_lin,'--')
    title('Lineal (--) vs no lineal (-) (con C.I equilibrio)')
end

%porcentaje de error respecto al sistema
max_lin=max(x_lin);
max_err=max(error);
relacion_err=(max_err./max_lin)*100;

%%

%funciones del arhcivo
%funcion de linealizacion
function [A, B]=funcion_linealizacion(params,p_eq,u_lin)
    %extraemos los parametros 
    kappa=params(1); g=params(2); m=params(3);
    delta=params(4); r=params(5); beta=params(6);
    k=params(7); gamma=params(8); rho=params(9);

    %definimos la funcion con variables simbolicas
    syms x1 x2 x3 u

    x=[x1; x2; x3];

    f=[
    kappa*x1*x2/(1+g*x1) - (m*x1);
    (delta*x2*x3)/(1+r*x2) - (beta*x2) - (k*x1*x2)/(1+g*x1);
    (-gamma*x3) - (rho*x2*x3/(1+r*x2)) + u
    ];

    %obtenemos las matrices lineales
    A_gen=jacobian(f,x);
    B_gen=jacobian(f,u);

    A=double(subs(A_gen,{x1,x2,x3,u},{p_eq(1),p_eq(2),p_eq(3),u_lin}));
    B=double(subs(B_gen,{x1,x2,x3,u},{p_eq(1),p_eq(2),p_eq(3),u_lin}));
end

%funcion: calculo punto eq
function p_eq=puntos_eq(params,u)
    %extraemos los parametros del vector
    kappa=params(1); g=params(2); m=params(3);
    delta=params(4); r=params(5); beta=params(6);
    k=params(7); gamma=params(8); rho=params(9);
    
    syms x1 x2 x3
    
    %ponemos ponemos el modelo
    f_eq=[
    kappa*x1*x2/(1+g*x1) - (m*x1);
    (delta*x2*x3)/(1+r*x2) - (beta*x2) - (k*x1*x2)/(1+g*x1);
    (-gamma*x3) - (rho*x2*x3/(1+r*x2)) + u
    ];

    %resolvemos las ecuaciones
    eqns=[f_eq(1)==0,f_eq(2)==0,f_eq(3)==0];
    sln=vpasolve(eqns,[x1 x2 x3]);
    sln_m=[sln.x1 sln.x2 sln.x3];
    
    %seleccionamos los equilibrios >=0
    filas_deseadas=all(sln_m>0 & isreal(sln_m),2);
    p_eq=sln_m(filas_deseadas,:);
end

