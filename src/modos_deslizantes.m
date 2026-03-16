clear;
close all;
clc;

%% definificion de parametros
%parametros del cultivo
kappa=0.99373; g=0.10390; m=1.07359;
delta=1.02885; r=1.22763; beta=1.3547;
k_param=0.01003; gamma=0.00001; rho=0.01146;
params=[kappa; g; m; delta; r; beta; k_param; gamma; rho];

%condiciones de operacion
u_eq = 0.08;          %entrada de eq
perturbacion = 0.05;  %perturbacion en la anrtada
tiempo_step = 19;     %momento de perturbacion
T_sim = 20;          %tiempo de simualcion
dt = 0.01;            %paso de integracion
ruido_std = 0.02;     %ruiedo de medicion

%puntos de eq y linealizacion
p_eq=puntos_eq(params, u_eq);
disp('Punto de equilibrio (x1, x2, x3):');
disp(p_eq);
[A, B]=funcion_linealizacion(params, p_eq, u_eq);

%% diseño del observador

%deficion de la salida medida
C = [1 0 0]; 

% El sistema se queda
% dot_y = A11*y + A12*x_rest + ...
% dot_x_rest = A21*y + A22*x_rest + ...

%medimos x1 (indice 1), la particion es:
A11 = A(1,1);
A12 = A(1, 2:3);      %enlace entre medidos y no medidos
A21 = A(2:3, 1);
A22 = A(2:3, 2:3);    %dinamica de los no medidos

%calculo de L2 (dinamica de los no medidos)
%(A22 - L2*A12) sea estable.
%usamos place
polos_deseados_restantes = [-1.5; -2.0]; % Polos rápidos y estables
L2_transpuesta = place(A22', A12', polos_deseados_restantes);
L2=L2_transpuesta'; 

disp('Matriz L2 calculada:');
disp(L2);

%seleccion de L1 (ganacia de conmutacion)
%mayor que cota de perturbacion y ruido para garantizar que deslice
L1=5; %l1 grande, mas robustex y mayor convergencia mucho chattering
      %le pequeño menos robustez y convergencia pero menos chattering

%matrices k
%K_obs = [L1; L2*L1]
K_sm = [L1; L2*L1]; 

disp('ganancia total (K_sm):');
disp(K_sm);

%% simualcion
time = 0:dt:T_sim;
steps = length(time);

%inicializacion de vectores
X_real = zeros(3, steps);
X_real(:,1) = p_eq'; %planta arranca en eq

%vector para observador
% z_hat es la desviacion estimada de x_hat total = z_hat + p_eq
z_hat = [3; -2; 3]; %condiciones iniciales del observador

X_est=zeros(3, steps);
X_est(:,1)=z_hat + p_eq';

%generacion de vectores para guardar los datos
hist_sign = zeros(1, steps); %para ver la señal de conmutacion

for k = 1:steps-1
    % planta no lineal
    u_t = u_eq;
    if time(k) >= tiempo_step
        u_t = u_eq + perturbacion; %perturbacion (desconocida por el obs)
    end
    
    x_k = X_real(:, k);
    dx = modelo_nolineal(x_k, u_t, params);
    X_real(:, k+1) = x_k + dt * dx;
    
    %medicion
    %salida con ruido
    y_med = C * X_real(:, k) + (0*ruido_std * randn);
    
    %desviacion medida (lo que va al observador)
    y_eq_val = C * p_eq';
    dy_med = y_med - y_eq_val; 
    
    % OBSERVADOR
    % dot_z_hat = A*z_hat + B*du + K_sm * sign(dy_med - C*z_hat)
    
    %salida estimada
    dy_est = C * z_hat;
    
    %error de salida y signo del error
    error_y = dy_med - dy_est;
    termino_signo = sign(error_y);
    
    %se guarda el dato para graficar el chattering
    hist_sign(k) = termino_signo;
    
    du=u_t-u_eq;

    %dinamica del observador
    dz_hat = A * z_hat + B * du + K_sm * termino_signo;
    
    %integracion por euler
    z_hat = z_hat + dt * dz_hat;
    
    %recuperacion de estado original
    X_est(:, k+1) = z_hat + p_eq';
end

%% graficas
nombres = {'Biomasa (x1)', 'Agua Planta (x2)', 'Agua Suelo (x3)'};

for i = 1:3
    figure('Color','w'); hold on; grid on;
    plot(time, X_real(i,:), 'b', 'LineWidth', 2, 'DisplayName', 'Real (NL)');
    plot(time, X_est(i,:), 'r--', 'LineWidth', 1, 'DisplayName', 'Estimado (OMD)');
    
    if i == 1 %solo mostramos medición en x1
        plot(time, X_real(i,:) + ruido_std*randn(1,steps), 'g.', 'MarkerSize',1, 'DisplayName','Medición ruidosa');
    end
    
    xline(tiempo_step, 'k:', 'DisplayName', 'Perturbación');
    title(['Estado ' nombres{i} ' - Observador Deslizante']);
    xlabel('Tiempo (días)'); ylabel(nombres{i});
    legend('show');
end

%grafica del conmutacion
figure('Color','w'); 
plot(time, hist_sign, 'k');
title('Señal de Conmutación (sign(e_y))');
xlabel('Tiempo (días)');
ylim([-1.5 1.5]);
grid on;

%calculo de errores
error_est = X_real - X_est;
e_x1 = error_est(1, :); %error biomasa
e_x2 = error_est(2, :); %error agua planta

%DIAGRAMA DE FASE DE LOS ERRORES (e_x2 vs e_x1)
%como viaja el error al origen por la superficie deslizante
figure('Color','w');
plot(e_x1, e_x2, 'b', 'LineWidth', 1);
hold on;

%se marca el origen
plot(0, 0, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
title('Diagrama de Fase del Error de Estimación');
xlabel('Error variable medida (e_{x1}) --> Superficie s=0');
ylabel('Error variable NO medida (e_{x2})');
grid on;

%supf deslizante una linea vertical donde e_x1=0 (pq es la que%medimos)
xline(0, 'r--', 'LineWidth', 1, 'DisplayName', 'Superficie Deslizante (e_{x1}=0)');
legend('Trayectoria del Error', 'Origen (Objetivo)', 'Superficie Deslizante');


%% funciones auxiliares
function dx = modelo_nolineal(x, u, params)
    kappa=params(1); g=params(2); m=params(3);
    delta=params(4); r=params(5); beta=params(6);
    k_param=params(7); gamma=params(8); rho=params(9);
    x1 = x(1); x2 = x(2); x3 = x(3);
    dx = zeros(3,1);
    dx(1) = kappa*x1*x2/(1+g*x1) - m*x1;
    dx(2) = (delta*x2*x3)/(1+r*x2) - beta*x2 - (k_param*x1*x2)/(1+g*x1);
    dx(3) = -gamma*x3 - (rho*x2*x3)/(1+r*x2) + u;
end

function [A, B] = funcion_linealizacion(params, p_eq, u_lin)
    kappa=params(1); g=params(2); m=params(3);
    delta=params(4); r=params(5); beta=params(6);
    k_param=params(7); gamma=params(8); rho=params(9);
    syms x1 x2 x3 u
    f = [
        kappa*x1*x2/(1+g*x1) - (m*x1);
        (delta*x2*x3)/(1+r*x2) - (beta*x2) - (k_param*x1*x2)/(1+g*x1);
        (-gamma*x3) - (rho*x2*x3/(1+r*x2)) + u
    ];
    A = double(subs(jacobian(f, [x1; x2; x3]), {x1,x2,x3,u}, {p_eq(1),p_eq(2),p_eq(3),u_lin}));
    B = double(subs(jacobian(f, u), {x1,x2,x3,u}, {p_eq(1),p_eq(2),p_eq(3),u_lin}));
end

function p_eq = puntos_eq(params, u_val)
    kappa=params(1); g=params(2); m=params(3);
    delta=params(4); r=params(5); beta=params(6);
    k_param=params(7); gamma=params(8); rho=params(9);
    syms x1 x2 x3
    f_eq = [
        kappa*x1*x2/(1+g*x1) - (m*x1);
        (delta*x2*x3)/(1+r*x2) - (beta*x2) - (k_param*x1*x2)/(1+g*x1);
        (-gamma*x3) - (rho*x2*x3/(1+r*x2)) + u_val
    ];
    sol=vpasolve([f_eq(1)==0, f_eq(2)==0, f_eq(3)==0], [x1, x2, x3]);
    sols = [double(sol.x1), double(sol.x2), double(sol.x3)];
    valid = all(sols > 0 & isreal(sols), 2);
    if isempty(valid) || sum(valid) == 0
        p_eq = [0 0 0];
    else
        p_eq = sols(valid, :);
        p_eq = p_eq(1, :); 
    end
end