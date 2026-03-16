clear;
close all;
clc;

%% definicion de params
kappa=0.99373; g=0.10390; m=1.07359;
delta=1.02885; r=1.22763; beta=1.3547;
k_param=0.01003; gamma=0.00001; rho=0.01146;
params=[kappa; g; m; delta; r; beta; k_param; gamma; rho];

%puntos de eq y linealizacion
u_eq=2.5; 
p_eq=puntos_eq(params, u_eq);
p_eq=p_eq';
[A_cont, B_cont]=funcion_linealizacion(params, p_eq, u_eq); 

%matrices de transcion Fk y Gk
dt=0.05; 
Fk=eye(3)+A_cont*dt;
Gk=B_cont*dt; %aproximadas por euler

%% parametros y diseño del filtro
T_sim=30;        %tiempo y stpes
time=0:dt:T_sim;
N_steps=length(time);

%matrices de Covarianza (Diseño Ingenieril)
P=diag([10, 100, 1]); 
factor_ruido=4;
Q=factor_ruido.*diag([0.01, 0.01, 0.05]); 

%variable medida (biomasa)
C=[1 0 0]; 

ruido_std=0.05;
R=factor_ruido*diag(ruido_std^2); 

%ENTRADA DE CONTROL (u) y PERTURBACIONES
u_history=u_eq*ones(1, N_steps);
perturbacion=1;
tiempo_step=19; 
u_history(time >= tiempo_step)=u_eq+perturbacion;

%incializacion de vectores
X_real=zeros(3, N_steps);       %real   
X_est_dev=zeros(3, N_steps);    %estimado
Y_med=zeros(size(C, 1), N_steps); %medido

%condiciones iniciales
X_real(:, 1)=p_eq;                  
X_est_dev(:, 1)=[200; -200; 180];  %estimado

%generacion de ruido para el modelo real
Noise_w=factor_ruido*(chol(Q)*randn(3, N_steps))'; 
Noise_v=factor_ruido*(chol(R)*randn(size(C, 1), N_steps))'; 

%% bucle del filtro de kalman
dx_hat=X_est_dev(:, 1); %inicializacion de primeros valores
P_k=P;                  
for k = 1:N_steps - 1
    u_prev = u_history(k);
    
    %SIMULACION DE LA PLANTA REAL (NO LINEAL)
    x_k=X_real(:, k);
    dx_real=modelo_nolineal(x_k, u_prev, params);
    X_real(:, k+1)=x_k+dt*dx_real+Noise_w(k, :)'; 
    %resuelta con el metodo de euler
    
    %MEDICION (en desviacion) 
    y_k_abs=C*X_real(:, k+1)+Noise_v(k, :)'; 
    y_eq_val=C*p_eq;       %en los puntos de equlibrio
    dy_k=y_k_abs-y_eq_val; % dy_k es es desviado
    Y_med(:, k+1) = y_k_abs;
    
    %PREDICCION (A PRIORI) 
    du_k=u_prev-u_eq;               %u(k-1) en desviacion
    dx_hat_minus=Fk*dx_hat+Gk*du_k; 
    P_minus=Fk*P_k*Fk'+Q; 
    
    %ACTUALIZACIonN (A POSTERIORI) 
    S=C*P_minus*C'+R; %calculo de las ganancias K
    Kk=(P_minus*C')*inv(S);    %
    
    Rk=R; %Rk constante porque incertidumbre de sensor no varia
    dy_hat_minus=C*dx_hat_minus; %mi y estimada
    dx_hat=dx_hat_minus + Kk * (dy_k - dy_hat_minus); %correcion de dx
    P_k=(eye(3)-Kk*C)*P_minus+(Kk*Rk*Kk'); %correcion de P
    
    %Almacenamiento de la estiamcion
    X_est_dev(:, k+1) = dx_hat;
end

%reconstruccion de estado absoluto
X_est_abs=X_est_dev+p_eq;

%% graficas
nombres_estados={'x_1 (Biomasa)', 'x_2 (Agua en Planta)', 'x_3 (Agua en Suelo)'};

%estimacion independiente para cada kalman
titulos_figuras = {'Figura 1: Estimación de X_1', 'Figura 2: Estimación de X_2', 'Figura 3: Estimación de X_3'};

%bucle para generar figuras INDEPENDIENTES para cada estado
for i = 1:3
    figure('Color','w');
    %se muestra tanto la realidad (NL) como la estimación (LKF)
    title([titulos_figuras{i} ' vs. Planta NL (Solo X1 Medida)']); 
    hold on; grid on;
    
    %planta Real (con ruido)
    plot(time, X_real(i, :), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Real (NL) con Ruido');
    
    %estimacion de Kalman
    plot(time, X_est_abs(i, :), 'r', 'LineWidth', 1, 'DisplayName', 'Estimado (LKF)');
    
    %linea de perturbacion
    xline(tiempo_step, 'k:', 'LineWidth', 1.5, 'HandleVisibility','off'); 
    
    %grafica de medicion (Solo si es X1)
    if i == 1 
        % Solo X1 tiene medicion (no es necesario estimar)
        plot(time, Y_med(1, :), 'g.', 'MarkerSize', 5, 'DisplayName', 'Medición Ruidosa x1');
    end
    
    ylabel(nombres_estados{i});
    xlabel('Tiempo (días)');
    legend('show', 'Location', 'best');
end

figure(4)
for i = 1:3
    subplot(3, 1, i); hold on; grid on;
    plot(time, X_real(i, :), 'b', 'LineWidth', 1.5, 'DisplayName', 'Real (NL) con Ruido');
    plot(time, X_est_abs(i, :), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Estimado (LKF)');
    
    if i == 1 
        % Solo hay medición para X1
        plot(time, Y_med(1, :), 'g.', 'MarkerSize', 2, 'DisplayName', 'Medición Ruidosa x1');
        xline(tiempo_step, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Perturbación');
        legend('show', 'Location', 'best');
    else 
        % Para X2 y X3, solo la estimación y el valor real.
        xline(tiempo_step, 'k:', 'LineWidth', 1.5, 'HandleVisibility','off'); 
        legend('show', 'Location', 'best');
    end
    ylabel(nombres_estados{i});
end
xlabel('Tiempo (días)');

%error de estimacion
nombres_estados = {'x_1 (Biomasa)', 'x_2 (Agua en Planta)', 'x_3 (Agua en Suelo)'};
error_abs = abs(X_real - X_est_abs); 

%bucle para las figuras
for k = 1:3
    figure('Color','w'); %nueva figura por cada grafica
    
    %k-esima fila del k-esimo estado
    plot(time, error_abs(k, :), 'LineWidth', 1.5, 'DisplayName', ['Error |x' num2str(k) '|']);
    
    %titulos y etiquetas
    titulo_error = ['FIGURA ' num2str(3+k) ': Error Absoluto de Estimación de ' nombres_estados{k}];
    title(titulo_error); 
    ylabel(['Error |x' num2str(k) ' - \hat{x}' num2str(k) '|']);
    xlabel('Tiempo (días)');
    
    legend('show', 'Location', 'best'); 
    grid on;
end

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
    p_eq = sols(valid, :);
    if ~isempty(p_eq), p_eq = p_eq(1, :); else, p_eq = [0 0 0]; end
end