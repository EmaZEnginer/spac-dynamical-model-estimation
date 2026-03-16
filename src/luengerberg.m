%observador luengerberg para cultivo
clear; close all; clc;

%parametros
kappa=0.99373; g=0.10390; m=1.07359;
delta=1.02885; r=1.22763; beta=1.3547;
k=0.01003; gamma=0.00001; rho=0.01146;

params=[kappa; g; m; delta; r; beta; k; gamma; rho];
u_eq=2.5; %entrada constante

%punto de equilibrio y linealizacion
%hallamos puntos de operacion
p_eq=puntos_eq(params, u_eq);
disp('Punto de equilibrio (x1, x2, x3):');
disp(p_eq);

%matrices A y B del sismtea linealizado
[A, B]=funcion_linealizacion(params, p_eq, u_eq);
%C se supone que mide biomasa
C=[1 0 0]; 

%dinamica de la planta en lazo abierto
polos_planta=eig(A);
disp('Polos de la planta (Lazo abierto):');
disp(polos_planta);

%tomamos el polo mas alejado de la planta y lo multiplicamos 3, 4 y 5 
%y asi sacamos la dinamica regida por el mas cerca
%ara simplificar y evitar errores numéricos, usaremos polos reales rápidos:
min_polo = min(real(polos_planta));
polos_deseados=[min_polo*3; min_polo*4; min_polo*5]; 
ts_establecimiento=4/max(polos_deseados);

disp('Polos deseados del observador:');
disp(polos_deseados);


%calculo de ganacia k
Ke=place(A', C', polos_deseados)';
disp('Matriz de ganancia Ke:');
disp(Ke);

%simulaciones
T_sim=500;           %tiempo de simulacion
dt=0.05;             %paso de integracion
time=0:dt:T_sim;
N=length(time);

%inicializacion de vectores
X_real=zeros(3, N);       %para la planta real
X_est_lin=zeros(3, N);    %para el estimado
Y_med=zeros(1, N);        %salidad medida con ruida

%condiciones iniciales (planta real) 
X_real(:, 1)=p_eq';          %en_eq
%X_real(:,1)=[0.2 1 0.5]';   %en otro punto

%el observador trabaja con las desviaciones 
%x_hat_deviacion=[0; 0; 0]; %empezando en 0
x_hat_desviacion=[2 -1 3]';

%ruido
ruido_std = 0.05; %desviacion del sensor
t_perturbacion=19;
magnitud_perturbacion=1;

%bucle de simualcion
for k=1:N-1
    % ---PLANTA REAL---
    x_k=X_real(:, k);
    
    %entrada de riego constante con perturbacion en t=t_preturbacion
    u_t=u_eq; 
    if time(k) > t_perturbacion, u_t = u_eq + magnitud_perturbacion; end
    
    %ecuaciones no lineales
    dx=modelo_nolineal(x_k, u_t, params);
    %integracion con euler
    X_real(:, k+1) = x_k + dt * dx;
    
    %medicion
    % Salida y=C*x+ Ruido
    y_k=C*X_real(:, k)+ruido_std * randn;
    Y_med(k)=y_k;
    
    % --- OBSERVADOR DE LUENBERGER (LINEAL) ---
    %Este estima respecto a la desviacion lineal (delta)
    %d(dx_hat)=A*dx_hat + B*du + Ke*(dy_medido - C*dx_hat)
    
    %entradas y salidas en desviación
    du=u_t-u_eq;              %desviacion de la entrada
    y_eq=C*p_eq';         %salida en el equilibrio
    dy_med=y_k-y_eq;      %desviacion de la salida
    
    %ecuacion del observador 
    dx_hat_dot = A*x_hat_desviacion + B * du + Ke * (dy_med - C * x_hat_desviacion);
    
    %solucion por euler
    x_hat_desviacion = x_hat_desviacion + dt * dx_hat_dot;
    
    %reconstruccion del estado: x_estimado=dx_estimado+p_eq
    X_est_lin(:, k+1) = p_eq' + x_hat_desviacion;
end

% %GRAFICAS
% nombres_estados = {'Biomasa (x1)', 'Agua Planta (x2)', 'Agua Suelo (x3)'};
% titulos_fig = {'Figura 1: Biomasa (x1) - Observador vs. Real',... 
%     'Figura 2: Agua Planta (x2) - Observador vs. Real',... 
%     'Figura 3: Agua Suelo (x3) - Observador vs. Real'};
% 
% for i = 1:3
%     figure('Color','w'); %nueva figura para cada grafica
%     hold on; grid on;
% 
%     %planta real conruido
%     plot(time, X_real(i, :), 'b', 'LineWidth', 1.5, 'DisplayName', 'Real (NL)');
% 
%     %estimacion del observador
%     plot(time, X_est_lin(i, :), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Estimado (Obs)');
% 
%     %medicion ruidosa
%     if i == 1 
%         %Y_med solo mide X1,solo en la primera grafica
%         plot(time, Y_med, 'g.', 'MarkerSize', 2, 'DisplayName', 'Medición Ruidosa');
%         % Ajustar la leyenda cuando la medición está presente
%         legend('show', 'Location', 'best');
%     else
%         %asjutamos leyenda cuando no está presente
%         legend('show', 'Location', 'best');
%     end
% 
%     %titulos y etiquetas
%     title(titulos_fig{i});
%     ylabel(nombres_estados{i});
%     xlabel('Tiempo (días)'); 
% end
% 
% %graficas de error
% %error de estimacion
% error_obs = X_real - X_est_lin; 
% 
% %bucle de graficas
% for k = 1:3
%     figure('Color','w'); %nueva figura por cada grafica
% 
%     %grafica k-esima fila del k-esimo estado
%     plot(time, error_obs(k, :), 'LineWidth', 1.5, 'DisplayName', ['Error ' nombres_estados{k}]);
% 
%     %titulos y etiquetas
%     titulo_error = ['FIGURA ' num2str(3+k) ': Error de Estimación (e = x - \hat{x}) de ' nombres_estados{k}];
%     title(titulo_error); 
%     ylabel(['Error e' num2str(k)]);
%     xlabel('Tiempo (días)');
% 
%     legend('show', 'Location', 'best'); 
%     grid on;
% end

%NOTA:ACA NOS ESTABA DANDO UN ERROR EN LA GRAFICACION PORQUE APARECIA
%DESPLAZADA LA GRAFICA, POR LO QUE IMPLEMENTAMOS LA SOLUCION DE QUITAR EL
%PRIMER PUNTO PARA CORREGIR ESA GRAFICACION Y EVITAR QUE ESE PRIMER 
%PUNTO ASOCIADO A LA GRAFICA SE QUEDARA EN DONDE NO DEBIA

% --- Recortar los dos primeros valores ---
time2     = time(2:end);
X_real2   = X_real(:, 2:end);
X_est_lin2 = X_est_lin(:, 2:end);
Y_med2    = Y_med(2:end);     % solo aplica a x1
error_obs2 = X_real2 - X_est_lin2;


%GRAFICAS
nombres_estados = {'Biomasa (x1)', 'Agua Planta (x2)', 'Agua Suelo (x3)'};
titulos_fig = {'Figura 1: Biomasa (x1) - Observador vs. Real',... 
    'Figura 2: Agua Planta (x2) - Observador vs. Real',... 
    'Figura 3: Agua Suelo (x3) - Observador vs. Real'};

for i = 1:3
    figure('Color','w'); %nueva figura para cada grafica
    hold on; grid on;
    
    %planta real conruido
    plot(time2, X_real2(i, :), 'b', 'LineWidth', 1.5, 'DisplayName', 'Real (NL)');
    %estimacion del observador
    plot(time2, X_est_lin2(i, :), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Estimado (Obs)');

    %medicion ruidosa
    if i == 1 
        %Y_med solo mide X1,solo en la primera grafica
        plot(time2, Y_med2, 'g.', 'MarkerSize', 2, 'DisplayName', 'Medición Ruidosa');
        % Ajustar la leyenda cuando la medición está presente
        legend('show', 'Location', 'best');
    else
        %asjutamos leyenda cuando no está presente
        legend('show', 'Location', 'best');
    end
    
    %titulos y etiquetas
    title(titulos_fig{i});
    ylabel(nombres_estados{i});
    xlabel('Tiempo (días)'); 
end

%bucle de graficas
for k = 1:3
    figure('Color','w'); %nueva figura por cada grafica
    
    %grafica k-esima fila del k-esimo estado
    plot(time2, error_obs2(k, :), 'LineWidth', 1.5, 'DisplayName', ['Error ' nombres_estados{k}]);
    
    %titulos y etiquetas
    titulo_error = ['FIGURA ' num2str(3+k) ': Error de Estimación (e = x - \hat{x}) de ' nombres_estados{k}];
    title(titulo_error); 
    ylabel(['Error e' num2str(k)]);
    xlabel('Tiempo (días)');
    
    legend('show', 'Location', 'best'); 
    grid on;
end


%FUNCIONES AUXILIARES
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
    
    A_sym = jacobian(f, [x1; x2; x3]);
    B_sym = jacobian(f, u);
    
    A = double(subs(A_sym, {x1,x2,x3,u}, {p_eq(1),p_eq(2),p_eq(3),u_lin}));
    B = double(subs(B_sym, {x1,x2,x3,u}, {p_eq(1),p_eq(2),p_eq(3),u_lin}));
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
    
    %filtra soluciones reales y positivas
    sols = [double(sol.x1), double(sol.x2), double(sol.x3)];
    valid = all(sols > 0 & isreal(sols), 2);
    p_eq = sols(valid, :);
    if isempty(p_eq)
        p_eq = [0 0 0]; % Fallback
    else
        p_eq = p_eq(1, :); %toma la primera valida
    end
end