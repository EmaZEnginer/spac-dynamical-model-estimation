clear;
close all;
clc;

%parametros
kappa=0.99373; g=0.10390; m=1.07359;
delta=1.02885; r=1.22763; beta=1.3547;
k_param=0.01003; gamma=0.00001; rho=0.01146;
params=[kappa; g; m; delta; r; beta; k_param; gamma; rho];

%condiciones de operacion
u_eq=2.5;            %entrada de eq
perturbacion=1;    % Magnitud de la entrada desconocida
tiempo_step=1.9;        % Momento donde inicia la perturbación desconocida
T_sim=2;            % Tiempo de simulación
dt=0.05;                % Paso de integración
ruido_std=0.5;          %desviacion del ruido    

% ---LINEALIZACION---
%calcular p_eq y linealizar
p_eq = puntos_eq(params, u_eq);
disp('Punto de equilibrio (x1, x2, x3):');
disp(p_eq);
[A, B]=funcion_linealizacion(params, p_eq, u_eq);

% --- OBSERVADOR DE ENTRADAS DESCONOCIDAS (OEDL) ---

%como mi modelo cambio y es X'=Ax+Bu+Dv pero u=0, y v=u entonces D=B
D=B; 

%condicion de enfrentamiento: C*D=dim(V)
C=[0 0 1]; 
disp('Matriz C utilizada (debe cumplir rank(C*D) != 0):');
disp(C);
rango_cd=rank(C*D);
rango=rank(D);

M_CD=C*D;

%diseño
%i) calcular E
E=-D*inv(M_CD);

%ii) calcular P
P=eye(3)+E*C;

%iii)observabilidad del par (PA, C) para ubicar polos
if rank(obsv(P*A, C)) < 3
    warning('El sistema equivalente (PA, C) no es completamente observable.');
end

% --- DISEÑO DEL OBSERVADOR DE ENTRADAS DESCONOCIDAS (OEDL) ---

%1.calcular los polos de la matriz de dinámica transformada P*A.
polos_PA = eig(P*A);

%como solo se puede cambiar un polo, dejamos dos estables de la
%planta en polos deseados y movemos el otro, la funcion place se
%encarga de calcular las k de tal forma que deje los polos fijos
%que no se pueden modificar de la planta y el que le pusimos
polos_fijos=polos_PA; 

%ordenamos por magnitud (solo por seguridad)
[~, idx] = sort(real(polos_fijos));
polos_fijos_ordenados=polos_fijos(idx);

%definimos los polos deseados: Mantenemos los dos fijos y movemos uno.
polos_deseados=[
    polos_fijos_ordenados(1); %polo fijo 1 (de la planta)
    polos_fijos_ordenados(2); %polo fijo 2 (de la planta)
    -1.5                     %polo que sí podemos mover lo ponemos
];                          %alejado para que sea rapido

ts_establecimiento=4/max(polos_deseados);

disp('Polos fijos de PA (modos no modificables):');
disp(polos_fijos_ordenados(1:2));
disp('Polo deseado para el modo observable:');
disp(polos_deseados(3));

%2)calculamos k con dos polos fijos y 1 deseado
K_gain=place((P*A)', C', polos_deseados)';

disp('Matriz de ganancia K_gain (ajustada):');
disp(K_gain);

%v)calcular matrices finales N, L, G
N=P*A-K_gain*C;
L=K_gain*(1 + C*E)-P*A*E;
G=P*B; % matriz para la parte CONOCIDA de la entrada (si la hubiera)

disp('Matriz N del observador:'); disp(N);

% --- SIMULACION ---
time=0:dt:T_sim;
steps=length(time);

%inicializar para modelo real
X_real=zeros(3, steps);
X_real(:,1)=p_eq'; %en eq

%el observador estima desviaciones (z), luego recuperamos x_hat
%inicializamos z cerca del cero (desviación nula)
z=zeros(3,1); 
X_est=zeros(3, steps);
%X_est(:,1)=p_eq'; %estimado inicial (asumimos equilibrio) NOTA:Aca
                    %no usamos la desviacion por comodidad
X_est(:,1)=[p_eq(1)-21,p_eq(2)-15,p_eq(3)+25]'; %estimado inicial en otro punto de desviación
%ruido


%bucle de simulacion
for i = 1:steps-1
    % --- PLANTA NO LINEAL (REALIDAD) ---
    %entrada real = u_eq + perturbacion desconocida
    u_t = u_eq;
    if time(i) >= tiempo_step   
        u_t=u_eq+perturbacion;
    end
    
    x_k = X_real(:, i);
    dx = modelo_nolineal(x_k, u_t, params);
    X_real(:, i+1) = x_k + dt * dx;
    
    % --- MEDICION ---
    %salida y (Linealizada/Desviación) para el observador
    %el observador trabaja con dy=y_medida-y_equilibrio
    y_absoluta=C*X_real(:, i)+(ruido_std*rand);    %y medida C*X_real (lo medido)
    y_eq_val=C*p_eq';             %y en p_eq (x_eq)
    dy=y_absoluta-y_eq_val;   %diferencia de y (dy)
    
    % --- OBSERVADOR OEDL ---
    %dz = N*z + L*dy + G*du_conocida
    %asumimos que la parte conocida u es 0
    %TODA la desviación es desconocida (v).
    du_conocida=0; 
    
    dz=N*z+L*dy+G*du_conocida;
    z=z+dt*dz;
    
    %recuperación del estado estimado: dx_hat = z - E*dy
    dx_hat=z-E*dy;
    
    %recuperacion estado absoluto estimado (sin desviacion)
    X_est(:, i+1) = dx_hat + p_eq';
end

% --- GRAFICAS ---
nombres = {'Biomasa (x1)', 'Agua Planta (x2)', 'Agua Suelo (x3)'};
titulos_fig = {'Figura 1: Biomasa (x1) - Real vs. Estimado',...
    'Figura 2: Agua Planta (x2) - Real vs. Estimado',...
    'Figura 3: Agua Suelo (x3) - Real vs. Estimado'};

%bucle para generar graficas independientes
for k=1:3
    figure('Color','w'); %nueva figura para cada estado
    hold on; grid on;
    
    % Graficar la respuesta real y la estimada
    plot(time, X_real(k,:), 'b', 'LineWidth', 1.5, 'DisplayName', 'Real (NL)');
    plot(time, X_est(k,:), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Estimado OEDL');
    
    % Títulos y Etiquetas
    title(titulos_fig{k});
    ylabel(nombres{k});
    xlabel('Tiempo (días)'); % Label de tiempo en cada gráfica
    
    % Línea de perturbación y Leyenda
    xline(tiempo_step, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Inicio Perturbación');
    legend('show', 'Location', 'best');
end

%graficas de error
error_abs = abs(X_real - X_est); 
error_nombres = {'Error Absoluto e1 (Biomasa)',... 
    'Error Absoluto e2 (Agua Planta)',... 
    'Error Absoluto e3 (Agua Suelo)'};

for k = 1:3
    figure('Color','w'); %nueva figura por cada error
    
    %grafica la k-esima fila del k-esimo estado
    plot(time, error_abs(k, :), 'LineWidth', 1.5);
    
    % Títulos y Etiquetas
    title(['Figura ' num2str(3+k) ': ' error_nombres{k}]); 
    ylabel(['|' nombres{k} ' Estimado - Real|']);
    xlabel('Tiempo (días)');
    
    legend(['Error ' num2str(k)], 'Location', 'best'); 
    grid on;
end


% --- FUNCIONES AUXILIARES ---
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