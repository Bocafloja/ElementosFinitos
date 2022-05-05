%% Técnicas de aproximación
%  Ricardo Nájera
%  Alberto Badel
%  2022

close all;
clear all;
clc;
syms a b x xi xii y;

dominio = [2 10];
y_analitica = (149*x^7 + 9980800)/(2499968*x^2 );
y_prueba = a*x^3 + b*x^2 + (5/8 - 124*a - 12*b)*x + (240*a + 20*b - 1/4);
R = 20*a*x^4 - 2976*a*x^2 + 4800*a*x + 24*b*x^3 - 288*b*x^2 + 400*b*x + 15*x^2 - 5*x;

%% Método de colocación
y_colocacion = subs(y_prueba, solve(subs(R,x,7)==0, subs(R,x,3)==0));

%% Método mínimos cuadrados
y_minimos = subs(y_prueba, solve(int(diff(R,a)*R, dominio)==0,int(diff(R,b)*R, dominio)==0));

%% Método de Galerkin
y_Galerkin = subs(y_prueba,solve(int(diff(y_prueba,a)*R, dominio)==0,int(diff(y_prueba,b)*R, dominio)==0));

%% Método de elementos finitos

% Función de forma del elemento lineal
H1 = (xii - x)/(xii-xi);
H2 = (x - xi)/(xii - xi);

% Función de peso y función de prueba
w = [H1;H2];
dw = diff(w,x);
y = [H1 H2];
dy = diff(y,x);

% Características del dominio
a = 2; %Valor x inicial
b = 10; %Valor x final
n = 16; % Número de elementos 
h = (b-a)/n; % Ancho de elemento

% Ecuación diferencial
funcion_K = int(((2*x^3)*dw*dy + (6*x^2)*w*dy + (4*x^2)*w*dy + (20*x)*w*y),x,xi,xii);

K_Global = zeros(n+1);
F_Global = zeros(n+1,1);

xi = a;
xii = h;

for t=1:n

    K_local = eval(funcion_K);

    K_Global(t,t) = K_Global(t,t) + K_local(1,1);
    K_Global(t+1,t) = K_Global(t+1,t) + K_local(2,1);
    K_Global(t,t+1) = K_Global(t,t+1) + K_local(1,2);
    K_Global(t+1,t+1) = K_Global(t+1,t+1) + K_local(2,2);

    xi = xi + h;
    xii = xii + h;

end

% Condiciones de frontera

K_Global(1,:) = zeros(n+1,1);
K_Global(1,1) = 1;
K_Global(n+1,:) = zeros(n+1,1);
K_Global(n+1,n+1) = 1;

F_Global(1,1) = 1;
F_Global(n+1,1) = 6;

xGraphLineal = 2:(10-2)/n:10;
yElemFinitos_lineal = K_Global\F_Global; % == yElemFinitos_lineal = inv(K_Global)*F_Global
 %% Método de elementos finitos (elemento cuadratico)

% Función de forma del elemento cuadrático
H1 = (((-xi^2)*(xii^2 - x^2))-2*xii^2)/(xi^2*(xi^2-xii^2));
H2 = (xi^2 - x^2)/(xi^2 + xii^2);

% Función de peso y función de prueba
w = [H1;H2];
dw = diff(w,x);
y = [H1 H2];
dy = diff(y,x);

% Características del dominio (elemento cuadrático)
n = 3; % Número de elementos 
h = (b-a)/n; % Ancho de elemento

% Ecuación diferencial
funcion_K = int(((2*x^3)*dw*dy + (6*x^2)*w*dy + (4*x^2)*w*dy + (20*x)*w*y),x,xi,xii);

K_Global = zeros(n+1);
F_Global = zeros(n+1,1);

xi = a;
xii = h;

for t=1:n

    K_local = eval(funcion_K);

    K_Global(t,t) = K_Global(t,t) + K_local(1,1);
    K_Global(t+1,t) = K_Global(t+1,t) + K_local(2,1);
    K_Global(t,t+1) = K_Global(t,t+1) + K_local(1,2);
    K_Global(t+1,t+1) = K_Global(t+1,t+1) + K_local(2,2);

    xi = xi + h;
    xii = xii + h;

end

% Condiciones de frontera

K_Global(1,:) = zeros(n+1,1);
K_Global(1,1) = 1;
K_Global(n+1,:) = zeros(n+1,1);
K_Global(n+1,n+1) = 1;

F_Global(1,1) = 1;
F_Global(n+1,1) = 6;

xGraphCuadr = 2:(10-2)/n:10;
yElemFinitos_Cuadratico = K_Global\F_Global; % == yElemFinitos_Cuadratico = inv(K_Global)*F_Global

%% Graficado
grafica1 = figure('Name','Comparación de metodos de aproximación','NumberTitle','off');
fplot(y_analitica, dominio, 'k--')
hold on
fplot([y_colocacion y_minimos y_Galerkin], dominio)
plot(xGraphLineal,yElemFinitos_lineal,xGraphCuadr,yElemFinitos_Cuadratico)
legend('Sol analítica', 'Colocación', 'Mínimos', 'Galerkin','Elementos Finitos Lineal','Elementos Finitos Cuadrático','Location','northwest')
title('Comparación de metodos de aproximación')
xlabel('x','FontSize',14,'Interpreter','latex')
ylabel('y', 'FontSize',14,'Interpreter','latex')

grafica2 = figure('Name','Comparación elementos finitos lineal y cuadrático','NumberTitle','off');
fplot(y_analitica, dominio, 'k--')
hold on
plot(xGraphLineal,yElemFinitos_lineal,xGraphCuadr,yElemFinitos_Cuadratico)
legend('Sol analítica','Elementos Finitos Lineal','Elementos Finitos Cuadrático','Location','northwest')
title('Comparación elementos finitos lineal y cuadrático')
xlabel('x','FontSize',14,'Interpreter','latex')
ylabel('y', 'FontSize',14,'Interpreter','latex')

%% Calculo de errores
x = 2:((10-2)/100):10;
y_analitica_err = (149*x.^7 + 9980800)./(2499968*x.^2 );
yColError = polyval(sym2poly(y_colocacion),x);
yMinError = polyval(sym2poly(y_minimos),x);
yGalError = polyval(sym2poly(y_Galerkin),x);

RMSE_col = sqrt(mean((y_analitica_err - yColError).^2))
RMSE_min = sqrt(mean((y_analitica_err - yMinError).^2))
RMSE_gal = sqrt(mean((y_analitica_err - yGalError).^2))

y_analitica_errLin = (149*xGraphLineal.^7 + 9980800)./(2499968*xGraphLineal.^2 );
RMSE_eflin = sqrt(mean((y_analitica_errLin - yElemFinitos_lineal').^2))

y_analitica_errCuad = (149*xGraphCuadr.^7 + 9980800)./(2499968*xGraphCuadr.^2 );
RMSE_efcuad = sqrt(mean((y_analitica_errCuad - yElemFinitos_Cuadratico').^2))










