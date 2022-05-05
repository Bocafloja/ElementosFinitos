%function ElementosFinitos1D

close all;
clear;
clc;

syms x xi xii y


%% Función de forma del elemento lineal
% if f_forma == 'lineal'
H1 = (xii - x)/(xii-xi);
H2 = (x - xi)/(xii - xi);

%% Función de forma del elemento cuadrático
% elseif f_forma == 'cuadrada'
%     H1 = ;
%     H2 = ;
% end
%% Función de peso y función de prueba
w = [H1;H2];
dw = diff(w,x);
y = [H1 H2];
dy = diff(y,x);

%% Características del dominio
a = 2; %Valor x inicial
b = 10; %Valor x final
n = 3; % Número de elementos 
h = (b-a)/n; % Ancho de elemento

%% Ecuación diferencial
% ter1 = (2*x^3)*dw*dy;
% ter2 = (6*x^2)*w*dy;
% ter3 = (4*x^2)*w*dy;
% ter4 = (20*x)*w*y;
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

%% Condiciones de frontera

K_Global(1,:) = zeros(n+1,1);
K_Global(1,1) = 1;
K_Global(n+1,:) = zeros(n+1,1);
K_Global(n+1,n+1) = 1;

F_Global(1,1) = 1;
F_Global(n+1,1) = 6;

yElemFinitos_lineal = inv(K_Global)*F_Global;

xGraph = 2:(10-2)/n:10;
plot(xGraph,yElemFinitos_lineal)



