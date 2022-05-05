clc;
clear all;
close all;

coordenada = load('CoordenadasMATLAB.txt');
conectividad = load('Conectividad.txt');
Propiedad = load('Propiedad.txt');
Fuerza = load('Fuerzas.txt');
Fronteras = load('Fronteras.txt');
resultANSYS = load('ResultAnsys.txt');

NumNodos = size(coordenada,1);
NumElementos = size(conectividad,1);

KGlobal = zeros(NumNodos*2);
FGlobal = zeros(NumNodos*2,1);

for i=1:NumElementos
    L = sqrt((coordenada(conectividad(i,2),1)-coordenada(conectividad(i,1),1))^2 + (coordenada(conectividad(i,2),2)-coordenada(conectividad(i,1),2))^2);
    l = (coordenada(conectividad(i,2),1)-coordenada(conectividad(i,1),1))/L ;
    m = (coordenada(conectividad(i,2),2)-coordenada(conectividad(i,1),2))/L ;
    Klocal = (Propiedad(i,1)*Propiedad(i,2)/L)*[l^2 l*m -l^2 -l*m; l*m m^2 -l*m -m^2; -l^2 -l*m l^2 l*m; -l*m -m^2 l*m m^2];

    KGlobal(conectividad(i,1)*2 - 1,conectividad(i,1)*2-1) = KGlobal(conectividad(i,1)*2 - 1,conectividad(i,1)*2-1) + Klocal(1,1);
    KGlobal(conectividad(i,1)*2 - 1,conectividad(i,1)*2) = KGlobal(conectividad(i,1)*2 - 1,conectividad(i,1)*2) + Klocal(1,2);
    KGlobal(conectividad(i,1)*2,conectividad(i,1)*2-1) = KGlobal(conectividad(i,1)*2,conectividad(i,1)*2-1) + Klocal(2,1);
    KGlobal(conectividad(i,1)*2,conectividad(i,1)*2) = KGlobal(conectividad(i,1)*2,conectividad(i,1)*2) + Klocal(2,2);    

    KGlobal(conectividad(i,1)*2 - 1,conectividad(i,2)*2 - 1) = KGlobal(conectividad(i,1)*2 - 1,conectividad(i,2)*2 - 1) + Klocal(1,3);
    KGlobal(conectividad(i,1)*2 - 1,conectividad(i,2)*2) = KGlobal(conectividad(i,1)*2 - 1,conectividad(i,2)*2) + Klocal(1,4);
    KGlobal(conectividad(i,1)*2 ,conectividad(i,2)*2 - 1) = KGlobal(conectividad(i,1)*2 ,conectividad(i,2)*2 - 1) + Klocal(2,3);
    KGlobal(conectividad(i,1)*2 ,conectividad(i,2)*2) = KGlobal(conectividad(i,1)*2 ,conectividad(i,2)*2) + Klocal(2,4);

    KGlobal(conectividad(i,2)*2 - 1 , conectividad(i,1)*2 - 1) = KGlobal(conectividad(i,2)*2 - 1 , conectividad(i,1)*2 - 1) + Klocal(3,1);
    KGlobal(conectividad(i,2)*2 - 1 , conectividad(i,1)*2) = KGlobal(conectividad(i,2)*2 - 1 , conectividad(i,1)*2) + Klocal(4,1);
    KGlobal(conectividad(i,2)*2 , conectividad(i,1)*2 - 1) = KGlobal(conectividad(i,2)*2 , conectividad(i,1)*2 - 1) + Klocal(3,2);
    KGlobal(conectividad(i,2)*2 , conectividad(i,1)*2) = KGlobal(conectividad(i,2)*2 , conectividad(i,1)*2) + Klocal(4,2);

    KGlobal(conectividad(i,2)*2 - 1, conectividad(i,2)*2 - 1) = KGlobal(conectividad(i,2)*2 - 1, conectividad(i,2)*2 - 1) + Klocal(3,3);
    KGlobal(conectividad(i,2)*2 - 1, conectividad(i,2)*2) = KGlobal(conectividad(i,2)*2 - 1, conectividad(i,2)*2) + Klocal(3,4);
    KGlobal(conectividad(i,2)*2, conectividad(i,2)*2 - 1) = KGlobal(conectividad(i,2)*2, conectividad(i,2)*2 - 1) + Klocal(4,3);
    KGlobal(conectividad(i,2)*2, conectividad(i,2)*2) = KGlobal(conectividad(i,2)*2, conectividad(i,2)*2) + Klocal(4,4);

end

NumFuerzas = size(Fuerza,1);

for j=1:NumFuerzas
   FGlobal(Fuerza(j,1)*2 - 1,1) = Fuerza(j,2);
   FGlobal(Fuerza(j,1)*2,1) = Fuerza(j,3);
end

%% Enfoque de penalización

C = 6000000000000;
KPenal = KGlobal;
FPenal = FGlobal;

NumFronteras = size(Fronteras,1);

for m = NumFronteras:-1:1
    KPenal(Fronteras(m,1),Fronteras(m,1)) = KPenal(Fronteras(m,1),Fronteras(m,1)) + C;
    FPenal(Fronteras(m,1),1) = FPenal(Fronteras(m,1),1) + C*Fronteras(m,2);
end

uPenal = inv(KPenal)*FPenal;
xDef = zeros(15,1);
yDef = zeros(15,1);
t=1;
s=1;
for i=1:30
    if rem(i,2)==0
        yDef(s,1) = uPenal(i,1);
        s= s+1;
    end
    if rem(i,2)~=0
        xDef(t,1) = uPenal(i,1);
        t = t+1;
    end
end

%% Calculo de error relativo

xDefSum = 0;
yDefSum = 0;
errAbsX = 0;
errAbsY = 0;
for i=1:15
    xDefSum = xDefSum + xDef(i,1);
    yDefSum = yDefSum + yDef(i,1);
    errAbsX = errAbsX + (xDef(i,1)-resultANSYS(i,1));
    errAbsY = errAbsY + (yDef(i,1)-resultANSYS(i,2));
end

Dx = errAbsX/xDefSum; % Error relativo respecto a X
Dy = errAbsY/yDefSum; % Error relativo respecto a Y




