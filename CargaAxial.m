function CargaAxial

clc
close all
clear all

Conect = load('Conectivi.txt');
CantElem = size(Conect,1);
Coord = load('Coordenadas.txt');
CantNodo = size(Coord,1);
Proper = load('Propiedades.txt');

KGlobal = zeros(CantNodo);
FGlobal = zeros(CantNodo, 1);

for t = 1:1:CantElem
    
    Le = ( Coord( Conect(t,2), 1) - Coord( Conect(t,1), 1) );
    KLocal = ( Proper(t,2)* Proper(t,1)/Le) * [1 -1; -1 1];
    
    KGlobal(Conect(t,1), Conect(t,1)) = KGlobal(Conect(t,1), Conect(t,1)) + KLocal(1,1);
    KGlobal(Conect(t,1), Conect(t,2)) = KGlobal(Conect(t,1), Conect(t,2)) + KLocal(1,2);
    KGlobal(Conect(t,2), Conect(t,1)) = KGlobal(Conect(t,2), Conect(t,1)) + KLocal(2,1);
    KGlobal(Conect(t,2), Conect(t,2)) = KGlobal(Conect(t,2), Conect(t,2)) + KLocal(2,2);
    
end

Fuerz = [3 240000];

FGlobal(Fuerz(1,1), 1) = Fuerz(1,2);

Front = load('CondiF.txt');
CantFront = size(Front,1);

C = 10000000000;

KPena = KGlobal;
FPena = FGlobal;

for d= 1:1:CantFront
    
    KPena(Front(d,1), Front(d,1)) = KPena(Front(d,1), Front(d,1)) + C;
    FPena(Front(d,1), 1) = FPena(Front(d,1), 1) + C*Front(d,2);
    
end

U = inv(KPena)*FPena














