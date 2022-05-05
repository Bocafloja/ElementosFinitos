function ElemenResorte
clear;
close all;
clc;

conect = load('Conectividad.txt');

cantElem = size(conect,1);
cantNodos = 6;

KGlobal = zeros(cantNodos);
FGlobal = zeros(cantNodos,1);

for p=1:1:cantElem

    Klocal = conect(p,1)*[1 -1;-1 1];

    KGlobal(conect(p,2),conect(p,2)) = KGlobal(conect(p,2),conect(p,2)) + Klocal(1,1);
    KGlobal(conect(p,2),conect(p,3)) = KGlobal(conect(p,2),conect(p,3)) + Klocal(1,2);
    KGlobal(conect(p,3),conect(p,2)) = KGlobal(conect(p,3),conect(p,2)) + Klocal(2,1);
    KGlobal(conect(p,3),conect(p,3)) = KGlobal(conect(p,3),conect(p,3)) + Klocal(2,2);

    FGlobal(conect(p,2),1) = FGlobal(conect(p,2),1) + conect(p,4);
    FGlobal(conect(p,3),1) = FGlobal(conect(p,3),1) + conect(p,5);

end

KElimina = KGlobal;
FElimina = FGlobal;
UElimina = zeros(cantNodos,1);
Frontera = [1   0;
            2   0;
            6   0];

cantFronteras = size(Frontera,1);

for t = cantFronteras:-1:1
    Kmult = KElimina(Frontera(t,1),:);
    UElimina = UElimina + (Kmult' * Frontera(t,2));

    KElimina(:,Frontera(t,1)) = [];
    KElimina(Frontera(t,1),:) = [];
    FElimina(Frontera(t,1),:) = [];
    UElimina(Frontera(t,1),:) = [];

end

UElimina = KElimina\(FElimina - UElimina)






