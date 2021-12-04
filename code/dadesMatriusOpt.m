%%%%%%%%%%% DADES MATRIUS OPTIMITZACIÓ %%%%%%%%%%%%
clear; clc;
%%%%%% Està fet bastant cutre, només per tenir una idea de com seria la
%%%%%% definició de les matrius que vam comentar l'altre dia (està fet pel
%%%%%% cas senzillet de 3 elements i 4 nodes que vam parlar però es pot
%%%%%% escalar sense problemes

nel         = 3;
nnodes      = 4;
nSections   = 37;
Seq         = sparse(nel,nSections*nel);
nValues     = ones(1,nSections);
z0          = linspace(-100,1000,nnodes); % inventadissim

% Definició vector x
x2          = sparse(nSections*nel,1);
x1          = z0';
k           = nnodes + 1;
x2ID        = randsample(37,nel);
for iElem = 1:nel
    n1     = (iElem-1)*nSections + x2ID(iElem);
    x2(n1) = 1;
end

% Definició matriu Aeq 
for iElem = 1:nel
    n1 = (iElem-1)*nSections + 1;
    n2 = (n1-1) + nSections;
    Seq(iElem,n1:n2) = nValues;
end

k   = sparse(nel,nnodes);
Aeq = [k,Seq];
x   = [x1;x2];

% Comprovació
b = Aeq*x;

