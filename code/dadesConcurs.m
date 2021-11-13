%%%%%%%%%%%
clear; clc;
% GENERAL PROBLEM DATA
nNodes = 76;
nNodx  = 4;                    % Number of nodes with same x coor
DOFs   = 3;
Nzk    = 10;
nel    = 260;
Ns     = nel;

% INITIAL GEOMETRIC DATA
zk   = randsample(28500,Nzk) - -25000; % Initial bottom nodes z coor (mm) (will change)
zsup = 4000;                           % Fixed upper nodes z coor (mm)
y12  = 0;                              % Right part nodes y coor (mm)
y34  = 4000;                           % Left part nodes y coor (mm)
x    = linspace(0,18000,nNodes/nNodx); % x coor vector (mm)

% ASSEMBLY OF COOR
nodes  = zeros(nNodes,DOFs);
for i = 1:length(x)
    m = 4*(i-1);
    nodes(m+1,1) = x(i);
    nodes(m+1,2) = y12;
    nodes(m+2,1) = x(i);
    nodes(m+2,2) = y12;
    nodes(m+2,3) = zsup;
    nodes(m+3,1) = x(i);
    nodes(m+3,2) = y34;
    nodes(m+3,3) = zsup;
    nodes(m+4,1) = x(i);
    nodes(m+4,2) = y34;
end

for k = 1:Nzk
    m = 4*(k-1);
    n = nNodes-3-m;
    nodes(m+1,3) = zk(k);
    nodes(n,3)   = zk(k);
    nodes(m+4,3) = zk(k);
    nodes(n+3,3) = zk(k);
end

data.nodes = nodes;

% Nodal connectivities 260 elements
% Aixo ja és una fumada grossa

