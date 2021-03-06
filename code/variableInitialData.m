% PROBLEM DATA 
Nz_el  = 10;                             % Number of shape variables
Nel    = 260;

% SHAPE VARIABLES
data.z_node = randsample(28500,Nz_el,true)-25000;    % Initial bottom nodes z coor (mm)
% data.z_node = ones(10,1)*(-1000);
% SIZE VARIABLES
data.Sel  = randsample(37,260,true);