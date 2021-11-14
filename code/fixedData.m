%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% FIXED DATA %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GENERAL PROBLEM DATA
nNodx  = 4;                    % Number of nodes with same x coor
DOFs   = 3;
nel    = 260;
Ns     = nel;

% DIMENSIONS
dim.nd   = 3;                % Problem dimension
dim.nel  = 260;              % Number of elements (bars)
dim.nnod = 76;               % Number of nodes (joints)
dim.nne  = 2;                % Number of nodes in a bar
dim.ni   = 3;                % Degrees of freedom per node
dim.ndof = dim.nnod*dim.ni;  % Total number of degrees of freedom


% INITIAL GEOMETRIC DATA
data.zsup = 4000;                             % Fixed upper nodes z coor (mm)
data.y12  = 0;                                % Right part nodes y coor (mm)
data.y34  = 4000;                             % Left part nodes y coor (mm)
data.x    = linspace(0,18000,dim.nnod/nNodx); % x coor vector (mm)


% FIXED NODES MATRIX
data.fixnod = [
    1,  1, 0; 
    1,  2, 0; 
    1,  3, 0;
    4,  1, 0; 
    4,  2, 0; 
    4,  3, 0;
    73, 1, 0;
    73, 2, 0;
    73, 3, 0;
    76, 1, 0;
    76, 2, 0;
    76, 3, 0;
];

% EXTERNAL FORCE DATA
Fext = zeros(dim.nnod - 4,dim.ni,3);

% CASE 1
% Nodes
Fext(1:2,1,1)   = [2;3];
Fext(3:70,1,1)  = [5:72]';
Fext(71:72,1,1) = [74;75];
% DOFs
Fext(:,2,1)    = ones(size(Fext,1),1);
% Force
Fext(:,3,1)    = 5e3*ones(size(Fext,1),1);

% CASE 2
% Nodes
Fext(1:2,1,2)   = [2;3];
Fext(3:70,1,2)  = [5:72]';
Fext(71:72,1,2) = [74;75];
% DOFs
Fext(:,2,2)     = 2*ones(size(Fext,1),1);
% Force
Fext(:,3,2)     = 1e3*ones(size(Fext,1),1);

% CASE 3
% Nodes
Fext(1:2,1,3)   = [2;3];
Fext(3:70,1,3)  = [5:72]';
Fext(71:72,1,3) = [74;75];
% DOFs
Fext(:,2,3)     = 3*ones(size(Fext,1),1);
% Force
Fext(:,3,3)     = -5e3*ones(size(Fext,1),1);

data.Fext       = Fext;

% MATERIAL AND SECTION DATA
% Materials section and inertia
data.materials.Si  = [
                        ];
data.materials.rho = 7.85e-6; % Density (kg/mm3)
data.materials.E   = 200e3;   % Young Modulus (MPa)

