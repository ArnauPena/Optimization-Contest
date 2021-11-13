%% Material
% Material properties
E = 200e9; % Pa

% Cross-section parameters
d = 110e-3; % m
tA = 5e-3; % m
a = 25e-3; % m
b = 75e-3; % m
tB = 3e-3; % m

% Mat parameters
Dext = d+tA/2;
Dint = d-tA/2;

AA = pi*((Dext/2)^2 - (Dint/2)^2);
IzA = pi/4*((Dext/2)^4 - (Dint/2)^4);

AB = 2*a*tB+tB*(b-tB);
IzB = 1/4*( ( (b-tB)^3*tB + 2*tB^3*a)/3 + 2*b^2*tB*a);


% Other data
T = 15e3; % N
Vto = 70; % m/s
dto = 65; % m
Mt = 15800; % kg
sig_rot = 230e6; % MPa
g = 9.81; % m/s2

%% Geometric

% Geometric data
H1 = 0.22; % m
H2 = 0.36; % m
H3 = 0.48; % m
H4 = 0.55; % m
H5 = 2.1; % m
L1 = 0.9; % m
L2 = 0.5; % m
L3 = 0.5; % m
L4 = 2.5; % m
L5 = 4.2; % m

% Nodal coordinates
data.nodes = [
      0,          H1; % Node 1
      0,       H1+H2; % Node 2
      0,    H1+H2+H3; % Node 3
      0, H1+H2+H3+H4; % Node 4
    -L1, H1+H2+H3+H4; % Node 5
];

% Nodal connectivities  
data.nodalconnec = [
    1, 2; % Element 1-2
    2, 3; % Element 2-3
    3, 4; % Element 3-4
    3, 5; % Element 3-5
];

% Material properties matrix
data.materials = [
               E,               AA,        IzA;  % Material 1
               E,               AB,        IzB;  % Material 2
];

% Material connectivities
data.matconnec = [
    1; % Element 1-2 / Material 1 (A)
    1; % Element 2-3 / Material 1 (A)
    1; % Element 3-4 / Material 1 (A)
    2; % Element 3-5 / Material 2 (B)
];


%% Boundary

Pmax = 5.8054e+05;
F1 =  -3.0527e+05;

% Cas 1 - Immediately before 1-6 breaks
data.fdata = [
    1, 1, F1;
    2, 1, Pmax/2;
    ];

% Cas 2 - After 1-6 breaks
data.fdata = [
    1, 1, 0;
    2, 1, Pmax;
    ];

%Forces i nodes fixats
data.fixnod = [
    4, 1, 0;
    4, 2, 0;
    4, 3, 0; % s'han d'imposar els tres DOFs...
    5, 1, 0;
    5, 2, 0;
    5, 3, 0; % s'han d'imposar els tres DOFs...
    ];

%% Problem data
% Dimensions
dim.nd = size(data.nodes,2);   % Problem dimension
dim.nel = size(data.nodalconnec,1); % Number of elements (bars)
dim.nnod = size(data.nodes,1); % Number of nodes (joints)
dim.nne = size(data.nodalconnec,2); % Number of nodes in a bar
dim.ni = 3;           % Degrees of freedom per node
dim.ndof = dim.nnod*dim.ni;  % Total number of degrees of freedom

% Solver to use
% solvertype = 'DIRECT';

% Expected results
results = [0.3591; 0.0024; 0.4474; 0.1980; 0.0024; 0.4474; 0.0242; 0.0024; 0.1917; 0; 0; 0; 0; 0; 0];
% results = [0.359082238819108;0.00237759060308813;0.447441449389525;0.198003317038879;0.00237759060308813;0.447441449389524;0.0241548862050407;0.00237759060308813;0.191669793932445;0;0;0;0;0;0];

connectivities = [1     2     3     4     5     6;
     4     5     6     7     8     9;
     7     8     9    10    11    12;
     7     8     9    13    14    15;];
 
load dades_stiffness.mat KComp;

s.dim         = dim;
s.data.fixnod = data.fixnod;
DOFManager = DOFFixer(s);
DOFManager.fix();