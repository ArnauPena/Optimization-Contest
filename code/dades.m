%% INPUT DATA

% Geometric data
H = 0.8;
W = 0.75;
B = 2.8;
d = 0.5;
D1 = 18e-3; 
d1 = 7.5e-3;
D2 = 3e-3;
I1 = [1e-8; 5e-9; 5e-9;];
I2 = [3e-10; 2e-10; 2e-10;];

% Mass
M = 135;

% Other
g = 9.81;

%% PREPROCESS

% Nodal coordinates matrix 
%  x(a,j) = coordinate of node a in the dimension j
data.nodes = [%X      Y      Z
         2*W,  -W/2,     0; % (1)
         2*W,   W/2,     0; % (2)
         2*W,     0,     H; % (3)
           0,     0,     H; % (4)
           0,    -B,     H; % (5)
           0,     B,     H; % (6)
           W,     0,     H; % (7)
];

% Nodal connectivities  
%  Tn(e,a) = global nodal number associated to node a of element e
data.nodalconnec = [%     a      b
           1,     2; % (1)
           3,     7; % (2)
           7,     4; % (3)
           3,     5; % (4)
           4,     5; % (5)
           3,     6; % (6)
           4,     6; % (7)
           5,     7; % (8)
           6,     7; % (9)
           1,     3; % (10)
           2,     3; % (11)
           1,     4; % (12)
           2,     4; % (13)
           1,     5; % (14)
           2,     6; % (15)
           7,     1; % (16)
           7,     2; % (17)
];

    
% Material properties matrix
%  mat(m,1) = Young modulus of material m
%  mat(m,2) = Section area of material m
%  mat(m,3) = Density of material m
%  --more columns can be added for additional material properties--
data.materials = [% Young M.        Section A.    Density   Mom. of inertia for buckling
            55e9, pi/4*(D1^2-d1^2),       2350,     min(I1);  % Material (1)
           147e9,        pi/4*D2^2,        950,     min(I2);% Material (2)
];

% Material connectivities
%  Tmat(e) = Row in mat corresponding to the material associated to element e 
L = 1.4306e+03;
T = 1.7510e+03;
D = 1.7510e+03;
data.matconnec = [% Mat. index
                   1; % (1)
                   1; % (2)
                   1; % (3)
                   1; % (4)
                   1; % (5)
                   1; % (6)
                   1; % (7)
                   1; % (8)
                   1; % (9)
                   1; % (10)
                   1; % (11)
                   2; % (12)
                   2; % (13)
                   2; % (14)
                   2; % (15)
                   2; % (16)
                   2; % (17)
];

data.fdata = [
    1, 3, -M*g/2;
    2, 3, -M*g/2;
    1, 1, T/2;
    2, 1, T/2;
    3, 3, L/5;
    4, 3, L/5;
    5, 3, L/5;
    6, 3, L/5;
    7, 3, L/5;
    3, 1, -D/5;
    4, 1, -D/5;
    5, 1, -D/5;
    6, 1, -D/5;
    7, 1, -D/5;
];
% Fixed nodes matrix
data.fixnod = [
    4, 1, 0; 
    4, 2, 0; 
    4, 3, 0;
    2, 1, 0; 
    2, 2, 0; 
    5, 3, 0;
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