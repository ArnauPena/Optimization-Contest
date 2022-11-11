%% Augmented Lagrange Data %%
tol       = 1e-15;
l_u       = 0;
l_sig     = 0;
rho_u     = 10;
rho_sig   = 10;
s0        = 0.6*ones(1,345);
% s0        = rand(1,345);
x0        = s0*(37-1) + 1;
x0        = round(x0);
cost      = 0;
const_u   = 0;
const_sig = 0;
lu        = 0;
lsig      = 0;
tau       = 0;
tauV      = 0;
dS        = 1;