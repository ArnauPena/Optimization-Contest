 %% Main
% Simplement per fer correr els tests

clc; clear; close all;

% file1 = 'dades_stiffness.m';
% test = TestSolver().create('STIFFNESS', file1);
% test.checkPassed(file1)

file2 = 'dades.m';
% test2 = TestSolver().create('FEM', file2);
% test2.checkPassed(file2)

z = [zeros(2,1); 0.8*ones(5,1)];
S = [ones(11,1); 2*ones(6,1)];
[displ, stress] = calculateResults(z, S, file2);

function [displ, stress] = calculateResults(zeds, sects, filename)
    run(filename)
    data.matconnec = sects;
    data.nodes(:,3) = zeds;
    s.dim        = dim;
    s.data       = data;
    s.solvertype = 'DIRECT'; % ITERATIVE
    FEM = FEMAnalyzer(s);
    FEM.perform();
    displ = FEM.displacement;
    stress = FEM.stress;
end