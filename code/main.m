 %% Main
% Simplement per fer correr els tests

clc; clear; close all;

% file1 = 'dades_stiffness.m';
% test = TestSolver().create('STIFFNESS', file1);
% test.checkPassed(file1)

% file2 = 'dades.m';
% test2 = TestSolver().create('FEM', file2);
% test2.checkPassed(file2)

% z = [zeros(2,1); 0.8*ones(5,1)];
% S = [ones(11,1); 2*ones(6,1)];

file1 = 'fixedData.m';
file2 = 'variableInitialData.m';
opt   = 1;
[displ, stress] = calculateResults(file1,file2,opt, 1);

function [displ, stress] = calculateResults(file1,file2,opt, graph)
    pre = PreProcessor(file1);
    pre.computeInitialData(file2);
    pre.compute();
    pre.computeForceCase(opt);

    s.data = pre.data;
    s.dim  = pre.dim;
    s.solvertype = 'DIRECT'; % ITERATIVE
    FEM = FEMAnalyzer(s);
    FEM.perform();
    displ = FEM.displacement;
    stress = FEM.stress;
    s.FEM = FEM;
    s.graph = graph;
    post = PostProcessor(s);
end