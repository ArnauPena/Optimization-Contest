 %% Main
% Simplement per fer correr els tests

clc; clear; close all;

file1 = 'fixedData.m';
file2 = 'variableInitialData.m';
opt   = 1;
[displ, stress] = calculateResults(file1,file2,opt, 'STRESS');
% que retorni M i 

function [displ, stress] = calculateResults(file1,file2,opt, plottingMode)
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
    s.plottingMode = plottingMode;
    post = PostProcessor(s);
end