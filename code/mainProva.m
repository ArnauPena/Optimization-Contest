%%%% PROVA %%%%%%%
clear; clc; close all;
a = PreProcessor('fixedData.m');
a.computeInitialData('variableInitialData.m');
a.compute();
s.data = a.data;
s.dim  = a.dim;
b = Visualizer(s);
b.compute();

