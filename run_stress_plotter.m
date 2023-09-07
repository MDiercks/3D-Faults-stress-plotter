% 3D-Fault Stress plotter app
% 
% tool to calculate cumulative stress from coseismic and interseismic stress (Coulomb3.4 Element Conditions files)
% - plots stress on fault network (using fault geometry exported by 3D-Faults v.2.4+)
% - automatically exports figures and stress statistics
% - combines different scenarios for earthquakes
% version 1.0 // 09/2023 // Manuel Diercks
% requirements: Coulomb v3.3, 3D-Faults v2.4, MATLAB R2020a or later versions
% NOTE: all stresses are converted to MPa
% 
% Inputs:
% - directory containing all coseismic .csv files calculated with Coulomb 3.3+
% - directory containing fault geometry files exported from 3D-Faults v2.4+
% - .csv Coulomb file with annual interseismic CST loading
%
% (all coseismic and interseismic CST must be modelled on the exact same fault network, check the number of elements in coulomb files)
%
% to run the app, press F5 or type run_stress_plotter in the command window
%##########################################################################
addpath source\
stress_plotter