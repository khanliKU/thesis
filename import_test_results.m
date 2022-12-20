%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: C:\Users\kutla\Documents\GitHub\thesis\test_results.csv
%
% Auto-generated by MATLAB on 09-Dec-2022 22:38:59

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 93);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Zaman_1_1_sec", "Kuvvet_1_1_N", "Uzama_1_1_mm", "Zaman_1_2_sec", "Kuvvet_1_2_N", "Uzama_1_2_mm", "Zaman_1_3_sec", "Kuvvet_1_3_N", "Uzama_1_3_mm", "Zaman_1_4_sec", "Kuvvet_1_4_N", "Uzama_1_4_mm", "Zaman_1_5_sec", "Kuvvet_1_5_N", "Uzama_1_5_mm", "Zaman_2_v20_sec", "Kuvvet_2_v20_N", "Uzama_2_v20_mm", "Zaman_no_sec", "Kuvvet_no_N", "Uzama_no_mm", "Zaman_2_v30_sec", "Kuvvet_2_v30_N", "Uzama_2_v30_mm", "Zaman_2_v40_sec", "Kuvvet_2_v40_N", "Uzama_2_v40_mm", "Zaman_2_v10_sec", "Kuvvet_2_v10_N", "Uzama_2_v10_mm", "Zaman_3_1_sec", "Kuvvet_3_1_N", "Uzama_3_1_mm", "Zaman_3_2_sec", "Kuvvet_3_2_N", "Uzama_3_2_mm", "Zaman_3_3_sec", "Kuvvet_3_3_N", "Uzama_3_3_mm", "Zaman_3_none_sec", "Kuvvet_3_none_N", "Uzama_3_none_mm", "Zaman_3_20v_sec", "Kuvvet_3_20v_N", "Uzama_3_20v_mm", "Zaman_4_10v_sec", "Kuvvet_4_10v_N", "Uzama_4_10v_mm", "Zaman_4_10v_sec1", "Kuvvet_4_10v_N1", "Uzama_4_10v_mm1", "Zaman_4_10v_sec2", "Kuvvet_4_10v_N2", "Uzama_4_10v_mm2", "Zaman_4_15v_sec", "Kuvvet_4_15v_N", "Uzama_4_15v_mm", "Zaman_4_15v_sec1", "Kuvvet_4_15v_N1", "Uzama_4_15v_mm1", "Zaman_5_25v_sec", "Kuvvet_5_25v_N", "Uzama_5_25v_mm", "Zaman_5_25v_sec1", "Kuvvet_5_25v_N1", "Uzama_5_25v_mm1", "Zaman_5_30v_sec", "Kuvvet_5_30v_N", "Uzama_5_30v_mm", "Zaman_5_35v_sec", "Kuvvet_5_35v_N", "Uzama_5_35v_mm", "Zaman_5_35v_sec1", "Kuvvet_5_35v_N1", "Uzama_5_35v_mm1", "Zaman_6_40v_sec", "Kuvvet_6_40v_N", "Uzama_6_40v_mm", "Zaman_6_40v_sec1", "Kuvvet_6_40v_N1", "Uzama_6_40v_mm1", "Zaman_6_45v_sec", "Kuvvet_6_45v_N", "Uzama_6_45v_mm", "Zaman_6_50v_sec", "Kuvvet_6_50v_N", "Uzama_6_50v_mm", "Zaman_6_55v_sec", "Kuvvet_6_55v_N", "Uzama_6_55v_mm", "Zaman_7_60v_sec", "Kuvvet_7_60v_N", "Uzama_7_60v_mm"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
testresults = readtable("C:\Users\kutla\Documents\GitHub\thesis\test_results.csv", opts);


%% Clear temporary variables
clear opts