%Meenakshi M Asokan's code
%August 2023

close all
clearvars
%Download ans save the 'Data' folder in the current working directory;
%The 'Data' folder should have subfolders containing data corresponding to
%each Figure
currentFolder = pwd;
data_filename = fullfile(currentFolder,'/Data/Figure1/behav_meas_allmice_3sessProtocol_22fm20pd.mat');
data = load(data_filename);
behav_meas_allmice = data.behav_meas_allmice_22fm20pd;
% temp = load('behav_meas_allmice_pseudocond.mat');
% behav_meas_allmice = temp.behav_meas_allmice;
fs = 30;
%%