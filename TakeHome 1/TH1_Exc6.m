clear all; close all; clc;
%% Definitions
% tau: event time,
% ev: event type (arrival 1, departure 0),
% i: lot arrival number,
% aw: number of lots just before arrival,
% xs: list that stores for each lot i and aw,
% s: variable used to store EPT start time,
% sw: number of lots just after EPT start,
% k: number of overtaken lots,

% ys: list that stores part of list xs,
% j: stores lot arrival number;

%% Reading and converting data from .txt file
fileID = fopen('group05.txt','r');
formatSpec = '%s';
long_str = fscanf(fileID, formatSpec);
long_str = strrep(long_str,'Created','a1a');
long_str = strrep(long_str,'W1D','a2a');
long_str = strrep(long_str,'W2D','a3a');
long_str = strrep(long_str,'W3D','a4a');
get_all_nrs = str2double(regexp(long_str,'\d+','match'));
tau = get_all_nrs(1:3:end);
i = get_all_nrs(2:3:end);
ev = get_all_nrs(3:3:end);
data = [tau;i;ev];

%% Extract individual workstation behaviour
% Workstation 1
Data_Workstation1_extracted = data(:,find(data(3,:)<3));
Data_Workstation1_extracted(3,find(Data_Workstation1_extracted(3,:)==2))=0;

% Workstation 2
Data_Workstation2_extracted = data(:,find(data(3,:)<4 & data(3,:)>1));
Data_Workstation2_extracted(3,find(Data_Workstation2_extracted(3,:)==2))=1;
Data_Workstation2_extracted(3,find(Data_Workstation2_extracted(3,:)==3))=0;

% Workstation 3
Data_Workstation3_extracted = data(:,find(data(3,:)>2));
Data_Workstation3_extracted(3,find(Data_Workstation3_extracted(3,:)==3))=1;
Data_Workstation3_extracted(3,find(Data_Workstation3_extracted(3,:)==4))=0;

%% choose workstation
% Choosing data source
data = Data_Workstation1_extracted;

%% determine inter arrival times
arrivingtimes = data(1,find(data(3,:)==1));
int_ar_times_list = [];
for i = 1:(size(arrivingtimes,2)-1)
    int_ar_time = arrivingtimes(i+1)-arrivingtimes(i);
    int_ar_times_list = [int_ar_times_list, int_ar_time];
end
c_of_var_int_ar_time = std(int_ar_times_list)/mean(int_ar_times_list)



















