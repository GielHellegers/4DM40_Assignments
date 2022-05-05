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

%% Determining initial wip
% Choosing data source
data = Data_Workstation1_extracted;
% Determining wip
departinglots = data(2,find(data(3,:)==0));
arrivinglots = data(2,find(data(3,:)==1));
initialwiplots = setdiff(departinglots,arrivinglots);
initialwip = size(initialwiplots,2);

%% Finding uncertainties
% Lots with uncertain EPTs
if isempty( setdiff(data(2,min(find(data(3,:)==0))),initialwiplots) ) == 1 % if first departure is in initialwip
    uncertainEPT = data(2,min(find(data(3,:)==0))); % EPT of entry is unknown due to unknown previous EPT start
else
    uncertainEPT = 0; % EPTs are certain for all lots arriving or departing the system
end

% Lots with uncertain Overtake
uncertainOvert = initialwiplots; % Overtake is uncertain for initial lots, since the order of arrival is unknown

%% Initialising xs
xs = [];
while isempty(initialwiplots) == 0
    xs = [xs; initialwiplots(1),-1]; % Uncertain order of arrivals
    initialwiplots(1)=[];
end

%% Distribution calculations
z=1; dEPT = []; dOvert = []; % Initialising
while z < size(data,2)+1
    datatemp = data(:,z); % read event data
    tau = datatemp(1); i = datatemp(2); ev = datatemp(3);
    if ev==1 % arrival
        if z==1 & isempty(xs)==0 % An initial wip exists, xs is not empty at first event
            [s,sw]=deal(tau,1+size(xs,1));
        elseif isempty(xs)==1
            [s,sw]= deal(tau,1); % Edgecase: system idle between lots
        end
        xs = [xs; [i, size(xs,1)]];
    elseif ev==0 % departure
        if isempty( setdiff(i,uncertainEPT) ) == 0 % If its EPT is uncertain, exclude from distribution
            dEPT = [dEPT; [tau-s sw]]; % EPT distribution
        end
        [xs, k, aw] = detOvert(xs, i, uncertainOvert);
        if isempty( setdiff(i,uncertainOvert) ) == 0 % If its Overtake is uncertain, exclude from distribution
            dOvert = [dOvert; [k aw]]; % Overtake distribution
        end
        if isempty(xs)==0
            s = tau;
            sw = size(xs,1);
        end
    end
    z=z+1; % Next event
end
% Print results
dEPT
dOvert

%% make distributions
max(dEPT(:,2));
min(dEPT(:,2));
EPT_d = [];
for w = min(dEPT(:,2)) : max(dEPT(:,2));
    w_spec_matrix = dEPT(find(dEPT(:,2) == w),:);
    w_spec_matrix = w_spec_matrix.';
    w_spec_array = w_spec_matrix(1,:);
    EPT_d = [EPT_d;w, {w_spec_array}] ;
end

max(dOvert(:,2));
min(dOvert(:,2));
Overt_d = [];
for w = min(dOvert(:,2)) : max(dOvert(:,2));
    w_spec_matrix = dOvert(find(dOvert(:,2) == w),:);
    w_spec_matrix = w_spec_matrix.';
    w_spec_array = w_spec_matrix(1,:);
    Overt_d = [Overt_d;w, {w_spec_array}];
end



EPT_d;
Overt_d;


rows_EPT_d = size(EPT_d,1);
for i = 1:rows_EPT_d
    array_to_inv = EPT_d{i,2};
    length_array = size(EPT_d{i,2},2);
    memory_list = [];
    for j = array_to_inv
        if ismember(j,memory_list) == 0;
            memory_list = [memory_list,j];
            count_j = sum(array_to_inv==j);
            distribution_j = count_j/length_array;
        end
        disp([num2str(j), ' with a distribution of ', num2str(distribution_j), ' for a wip of ', num2str(i)])
    end
end

rows_Overt_d = size(Overt_d,1);
for i = 1:rows_Overt_d
    array_to_inv = Overt_d{i,2};
    length_array = size(Overt_d{i,2},2);
    memory_list = [];
    for j = array_to_inv
        if ismember(j,memory_list) == 0;
            memory_list = [memory_list,j];
            count_j = sum(array_to_inv==j);
            distribution_j = count_j/length_array;
        end
        disp([num2str(j), ' with a distribution of ', num2str(distribution_j), ' for a wip of ', num2str(i)])
    end
end

%% determine mean EPT and overtake
means = mean(dEPT);
mean_dEPT = means(1)
stds = std(dEPT);
std_dEPT = stds(1)
c_of_var_dEPT = std_dEPT/mean_dEPT


means = mean(dOvert);
mean_dOvert = means(1)
stds = std(dOvert);
std_dOvert = stds(1)
c_of_var_dOvert = std_dOvert/mean_dOvert

%% plotting distributions
%Vars to change
mu = mean_dOvert;                                %// Mean
sigma = std_dOvert;                            %// Standard deviation

%// Plot curve
x = (-5 * sigma:0.01:5 * sigma) + mu;  %// Plotting range
y = exp(- 0.5 * ((x - mu) / sigma) .^ 2) / (sigma * sqrt(2 * pi));
plot(x, y)

%// Hide ticks
set(gca, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', [])

%% functions
function [xs1,k1,aw1]=detOvert(xs, i, uncertainOvert)
    ys= []; kc = 0;
    if isempty( setdiff(i,uncertainOvert) ) == 1 % If its Overtake is uncertain, return
        xs(1,:)=[];
        xs1 = [ys;xs];
        k1 = size(ys,1);
        aw1 = -1;
        return
    end
    while isempty(xs) == 0
        j= xs(1,1);
        aw = xs(1,2);
        xs(1,:)=[];
        if j<i
            ys = [ys;[j, aw]];
        elseif j == i
            xs1 = [ys;xs];
            k1 = size(ys,1)-kc;
            aw1 = aw;
            return
        elseif aw == -1
            ys = [ys;[j, aw]]; % Initial lot overtaken by lot i
        elseif j>i
            ys = [ys; [j,aw]];   % No overtake, lot arrived later
            kc = kc+1; % kc is a variable compensating for j>i producing no overtake
        else
            return
        end
    end
end