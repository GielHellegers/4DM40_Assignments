clear all; clc;
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
fileID = fopen('TH1_Exc1.txt','r');
formatSpec = '%s';
long_str = fscanf(fileID, formatSpec);
long_str = strrep(long_str,'A','a1a');
long_str = strrep(long_str,'D','d0d');
get_all_nrs = str2double(regexp(long_str,'\d+','match'));
tau = get_all_nrs(1:3:end);
i = get_all_nrs(2:3:end);
ev = get_all_nrs(3:3:end);
data = [tau;i;ev];

%% Distribution calculations
z=1; xs = []; dEPT = []; dOvert = []; % Initialising
while z < size(data,2)+1
    datatemp = data(:,z); % read event data
    tau = datatemp(1); i = datatemp(2); ev = datatemp(3);
    if ev==1 % arrival
        if isempty(xs)==1
            [s,sw]= deal(tau,1); % Edgecase: system idle between lots
        end
        xs = [xs; [i, size(xs,1)]];
    elseif ev==0 % departure
        dEPT = [dEPT; [tau-s sw]]; % EPT distribution
        [xs, k, aw] = detOvert(xs, i);
        dOvert = [dOvert; [k aw]]; % Overtake distribution
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

function [xs1,k1,aw1]=detOvert(xs, i)
    ys= [];
    while isempty(xs) == 0
        j= xs(1,1);
        aw = xs(1,2);
        xs(1,:)=[];
        if j<i
            ys= [ys;[j, aw]];
        elseif j== i
            xs1 = [ys;xs];
            k1 = size(ys,1);
            aw1 = aw;
            return
        else
            return
        end
    end
end
        