% COMPUTE STATIONARITY METHOD AND STATISTICS

%(c) 2015 Claudio S. Quiroga-Lombard, Bernstein Center for Computational
% Neuroscience Heidelberg-Mannheim

% For more information please refer to
% Quiroga-Lombard, Claudio S., Joachim Hass, and Daniel Durstewitz. 
%"Method for stationarity-segmentation of spike train data with application to the Pearson cross-correlation." 
% Journal of neurophysiology 110.2 (2013): 562-572.

% Load Data

clear
L = cell([],[]);

% Load example, spike times in seconds.
% In the example STMtx1 are the spike times in seconds.
% t.mat is the time
% rate.mat is the rate profile

L = ('/home/claudio.quiroga/Claudio/Upload_ST_Method_22_04_15/Stationary_Method_&_Statistics/STMtx1.mat');
load(L);

Time_vector = [];
Time_vector = [300 700]; % Defines the recording time of interest.
    
% STATIONARY CONDITIONS FOR EXP DATA    
    
data_original = [];    
data_original = STMtx1; % STMtx1 is the name of the data loaded in line 10.
clear STMtx1;
                                             
% Apply Stationary Conditions 

numberISI = 10; % length of the sliding window (number of ISI contained in the window)
percentage = .96; % applied constraints, 1 = no constraints.

% For visualization purposes (see Plot below), choose for example, numberISI 10 and 20
% with same percentage, .96. 

T_centered_corr =[];   
ISI_0 = [];
nst_spike_times = [];
    
[T_centered_corr, ISI_0, nst_spike_times] = stationary_conditions(data_original,numberISI,percentage, Time_vector);
        
% ISI_0 are the original inter-spike intervals, without applying the
% method.

% T_centered_corr saves the spike times considered as stationary by the
% method. 

% nst_spike_times are the original spike times (that are elicited during
% Time_vector)
    
%% Statistics
    
stat = estadisticas(ISI_0, T_centered_corr);

% stat contains all the relevant statistics, nst_ISI, nst_Cv, nst_Lv and matrix_nst_ISIs
% are non-stationary. The statistics labeled with "st" at the beginning are the stationary
% ones.

%% Plot that shows the spike trains in both regimes

% Raster plot with rate profile

% Non-Stationary

subplot(2,1,1)

% Load time and rate variables

load('t.mat')
load('rate.mat')

% Convert rate for plotting purposes.

index = [];
index = find(rate == 50);

rate(index) = 2;

index = [];
index = find(rate == 5);

rate(index) = 1.1;

plot(t,rate,'r')
hold on

plot(nst_spike_times{1},1,'ok','MarkerFaceColor','k')
hold on

set(gca,'FontSize',24)
set(gca,'LineWidth',2)
set(gcf,'Color',[1 1 1])
axis([385 450 0.8 2])
box off

xlabel('Time (sec)','FontSize',24)
ylabel('neuron','FontSize',24)
title('Non-Stationary','Fontsize',24)
set(gca,'ytick',[])

% Stationary

subplot(2,1,2)

plot(t,rate,'r')
hold on

for i=1:size(T_centered_corr{1},2)
        
    plot(T_centered_corr{1}{i},1,'ok','MarkerFaceColor','k')
    hold on
    
end

set(gca,'FontSize',24)
set(gca,'LineWidth',2)
set(gcf,'Color',[1 1 1])
axis([385 450 0.8 2])
box off

xlabel('Time (sec)','FontSize',24)
ylabel('neuron','FontSize',24)
title(['Stationary - ' 'NumberISI = ' num2str(numberISI) ' Percentage = ' num2str(percentage)],'Fontsize',24)
set(gca,'ytick',[])