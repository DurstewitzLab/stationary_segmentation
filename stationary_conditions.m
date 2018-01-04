function [T_centered_corr, ISI_0, data_no_delay] = stationary_conditions(data,numberISI,percentage, vector)

%(c) 2014 Claudio S. Quiroga-Lombard, Bernstein Center for Computational
% Neuroscience Heidelberg-Mannheim

% Program that computes stationary conditions

T_centered_corr = [];
data2 = [];

%% Data Conversion

% Convert ST vector or matrix into cell array

if ~iscell(data)
    ST0=data; data=[];
    if size(ST0,2)==1, data{1}=ST0;
    else
        for i=1:size(ST0,2), data{i}=ST0(:,i); end;
    end;
end;

% Eliminates NaN from original data. The output is data2

data2 = [];

for i= 1:size(data,2)
    
    b = 1;
    
    for k= 1:length(data{i})
        
        c= isnan(data{i}(k));   
        
        if (c == 0 && data{i}(k) >0) % rat1 cell 9 has a zero! impossible
         
            data2{i}(b)= data{i}(k);
            
            b = b + 1;
                                    
        end 
        
    end
    
    if isempty(data{i})
        
        data2{i} = [];
               
    end
           
end

% Check dimensions

if size(data,2) ~= size(data2,2)
    
    error('myApp:argChk', 'Error de dimension');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This part of the code selects from all
% the spike times only the ones that are elicited during the defined
% recording period.

%% Define "recording times" boundaries

[data_no_delay, ISI_no_delay, ISI_0] = delay_division(data2, vector);

% Data_no_delay is a cell which containes all the spikes contained in
% the recording time.
% ISI_no_delay are the ISIs of data_no_delay.
% This part of the code can be improved since it was designed originally
% for more general purposes.
 
% Error message

if size(data2,2) ~= size(data_no_delay,1)
    
 error('myApp:argChk', 'Wrong_1')
        
end    
    
neurons = size(data_no_delay,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Maximum likelihood to find the best transformation parameter

int = 0.05;
ini = 0.05;
final = 1;

[ISI_NDT,ISI_2, L_max, best_lambda, exp] = Lmax(ISI_no_delay,ISI_0,int,ini,final);

clear int
clear ini
clear final
 
% ISI_2 is the best approximation to a Gaussian distribution with zero mean and std = 1
% In this situation ISI_NDT has the same information that ISI_2. 

% Compute the ISI for sliding window(numberISI)

g_ISI = cell([],[]);

t_ISI =  cell([],[]);

chi2_ISI = cell([],[]);

T = cell([],[]);

T_centered = cell([],[]);

for k=1:neurons   
      
    a = [];  
    a = sprintf('matrix_form_%d.txt',k);
            
    for l=1:size(ISI_NDT,2)        
        
        for c=1:length(ISI_NDT{k,l})-numberISI+1
            
            % An error can appear here if the dimension of ISI_NDT is less
            % than 9

            T{k,l}(c) = data_no_delay{k,l}(c);

            T_centered{k,l}(c) = data_no_delay{k,l}(c) + (data_no_delay{k,l}(c+numberISI) - data_no_delay{k,l}(c))/2; 

            % Final final spike time - initial initial spike time of the 10
            % ISIs == 11 initial spike times

            g_ISI{k,l}{c} = ISI_NDT{k,l}(c:numberISI-1+c);
            
            % Saves the info of all the transformed ISIs per sliding window
            % As always "k" denotes de neuron, "l" the "no delay interval"
            % and "c" the sliding window.
            
            % Compute the t_value and the chi2 over each sliding window
            
            t_ISI{k,l}(c) = mean(ISI_NDT{k,l}(c:numberISI-1+c))*sqrt(numberISI);
            
            chi2_ISI{k,l}(c) = sum(ISI_NDT{k,l}(c:numberISI-1+c).^2);
            
            % Remember that the process comes from ISI_2 which has zero
            % mean and std = 1.
            
            % Saves these values "in matrix form"

%             fid = fopen(num2str(a), 'a');
%             fprintf(fid, '%6.8f %6.8f\n', t_ISI{k,l}(c), chi2_ISI{k,l}(c));
%             fclose(fid)           
                                
        end
    
    end

end

% Theoretical Threshold for T_k Gaussian (first condition)

Thres_mean_down = norminv((1-percentage)/2,0,1);
Thres_mean_up = norminv(1-(1-percentage)/2,0,1);

% Theoretical Threshold for Chi2 (second condition)

Thres_chi2_down = chi2inv((1-percentage)/2,numberISI);
Thres_chi2_up = chi2inv(1-(1-percentage)/2,numberISI);
    
%% Apply stationarity

B = cell([],[]);

T_centered_corr = cell([],[]);

for k=1: neurons
                
        for l=1:size(t_ISI,2) 

            index_1 = [];
            
            index_1 = find(t_ISI{k,l} <= Thres_mean_up & t_ISI{k,l} >= Thres_mean_down); % First Condition
            
            index_2 = [];
            
            index_2 = find(chi2_ISI{k,l} <= Thres_chi2_up & chi2_ISI{k,l} >= Thres_chi2_down); % Second Condition
            
            int_index = [];
            
            int_index = intersect(index_1,index_2); % ISIs which meet both conditions
            
            if ~isempty(int_index)
                
                % It only makes sense to define stationary intervals if
                % there is at least one ISI that is stationary.
            
                % Divide in ST_intervals (jump in int_index)

                jumps = [];

                jumps = [0 find(diff(int_index) >1) length(int_index)];

                ccc = 1;

                for cc=1:length(jumps)-1

                    B{k,l}{cc} = [t_ISI{k,l}(int_index(1+jumps(cc):jumps(cc+1)))', chi2_ISI{k,l}(int_index(1+jumps(cc):jumps(cc+1)))', T_centered{k,l}(int_index(1+jumps(cc):jumps(cc+1)))'];

                % As always, "k" defines the neuron.
                % Each column of B{k,l} contains the values (t,chi2) for each
                % of the sliding windows wich are Stationary and the third
                % column is the "centered stationary time".

                % For each k and each l I defined the stationary spike times as
                % those which lie between the defined stationary intervals
                % (third column of B)

                index = [];

                index = find(data_no_delay{k,l} >= B{k,l}{cc}(1,3) & data_no_delay{k,l} <= B{k,l}{cc}(end,3));

                    if ~isempty(index)

                        T_centered_corr{k,l}{ccc} = data_no_delay{k,l}(index);

                        ccc = ccc + 1;

                    end

                end
                
            end
            
            
        end
             
end    

