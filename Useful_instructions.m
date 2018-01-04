
% USEFUL INSTRUCTIONS

% 1. The cell L is used to load the data. The data should be in the same format as
% for the STMtx matrix, where each neuron is a column and the spike times are in
% rows.

% 2. Time vector refers to the recording times of the experiment.

% 3. NumberISI = 10 and percentage = 0.96 are the parameters of the method.
% NumberISI is the number of ISIs consider in each sliding window.
% Percentage [0 1] are the confidence limits for the Gaussian and X^2
% distributions. Lower NumberISI and/or lower percentage increase the
% constraints to the spike times.

% 4. The cell T_centered_corr is the one that has the stationary (ST) spike times.
% T_centered_corr{i}{j} gives the ST spike times of neuron "i" in the ST
% segment "j". The ST spike times of a given neuron are the collection of all the 
% spike times ocurring in the {j} ST segments.
