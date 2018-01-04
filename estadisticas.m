function stat = estadisticas(ISI_0, T_centered_corr)

%(c) 2014 Claudio S. Quiroga-Lombard, Bernstein Center for Computational
% Neuroscience Heidelberg-Mannheim

% T_centered_corr is the most important function of
% stationary_conditions_pop_max_like.m, this is because it contains the
% stationary "spike times".
% T_centered_corr is of the form T{k,j}{i}, k is the neuron, j is the stationary interval
% and i is the cell where the spike times are saved.

% Statistics for the NST case

stat = [];

for k=1: size(ISI_0,2)
    
    stat.nst_ISI{k} = ISI_0{k};
    
    stat.nst_Cv(k) = std(stat.nst_ISI{k})/mean(stat.nst_ISI{k});
    
    x1= stat.nst_ISI{k}(1:end-1) - stat.nst_ISI{k}(2:end);
    x2= stat.nst_ISI{k}(1:end-1) + stat.nst_ISI{k}(2:end);

    stat.nst_Lv(k) = sum((3*x1.^2)./x2.^2)/(length(stat.nst_ISI{k})-1);
    
    clear x1
    clear x2
    
end

% Statistics for the ST case

% The first step is to compute the ISIs for each of the cell entries.
% and get all of them together to form a vector

a = [];
a = cellfun(@length,T_centered_corr);

for k=1: size(T_centered_corr,1)
    
    b = sprintf('st_ISI_%d.txt',k);
    
    fid = fopen(num2str(b), 'a');
    fprintf(fid, '%6.6f\n',[]);
    fclose(fid)
    
    c = [];
    
    for l=1: size(T_centered_corr,2)
    
        for i=1: a(k,l)
            
            ISI = diff(T_centered_corr{k,l}{i});
        
            fid = fopen(num2str(b), 'a');
            fprintf(fid, '%6.6f\n', ISI);
            fclose(fid)
            
            c(l,i) = size(T_centered_corr{k,l}{i},2);
            
        end
        
    end
    
    % Check for dimension      
    
    fid = fopen(num2str(b));
    stat.st_ISI{k} = fscanf(fid, '%g', [1 inf]);
    fclose(fid);

    delete(num2str(b))
    
    index = [];
    
    index = find( c(1:end) > 0);
    
    if size(stat.st_ISI{k},2) ~= sum(c(1:end)) - length(index)
        
          error('myApp:argChk', 'Wrong_1!!!');
           
    end
    
    % Compute the statistics
        
    stat.st_Cv(k) = std(stat.st_ISI{k})/mean(stat.st_ISI{k});

    x1= stat.st_ISI{k}(1:end-1) - stat.st_ISI{k}(2:end);
    x2= stat.st_ISI{k}(1:end-1) + stat.st_ISI{k}(2:end);

    stat.st_Lv(k) = sum((3*x1.^2)./x2.^2)/(length(stat.st_ISI{k})-1);

    clear x1
    clear x2        
            
end

% Convert cells into matrix to save the ISIs in a correct way

a = [];
a = cellfun(@length,stat.nst_ISI);  
    
matrix_nst_ISIs = NaN(max(a) , size(stat.nst_ISI,2));

b = [];
b = cellfun(@length,stat.st_ISI);

matrix_st_ISIs = NaN(max(b) , size(stat.st_ISI,2));

for k=1:size(stat.nst_ISI,2)
    
    matrix_nst_ISIs(1:size(stat.nst_ISI{k},2),k) = stat.nst_ISI{k};
    
    matrix_st_ISIs(1:size(stat.st_ISI{k},2),k) = stat.st_ISI{k};
    
end

index1 = [];

index1 =  find(~isnan(matrix_nst_ISIs(1:end)));

index2 = [];

index2 = find(~isnan(matrix_st_ISIs(1:end)));

if length(index1) ~= sum(a) || length(index2) ~= sum(b)  

    error('myApp:argChk', 'Wrong_1!!!');
           
end

stat.matrix_nst_ISIs = matrix_nst_ISIs;
stat.matrix_st_ISIs = matrix_st_ISIs;
    
    
    
        