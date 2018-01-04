function [data_no_delay, ISI_no_delay, ISI_0] = delay_division(data2, vector)

%(c) 2014 Claudio S. Quiroga-Lombard, Bernstein Center for Computational
% Neuroscience Heidelberg-Mannheim

data_no_delay = [];
        
for l = 1:size(vector,1) 
                
        for s = 1: size(data2,2)
                  
                % I will only conserve the spike times which lie between
                % the limits of the vector.
            
                index = [];
                
                index = find(data2{s} >= vector(l,1) & data2{s} <= vector(l,2));
                                                  
                data_no_delay{s,l} = data2{s}(index);
               
        end
        
end
                    
if size(data2,2) ~= size(data_no_delay,1)
            
   error('myApp:argChk', 'Wrong!!!');
           
end

%% Compute the ISIs per segment

a = [];
a = cellfun(@length, data_no_delay);

ISI_no_delay = [];
ISI_0 =cell([],[]);

for k = 1:size(a,1)
    
    b = sprintf('ISI_0_acumul_%d.txt',k);
        
    for l=1:size(a,2)
                    
        ISI_no_delay{k,l} = diff(data_no_delay{k,l});
                    
        fid = fopen(num2str(b), 'a');
        fprintf(fid, '%6.8f\n', ISI_no_delay{k,l});
        fclose(fid)           
        
    end       
    
    ISI = [];
    
    fid = fopen(num2str(b));
    ISI = fscanf(fid, '%g', [1 inf]);
    fclose(fid);
          
    ISI_0{k} = ISI;
    
    delete(num2str(b));

end

b = [];
b = cellfun(@length, ISI_no_delay); 

for ii=1:size(a,1)

    index = [];
    
    index = find(a(ii,:) == 0);
    
    if sum(a(ii,:))-size(a,2)+length(index)~= sum(b(ii,:))

       error('myApp:argChk', 'Wrong_1')

    end
    
end

c = [];
c = cellfun(@length, ISI_0);

if ~isequal(sum(b,2), c')

   error('myApp:argChk', 'Wrong_2')
        
end
