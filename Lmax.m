function [ISI_NDT,ISI_2, L_max, best_lambda, exp1] = Lmax(ISI_no_delay,ISI_0,int,ini,final)

% Box and Cox - An Analysis of Transformations 1964
%(c) 2014 Claudio S. Quiroga-Lombard, Bernstein Center for Computational
% Neuroscience Heidelberg-Mannheim

% Pick the exponent for the nonlinear transformation into approx Gaussian distribution
% Non linear transformation to make ISIs more Gaussian distributed

exp = ini:int:final;

exp1 = exp;

lexp = length(exp);

ISI_w=[];

neurons = size(ISI_no_delay,1);

for k=1:neurons
    
    %index_conserved = [];
    
    %index_conserved = find(k == conserved); % selecccion las neuronas con las que se trabaja
    
    %if ~isempty(index_conserved)  % Si la neurona pertenece a la conservada
        
        for w= 1:lexp

            ISI_w{k,w} = ((ISI_0{k}).^exp(w) - 1)./exp(w);

            dim=0;

        end

        % Include log in the last point

        if (w == lexp)

            ISI_w{k,w+1} = log(ISI_0{k});

            dim=1; %considering log adds a virtual dimensio to the exponents

        end
        
    %end
           
end

% Maximum likelihood

L_max = [];

for k=1:neurons

    %index_conserved = [];
    
    %index_conserved = find(k == conserved); % selecccion las neuronas con las que se trabaja
    
    %if ~isempty(index_conserved)  % Si la neurona pertenece a la conservada
        
        for i=1:size(ISI_w,2)

            n = length(ISI_w{k,i});

            variance = var(ISI_w{k,i});

            if i < size(ISI_w,2)

                L_max(k,i) =  -(n/2)*log(variance) + (exp(i) -1)*sum(log(ISI_0{k}));

           else
                L_max(k,i) = -(n/2)*log(variance) + (-1)*sum(log(ISI_0{k}));

           end

            clear n
            clear variance

        end
                
    %end
    
end
% 
% 
% k=2;
% 
% plot([exp exp(end)+0.1],L_max(k,:),'.-k')
% plot(exp,L_max(k,1:end-1),'.-k')
% hold on
% plot(exp(end)+0.1,L_max(k,end),'ok')
% xlabel('l','FontName','Symbol','FontSize',24)
% ylabel('Lmax','FontName','Arial','FontSize',24)

% Chooose the maximum and standardize

ISI_2_pre = [];
ISI_2 = [];
best_lambda = [];

exp = [exp 0];

for k = 1:neurons

    %index_conserved = [];
    
    %index_conserved = find(k == conserved); % selecccion las neuronas con las que se trabaja
    
    %if ~isempty(index_conserved)  % Si la neurona pertenece a la conservada
        
        index =  find(L_max(k,:) == max(L_max(k,:)));

        best_lambda(k) = exp(index);

        ISI_pre{k} = ISI_w{k,index(1)};

        ISI_2{k} = (ISI_pre{k} - mean(ISI_pre{k}))./std(ISI_pre{k});

        clear index
               
    %end
    
end

% Transfer data to data_no_delay

a = [];
a = cellfun(@length, ISI_no_delay); 

ISI_NDT = [];

for k = 1: size(a,1)
    
    c = 0;
    
    for l=1:size(a,2)
        
        ISI_NDT{k,l} = ISI_2{k}(1+c:a(k,l)+c);
        
        c = c + a(k,l);
        
    end
    
end

b = [];
b = cellfun(@length, ISI_NDT); 

if ~isequal(b, a)

   error('myApp:argChk', 'Wrong_1')
        
end
