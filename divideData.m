function [data, labels, clust] = divideData(input_data)

% [data,labels,clust]=divide_data(input_data)
% divide_data computes the number of clusters in a dataset and identify the
% subgroups of the data
% If you use this code in your work, please cite:
%
% Köktürk, Ba?ak Esin, and Bilge Karaçal?. "Annealing-based model-free expectation maximisation for
% multi-colour flow cytometry data clustering." International Journal of Data Mining and Bioinformatics 14,
% no. 1 (2016): 86-99.


tempA = 0;

input_data_matris{1} = input_data;
ind_vector = [1:1:length(input_data)];

tempA = tempA+1;

input_data = input_data_matris{tempA};


sprintf('Division: root1')
[root1] = abmfem(input_data, ind_vector);
sprintf('Division: root1.data1')
[root11] = abmfem(root1.data1, root1.ind1);
sprintf('Division: root1.data2')
[root12] = abmfem(root1.data2, root1.ind2);


division = [1];
level = [1];
tempB = 0;

iteration = 0;
while(sum(division)~=0)
    iteration = iteration + 1;
    
    if( iteration == 1 )
        mothername = 'root1';
        division = [0];
    else
        tempC = find(division == 1);
        mothername=sprintf('root%d',level(tempC(1))) ;
        tempD = find(level==level(tempC(1)));
        division(tempD) = 0;
    end
    
    for i = 1:2
        childname = sprintf('%s%d', mothername,i);
        error1 = eval(sprintf('%s.error', childname));
        error = eval(sprintf('%s.error', mothername));
        
        if( error1 <= error+0.01 )
            parent1 = eval(sprintf('%s.data1', childname));
            ind1 = eval(sprintf('%s.ind1', childname));
            sprintf('Division: %s.data1', childname)
            [x] = abmfem(parent1, ind1);
            eval(sprintf('root%d1 = x', str2num(childname(5:end))));
            level = [level str2num(childname(5:end))];
            if( x.error <= error1+0.01 && (length(x.data1)+length(x.data2)) > 300 )
                division=[division 1];
            else
                division = [division 0];
                tempB = tempB + 1;
                clust{tempB} = parent1;
                label{tempB} = ind1;
            end
            
            parent2 = eval(sprintf('%s.data2', childname));
            ind2 = eval(sprintf('%s.ind2', childname));
            sprintf('Division: %s.data2', childname)
            [x] = abmfem(parent2, ind2);
            
            eval(sprintf('root%d2 = x', str2num(childname(5:end))));
            
            if( x.error <= error1+0.01 && (length(x.data1)+length(x.data2)) > 300 )
                division=[division 1];
            else
                division = [division 0];
                tempB = tempB+1;
                clust{tempB} = parent2;
                label{tempB} = ind2;
            end
            level = [level str2num(childname(5:end))];
            
        else
            xx = find(level == str2num(sprintf('%s%d', mothername(5:end), i)));
            tempE = sprintf('%s%d', mothername(5:end), i);
            if( ~isempty(xx) == 0 && length(tempE) < 3 )
                tempB = tempB + 1;
                clust{tempB} = eval(sprintf('%s.data%d', mothername, i));
                label{tempB} = eval(sprintf('%s.ind%d', mothername, i));
            end
        end
    end
end
data = [];
labels = [];

for k = 1:length(clust)
    data(label{k},:) = clust{k} ;
    labels(label{k}) = k;
end
end
