function TableToTex(data, formatAll)
    for i=1:size(data,1) 
        for j=1:size(data,2)
            fprintf(formatAll{i}, data(i,j)) ;            
            if j < size(data,2)
                fprintf(' & ');
            else
                fprintf(' \\\\\n');
            end
        end
    end
end