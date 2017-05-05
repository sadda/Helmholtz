function TableToTex(data, formatAll)
    for i=1:size(data,1)
        for j=1:size(data,2)
            value      = data(i,j);
            formatType = formatAll{i};
            switch formatType(end)
                case 'e'
                    if value == 0
                        error('I should really do something about this');
                    end
                    orderOfMag         = floor(log10(abs(value)));
                    valueToPrint       = value*10.^(-orderOfMag);
                    formatTypeMod      = formatType;
                    formatTypeMod(end) = 'f';

                    fprintf('$');
                    fprintf(formatTypeMod, valueToPrint);
                    fprintf('\\sdot 10^{%d}', orderOfMag);
                    fprintf('$');
                case 's'
                    error('Not implemented. Need to pass to cell arrays.');
                    % fprintf(formatAll{i}, value);
                otherwise
                    fprintf('$');
                    fprintf(formatAll{i}, value);
                    fprintf('$');
            end
            if j < size(data,2)
                fprintf(' & ');
            else
                fprintf(' \\\\\n');
            end
        end
    end
end