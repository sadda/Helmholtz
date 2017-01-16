function options = SetField(options, fieldName, fieldValue)
    
    if ~isfield(options, fieldName)
        options.(fieldName) = fieldValue;
    end
end