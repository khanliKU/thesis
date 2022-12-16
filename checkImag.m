function result = checkImag(val)
    if imag(val) == 0 && ~isnan(val)
        result = val;
    else
        result = val;
        %throw(MException(999,'Imaginary Value!'))
    end
end