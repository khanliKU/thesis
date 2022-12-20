function result = checkImag(val)
    if imag(val) == 0 && ~isnan(val)
        result = val;
    else
        throw(MException(999,'Imaginary Value!'))
    end
end