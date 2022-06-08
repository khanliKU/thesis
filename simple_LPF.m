function filtered = simple_LPF(input,fac)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    filtered(1) = input(1);
    for i=2:length(input)
        filtered(i) = filtered(i-1) + fac * (input(i) - filtered(i-1));
    end
end