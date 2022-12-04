function m = ideal_gas_mass(P,V,T,r)
%
% Calculate mass from ideal gas equation
%
%   ideal_gas_mass(P,V,T,r)
%
%   P: pressure in Pa
%   V: volume in m3
%   T: temperature in K
%   r: individual gas constant J/KgK
    m = P * V / (r * T);
end