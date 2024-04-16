function [f] = AGA(Re_vec,f_curve,Re)
%AGA interpolate the proper friction coefficient for a given Reyonld
%number. Need to be used in combinationwith AGA_vec which generates the
%lookup table of friction coefficient function for the whole gas network
%   Re_vec: Reynold number of sample points in the look-up table f-curve
%   f_curve: lookup table of friction coefficient values for the reynold
%   numbers in Re_vec
%   Re: reyonld number to be interpolated

%   f: friction coefficient for Re
f = zeros(length(Re),1);
for k=1:length(Re)
    f(k) = interp1(Re_vec,f_curve(:,k),Re(k));
end
end

