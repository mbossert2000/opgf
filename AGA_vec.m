function [f_curve, Re] = AGA_vec(D,roughness)
%AGA_vec: generates the look-up table of friction coefficients for all pipe in the network using the American Gas Association correlation formula
%   D: list of pipe diameters [m]
%   roughness: list of pipe roughness [m]
%   f_curve: sample points of the AGA friction coefficient function
%   Re: reyonld number of the sample points
n_pipes = length(D);
Re_min = 1e3;
Re_max = 1e8;
n_points = 100;
Re = logspace(3, 8, n_points);
f_curve = zeros(n_pipes,n_points);
Re_cr = 35.235 * (roughness./D).^-1.1039;
f0 = 1e-1; %order of magnitude
options = optimset('Display','off');
for pipe=1:n_pipes
    for i = 1:n_points
        if Re(i) <= Re_cr(pipe)
            aga_implicit = @(f) (1./sqrt(f)) + 2*log10(2.825./(Re(i).*sqrt(f)));
            f_curve(pipe,i) = fsolve(aga_implicit,f0,options);
        else
            f_curve(pipe,i) = (-2*log10(roughness(pipe)./(D(pipe)*3.7))).^-2;
        end
    end
end
f_curve = f_curve';
end

