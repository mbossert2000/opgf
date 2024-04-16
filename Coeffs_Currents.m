function [KI_p, KI_q, I0]=Coeffs_Currents(YYL,YYT,E,COEFFcomplex_alph,COEFFqcomplex_alph,nph,lines)
% returns current sensitivity coefficient in matrix form
% returns vector of current in each branch, seen from the two side
% !!!! YYT = A'*Yl*A
    Currents = cell(size(YYL)/nph);
    Icoeff = cell(1,size(COEFFcomplex_alph,2));
    Icoeffq = cell(1,size(COEFFcomplex_alph,2));
    Icoeff_complex = cell(1,size(COEFFcomplex_alph,2));
    Icoeffq_complex = cell(1,size(COEFFcomplex_alph,2));
    
    % For all active/reactive power injections
    for l = 1:size(COEFFcomplex_alph,2)
        % Re-initialize Temp Matrices
        temp_p = cell(size(YYL)/nph);
        temp_q = cell(size(YYL)/nph); 
        temp_p_complex = cell(size(YYL)/nph);
        temp_q_complex = cell(size(YYL)/nph);
        
        % For all rows of Admittance Matrix
        for i = 1:nph:size(YYL,1)
            for j = 1:nph:size(YYL,2)
                
                % Compute Currents only once
                if (l == 1)
                    Currents{ceil(i/nph),ceil(j/nph)} = YYL(i:i+(nph-1),j:j+(nph-1))*(E(i:i+(nph-1)) - E(j:j+(nph-1))) + ...
                                    YYT(i:i+(nph-1),j:j+(nph-1))*E(i:i+(nph-1));
                end
                
                if( max(abs(Currents{ceil(i/nph),ceil(j/nph)})) ~= 0)
                    % Compute Complex Coeffs
                    dIij_dPl = YYL(i:i+(nph-1),j:j+(nph-1))*(COEFFcomplex_alph(i:i+(nph-1),l) - COEFFcomplex_alph(j:j+(nph-1),l)) + ...
                               YYT(i:i+(nph-1),j:j+(nph-1))*COEFFcomplex_alph(i:i+(nph-1),l);
                    dIij_dQl = YYL(i:i+(nph-1),j:j+(nph-1))*(COEFFqcomplex_alph(i:i+(nph-1),l) - COEFFqcomplex_alph(j:j+(nph-1),l)) + ...
                               YYT(i:i+(nph-1),j:j+(nph-1))*COEFFqcomplex_alph(i:i+(nph-1),l);

                    % Compute Coeffs
                    temp_p{ceil(i/nph),ceil(j/nph)} =  ((1./abs(Currents{ceil(i/nph),ceil(j/nph)})).*real(conj(Currents{ceil(i/nph),ceil(j/nph)}).*dIij_dPl));
                    temp_q{ceil(i/nph),ceil(j/nph)} =  ((1./abs(Currents{ceil(i/nph),ceil(j/nph)})).*real(conj(Currents{ceil(i/nph),ceil(j/nph)}).*dIij_dQl));
                    temp_p_complex{ceil(i/nph),ceil(j/nph)} =  dIij_dPl;
                    temp_q_complex{ceil(i/nph),ceil(j/nph)} =  dIij_dQl;
                else
                    temp_p{ceil(i/nph),ceil(j/nph)} =  (zeros(nph,1));
                    temp_q{ceil(i/nph),ceil(j/nph)} =  (zeros(nph,1));
                    temp_p_complex{ceil(i/nph),ceil(j/nph)} =  (zeros(nph,1));
                    temp_q_complex{ceil(i/nph),ceil(j/nph)} =  (zeros(nph,1));
                end
                
            end
        end
        
        Icoeff{l} = temp_p;
        Icoeffq{l} = temp_q;
        Icoeff_complex{l} = temp_p_complex;
        Icoeffq_complex{l} = temp_q_complex; 
    end
    
    n_buses = size(YYT,1);
    n_lines = height(lines);
    I0 = zeros(2*n_lines,1);
    KI_p = zeros(2*n_lines,n_buses-1);
    KI_q = zeros(2*n_lines,n_buses-1);
    for l=1:n_lines
        from = lines(l,1);
        to = lines(l,2);
        I0(2*l-1) = abs(Currents{from,to});
        I0(2*l) = abs(Currents{to,from});
        for b=1:n_buses-1
            Icoeff_b = Icoeff{b};
            KI_p(2*l-1,b) = Icoeff_b{from,to};
            KI_p(2*l,b) = Icoeff_b{to,from};
            Icoeffq_b = Icoeffq{b};
            KI_q(2*l-1,b) = Icoeffq_b{from,to};
            KI_q(2*l,b) = Icoeffq_b{to,from};
        end
    end
end