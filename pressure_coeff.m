function [K_inj] = pressure_coeff(M,C,p)
%PRESSURE_COEFF Compute the pressure sensitivity coefficinets
%   M: incidence matrix of gas network
%   C: vector of gas pipe coefficients for weymouth equation
%   p: pressure of the nodes in the gas network
%   K_inj: matrix of pressure sentivity coefficients
[n_pipes,n_nodes] = size(M);

K_inj = zeros(n_nodes,n_nodes);

%build the matrix:
A = zeros(n_nodes,n_nodes);
for n=2:n_nodes %lines (dQ_n/dQ_i)
    %build equations:
    for pipe=1:n_pipes
        if M(pipe,n) == 1
            %pipe is departing
            from = n;
            sgn = 1;
            for to_=1:n_nodes
                if M(pipe,to_) == -1
                    to = to_;
                    if p(from) < p(to)
                        %switch from and to
                        temp = from;
                        from = to;
                        to = temp;
                        sgn = -1*sgn;
                    end
                    A(n,from) = A(n,from) + sgn*C(pipe)*p(from)/(sqrt(p(from)^2 - p(to)^2));
                    A(n,to) = A(n,to) - sgn*C(pipe)*p(to)/(sqrt(p(from)^2 - p(to)^2));
                end
            end
        elseif M(pipe,n) == -1
            %pipe is arriving
            to = n;
            sgn = -1;
            for from_=1:n_nodes
                if M(pipe,from_) == 1
                    from = from_;
                    if p(from) < p(to)
                        %switch from and to
                        temp = from;
                        from = to;
                        to = temp;
                        sgn = -1*sgn;
                    end
                    A(n,from) = A(n,from) + sgn*C(pipe)*p(from)/(sqrt(p(from)^2 - p(to)^2));
                    A(n,to) = A(n,to) - sgn*C(pipe)*p(to)/(sqrt(p(from)^2 - p(to)^2));
                end
            end
        end
    end
end
A = A(2:end,2:end); %disregard the slack

for m=2:n_nodes %injections
    b = zeros(n_nodes-1,1);
    b(m-1) = 1;
    
    K_inj(2:end,m) = A\b;
end
K_inj = K_inj(:,2:end);
end

