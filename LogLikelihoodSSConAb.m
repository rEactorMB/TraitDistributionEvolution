function log_sum = LogLikelihoodSSConAb(trait_values, a, D, L, x0, t_values, varargin)
% Log likelihood calculation for data whose probability is given by the solution
% to the partial differential equation dy/dt = -a dy/dx + D/2 d^2y/dx^2, (y=y(x,t))
% subject to a delta function initial condition: y(x,0) = delta(x - x0)
% and Dirichlet boundaries at x=0,L: y(0,t) = y(L,t) = 0
% trait_values represent x and t_values represent t
if nargin == 6
    weight = 100;
elseif nargin == 7
    weight = cell2mat(varargin(1));
end
log_sum = 0;
lt = length(t_values);
if lt == 1
    if size(trait_values, 1) > size(trait_values, 2)
        trait_values = trait_values.';
    end
    log_sum = sum(log(abs(SSConAb(trait_values(logical((trait_values > 0) .* (trait_values < L))), a, D, L, x0, t_values))));
    log_sum = log_sum - (weight * nnz((trait_values <= 0) + (trait_values >= L)));
elseif size(trait_values, 2) == lt
    domain = (trait_values > 0) .* (trait_values < L);
    penalty = weight * (lt - nnz(domain));
    t_values = t_values(domain);
    trait_values = trait_values(domain);
    for i = 1:nnz(domain)
        log_sum = log_sum + log(abs(SSConAb(trait_values(i), a, D, L, x0, t_values(i))));
    end
    log_sum = log_sum - penalty;
elseif size(trait_values, 1) == lt
    if ismatrix(trait_values)
        for i = 1:lt
            log_sum = log_sum + sum(log(abs(SSConAb(trait_values(i, logical((trait_values(i, :) > 0) .* (trait_values(i, :) < L))), a, D, L, x0, t_values(i)))));
            log_sum = log_sum - (weight * nnz((trait_values(i, :) <= 0) + (trait_values(i, :) >= L)));
        end
    else
        log_sum = LogLikelihoodSSConAb(trait_values.', a, D, L, x0, t_values, weight);
    end
end
end
