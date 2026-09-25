function [coeff,condition]=allocate_directional_record(prior,observed)
% prior: frequency x direction complex amplitudes at the probe.
% Preserves each prior relative directional complex shape; matches observed sum.
observed=observed(:);assert(size(prior,1)==numel(observed));
total=sum(prior,2);mass=sum(abs(prior),2);
assert(all(abs(total)>64*eps*mass & mass>0), ...
    'Initial directional amplitudes cancel to numerical precision; no division floor is applied.');
W=prior./total;coeff=W.*observed;
condition=struct('per_frequency',sum(abs(W),2), ...
    'max',max(sum(abs(W),2)), ...
    'observed_energy_weighted',sum(abs(observed).^2.*sum(abs(W),2))/sum(abs(observed).^2), ...
    'sum_relative_error',norm(sum(coeff,2)-observed)/max(norm(observed),realmin));
assert(condition.sum_relative_error<1e-11);
end
