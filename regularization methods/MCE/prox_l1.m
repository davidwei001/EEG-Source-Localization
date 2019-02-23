
function [prox_v] = prox_l1(v, beta)

% Proximity operator of the following form:
% prox(v) = argmin{ 1/2*||u-v||^2 + beta*||u||^2}
% by Ying Li, 10/17/14

prox_v = sign(v).* max( abs(v)-beta, 0 );

% prox_v = v * ( 1- min( 1, beta/abs(v)) );   % equivalent to the above

