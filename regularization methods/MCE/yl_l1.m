
function [u, z, zt, history] = yl_l1(A, b, pm, ss, DISPLAY_FLAG)

% Calculate the inverse solution using TV regularization:
% min 1/2*||A*u-b||^2 + alpha2*||z||^1
% s.t. u=z;

% The augmented Lagrangian function:
% L(u,z,zt)
%    = 1/2*||A*u-b||^2 + alpha2*||z||^1
%      + rho/2*(||u-z+zt||^2)

% Input:
%       A: lead field matrix
%       b: measurement of the sensors
%       V: gradient matix
%       rho, alpha1, alpha2: regularization parameter
%       DISPLAY_FLAG: 1-display the iteration info; 0-don't display

% Output:
%       u: estimated current denisty

% Assumption: dipoles are perpendicular to the cortex surface
% by Ying Li, 11/03/14

%% Initialization

%tol = 10^(-12);   % 
% tol_rel = 10^(-4); % relative tolerance
% MAX_ITERS = 10^4;  % can be changed in the main function,500 is too small; maybe good for TGV, but not enough for TV
% rho = 10^(-11);
% alpha2 = 5*rho;

n = size(A,1);    % number of sensors
m = size(A,2);    % number of dipoles
   
if isfield(pm,'tol_rel'); tol_rel = pm.tol_rel; end
if isfield(pm,'MAX_ITERS'); MAX_ITERS = pm.MAX_ITERS; end
if isfield(pm,'rho'); rho = pm.rho; end
if isfield(pm,'alpha2'); alpha2 = pm.alpha2; end

if isfield(pm,'Lp')
    Lp = pm.Lp;
else
    P = A'*A  + rho*( eye(m) );   % not sparse at all. Out of memory
    Lp = chol(P,'lower');   % much faster than pinv. not sparse
    %[Lp, p]=chol(P);
end

u = zeros(m,1);

z = zeros(m,1);
zt = zeros(m,1);

%% Iteration
iter = 0;
uchange = [];

while(1)
    
    iter = iter+1;
    u0 = u;    % old u
    u = Lp'\ (Lp \ ( A'*b + rho*(z-zt) ));  % much faster than P \ ( )
    %u = P'\ (P \ ( A'*b + rho*(z-zt) ));
    z = prox_l1 (u + zt, alpha2/rho );
    zt = zt + (u - z);
    
    %Err_nor = norm(u-u0)/norm(u0);
    Err = norm(u-u0)/norm(u0);   % relative error
    if (DISPLAY_FLAG) disp(['iters: ',num2str(iter),'   u change: ',num2str(Err)]); end;
    %[history.Error(iter) history.Error_rel(iter) history.dis_maxi(iter) history.df(iter) history.Half_radius(iter) history.Half_area(iter) history.sparsity(iter) history.SD(iter) history.ROC(iter) history.AUC_unbiased(iter)] = yl_criteria(ori ,u, ss,0); 
    
    % diagnostics, reporting, termination checks
    history.term1(iter)  = 1/2*sum((A*u - b).^2);
    history.term3(iter)  = alpha2*norm(u,1);
    %history.term2_noalpha(iter)  = norm(V*u,1);
    %history.term3_noalpha(iter)  = norm(u,1);
    history.objective(iter)  =  history.term1(iter)+ history.term3(iter);
    history.relerr(iter) = Err;
    
    if (Err < tol_rel) || (iter > MAX_ITERS)    % stopping criteria
        disp(['TV:  u change: ',num2str(Err),'   Iteration: ',num2str(iter)]); 
        %save('TV_uchange','uchange');
        %filename = strcat('TV_iter',num2str(iter,'%06d'));
        %save(filename,'u')
        break;
    end
    
%     if mod(iter,100)==0
%         uchange = [uchange Err];
%         filename = strcat('TV_iter',num2str(iter,'%06d'));
%         save(filename,'u')
%     end
    
end

end


