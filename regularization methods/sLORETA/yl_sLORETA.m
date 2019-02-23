%*************************************************************************
%   Copyright (C)  Bioengineering, UCLA
%   All Rights Reserved.
%
%   Created     :   2014-4-17 
%   Author      :   Ying Li (yingli@ucla.edu)
%   Advisor     :   Prof. Wentai Liu
%
%   File Name   :   yl_sLORETA.m
%   Abstract    :   calculate the inverse solution use sLORETA method

%   Inputs
%      'noisecov'          % can't be empty, otherwise will not use regularization
%      'sourcecov'         % can be empty
%      'LFM'               % lead field matrix
%      'snr'               % default = 5;
%      'dowhiten'          % default = 1;
%      'doscale'           % default = 1; scale R, so that trace(ARA')/trace(C) = 1

%*************************************************************************


function [Jy] = yl_sLORETA(K, U, s, V, pot, sel_regu, lambda)

b = pot;

if (nargin==6)   % when no lambda is specified
    
    if sel_regu == 1
        % regu method1: tikhonov + l-curve
        lambda_l = l_curve(U,s,b);
        [Jy T]  = yl_tikhonov(U,s,V,b,lambda_l);     % 3rd, seems ok
    elseif sel_regu ==2
        % regu method2: tsvd + l-curve
        k_l = l_curve(U,s,b,'tsvd');
        if isnan(k_l)
            Jy = zeros(size(K,2),1); % Spline Toolbox not available.
            fprintf('Spline Toolbox not available!');
            return;
        else
            [Jy T] = yl_tsvd(U,s,V,b,k_l);           % 4th, seems worst
        end
    elseif sel_regu == 3
        % regu method3: tikhonov + gcv
        lambda_gcv = gcv(U,s,b);
        [Jy T] = yl_tikhonov(U,s,V,b,lambda_gcv);   % 2nd, seems good
    elseif sel_regu ==4
        % regu method4: tsvd + gcv
        k_gcv = gcv(U,s,b,'tsvd');
        %Jy = tsvd(U,s,V,b,k_gcv);            % 1st, seems best?
        [Jy T] = yl_tsvd(U,s,V,b,k_gcv);
    end
    
elseif (nargin==7)   % when lambda is specified
    if (sel_regu == 1 | sel_regu == 3)
        lambda = 0;
        [Jy T] = yl_tikhonov(U,s,V,b,lambda);   % 2nd, seems good
    elseif (sel_regu == 2 | sel_regu == 4)
        k = size(U,1);
        [Jy T] = yl_tsvd(U,s,V,b,k);
    end
    
end
%standardized
for i=1:size(T,1)
    S(i) = T(i,:)*K(:,i);
    %T(i,:)=sqrt(1/S(i))*T(i,:);
    Jy(i)= Jy(i)/sqrt(S(i));
end

end