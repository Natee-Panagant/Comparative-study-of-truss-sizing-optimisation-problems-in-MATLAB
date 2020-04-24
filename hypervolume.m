function [hv] = hypervolume(f,W)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of Hypervolume (HV) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% hv = hypervolume (higher is better)
% f = Objectives function vector [f1;f2;...,fn]
% W = Reference points [w1;w2;...,wn]
%
if numel(f)==0
    hv=0;
else
    [~,ns]=sort(f(1,:));% Sort the first objective value in ascending order
    f=f(:,ns);
    hv=abs(W(1)-f(1,1))*abs(W(2)-f(2,1))+sum(abs(W(1)-f(1,2:end)).*abs(f(2,1:end-1)-f(2,2:end)));
end
end

