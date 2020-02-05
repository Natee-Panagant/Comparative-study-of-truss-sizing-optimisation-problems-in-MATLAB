function ste = spacing_to_extent(fpareto)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of Spacing-to-Extent (STE) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ste = spacing-to-extent (lower is better)
% fpareto = obtained Pareto front
%
N=size(fpareto,2);
d=zeros(1,N);
% Front spacing
if N > 2
    for i=1:N
        di=abs(fpareto(:,i)-fpareto(:,[1:(i-1) (i+1):N]));
        di=sum(di,1);
        d(i)=min(di);
    end
    dmean=mean(d);
    spacing=sqrt(sum((d-dmean).^2)/(N-1));
else
    spacing=1e10;
end

% Front extent
fmax=max(fpareto,[],2);
fmin=min(fpareto,[],2);
extent=max(sum(fmax-fmin),1e-10);

% STE
ste=spacing/extent;

