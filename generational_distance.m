function gd=generational_distance(fpareto,RefFront)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of General Distance (GD) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% gd = generational distance (lower is better)
% fpareto = obtained Pareto front
% RefFront = reference Pareto front
%
if numel(fpareto)==0
   gd=1e10;
   return;
end
N=size(fpareto,2);
gd_2=zeros(1,N); %gd^2
for i=1:N
    d_2=(fpareto(:,i)-RefFront(:,:)).^2;
    d_2=sum(d_2,1);
    gd_2(i)=min(d_2);
end
gd=sqrt(sum(gd_2))/N;
