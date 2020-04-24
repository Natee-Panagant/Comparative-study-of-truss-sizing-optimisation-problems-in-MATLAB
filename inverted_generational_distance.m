function igd=inverted_generational_distance(fpareto,RefFront)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of Inverted General Distance (IGD) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% igd = inverted generational distance (lower is better)
% fpareto = obtained Pareto front
% RefFront = reference Pareto front
%
if numel(fpareto)==0
   igd=1e10;
   return;
end
N=size(RefFront,2);
igd_2=zeros(1,N); %igd^2
for i=1:N
    d_2=(RefFront(:,i)-fpareto(:,:)).^2;
    d_2=sum(d_2,1);
    igd_2(i)=min(d_2);
end
igd=sqrt(sum(igd_2))/N;