function rst = MOMVO(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)
%______________________________________________________________________________________
%  Multi-Objective Multi-Verse Optimization (MOMVO) algorithm source codes version 1.0
%
%  Developed in MATLAB R2016a
%
%  Author and programmer: Seyedali Mirjalili
%
%         e-Mail: ali.mirjalili@gmail.com
%                 seyedali.mirjalili@griffithuni.edu.au
%
%       Homepage: http://www.alimirjalili.com
%
%   Main paper:
%   S. Mirjalili, P. Jangir, S. Z. Mirjalili, S. Saremi, and I. N. Trivedi
%   Optimization of problems with multiple objectives using the multi-verse optimization algorithm, 
%   Knowledge-based Systems, 2017, DOI: http://dx.doi.org/10.1016/j.knosys.2017.07.018
%______________________________________________________________________________________

% function [Best_universe_Inflation_rate,Best_universe,Archive_F]=MOMVO(Max_time,N,ArchiveMaxSize)
rand('state',sum(100*clock));
tic 

Max_time=nloop;
N=nsol;
ArchiveMaxSize=narchive;

fobj=fun;
    dim=nvar;
    lb=a';
    ub=b';
    obj_no=numel(feval(fun,a));
    

Best_universe=zeros(1,dim);
[ff,gg]=feval(fun,Best_universe');
Best_universe_Inflation_rate=inf*ones(1,obj_no); %change this to -inf for maximization problems
Best_universe_G=inf(1,numel(gg));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Archive_X=zeros(ArchiveMaxSize,dim);
Archive_F=ones(ArchiveMaxSize,obj_no)*inf;
Archive_G=ones(ArchiveMaxSize,numel(gg))*inf;
Archive_member_no=0;
%Minimum and maximum of Wormhole Existence Probability (min and max in
% Eq.(3.3) in the paper
WEP_Max=1;
WEP_Min=0.2;
%Initialize the positions of universes
Universes=initialization(N,dim,ub,lb);
%Iteration(time) counter
Time=1;

%Main loop
while Time<Max_time+1
    
    %Eq. (3.3) in the original MVO paper
    WEP=WEP_Min+Time*((WEP_Max-WEP_Min)/Max_time);
    
    %Travelling Distance Rate (Formula): Eq. (3.4) in the original MVO paper
    TDR=1-((Time)^(1/6)/(Max_time)^(1/6));
    
    %Inflation rates (I) (fitness values)
    %Inflation_rates=zeros(1,size(Universes,1));
    
    for i=1:size(Universes,1)
        
        %Boundary checking (to bring back the universes inside search
        % space if they go beyoud the boundaries
        Flag4ub=Universes(i,:)>ub;
        Flag4lb=Universes(i,:)<lb;
        Universes(i,:)=(Universes(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        %Calculate the inflation rate (fitness) of universes
%         Inflation_rates(i,:)=fobj(Universes(i,:));
        
        % Objective Function
        [fi,gi]=feval(fobj,Universes(i,:)');
        % Penalty Function
        if max(gi)>0
            Inflation_rates(i,:)=fpenal_1(fi',gi');
        else
            Inflation_rates(i,:)=fi';
        end
        G(i,:)=gi';
        %
        %
        
        %Elitism
        if dominates(Inflation_rates(i,:),G(i,:),Best_universe_Inflation_rate,Best_universe_G)
            Best_universe_Inflation_rate=Inflation_rates(i,:);
            Best_universe=Universes(i,:);
            Best_universe_G=G(i,:);
        end
        
    end
    
    [sorted_Inflation_rates,sorted_indexes]=sort(Inflation_rates);
    
    for newindex=1:N
        Sorted_universes(newindex,:)=Universes(sorted_indexes(newindex),:);
    end
    
    %Normaized inflation rates (NI in Eq. (3.1) in the original MVO paper)
    normalized_sorted_Inflation_rates=normr(sorted_Inflation_rates);
    
    Universes(1,:)= Sorted_universes(1,:);
    
    [Archive_X, Archive_F, Archive_G, Archive_member_no]=UpdateArchive(Archive_X, Archive_F, Archive_G, Universes, Inflation_rates, G, Archive_member_no);
    if Archive_member_no>ArchiveMaxSize
        Archive_mem_ranks=RankingProcess(Archive_F, Archive_G, ArchiveMaxSize, obj_no);
        [Archive_X, Archive_F, Archive_G, Archive_mem_ranks, Archive_member_no]=HandleFullArchive(Archive_X, Archive_F, Archive_G, Archive_member_no, Archive_mem_ranks, ArchiveMaxSize);
    else
        Archive_mem_ranks=RankingProcess(Archive_F, Archive_G, ArchiveMaxSize, obj_no);
    end
    Archive_mem_ranks=RankingProcess(Archive_F, Archive_G, ArchiveMaxSize, obj_no);
    % to improve coverage
    index=RouletteWheelSelection(1./Archive_mem_ranks);
    if index==-1
        index=1;
    end
   Best_universe_Inflation_rate=Archive_F(index,:);
   Best_universe=Archive_X(index,:); 
   
  
    
    %Update the Position of universes
    for i=2:size(Universes,1)%Starting from 2 since the firt one is the elite
        Back_hole_index=i;
        for j=1:size(Universes,2)
            r1=rand();
            if r1<normalized_sorted_Inflation_rates(i)
                White_hole_index=RouletteWheelSelection(-sorted_Inflation_rates);% for maximization problem -sorted_Inflation_rates should be written as sorted_Inflation_rates
                if White_hole_index==-1
                    White_hole_index=1;
                end
                %Eq. (3.1) in the paper
                Universes(Back_hole_index,j)=Sorted_universes(White_hole_index,j);
            end
            
            if (size(lb',1)==1)
                %Eq. (3.2) in the original MVO paper if the boundaries are all the same
                r2=rand();
                if r2<WEP
                    r3=rand();
                    if r3<0.5
                        Universes(i,j)=Best_universe(1,j)+TDR*((ub-lb)*rand+lb);
                    end
                    if r3>0.5
                        Universes(i,j)=Best_universe(1,j)-TDR*((ub-lb)*rand+lb);
                    end
                end
            end
            
            if (size(lb',1)~=1)
                %Eq. (3.2) in the original MVO paper if the upper and lower bounds are
                %different for each variables
                r2=rand();
                if r2<WEP
                    r3=rand();
                    if r3<0.5
                        Universes(i,j)=Best_universe(1,j)+TDR*((ub(j)-lb(j))*rand+lb(j));
                    end
                    if r3>0.5
                        Universes(i,j)=Best_universe(1,j)-TDR*((ub(j)-lb(j))*rand+lb(j));
                    end
                end
            end
            
        end
    end   
%  display(['At the iteration ', num2str(Time), ' there are ', num2str(Archive_member_no), ' non-dominated solutions in the archive']);

    % Save results
    [ppareto,fpareto,gpareto,~]=pbil_selection0([],[],[],Archive_X',Archive_F',Archive_G',narchive);
    
    rst.ppareto{Time}=ppareto;
    rst.fpareto{Time}=fpareto;
    rst.gpareto{Time}=gpareto;
    rst.timestamp=datetime('now');
        
    Time=Time+1;
    
end

%     ppareto=Archive_X';
%     for i=1:size(Archive_X,1)
%         [fpareto(:,i),gpareto(:,i)]=feval(fun,ppareto(:,i));
%     end
%     
% %     [ppareto,fpareto,gpareto,~]=pbil_selection0([],[],[],ppareto,fpareto,gpareto,narchive);
%     save(fout,'ppareto','fpareto','gpareto')
end

%%%%%%%%%%%%%%%%
% Sub-Function %
%%%%%%%%%%%%%%%%
function DELTA= Diversity_metric_delta(pareto_fun,Factual)

%-----------------2D--------------
[~,index]=sort(pareto_fun(:,1));

B=pareto_fun(index,:);
row=size(B,1);

for i=1:row-1
      Distance(i)=norm(B(i,:)-B(i+1,:));    
end

d_average= (1/row)*sum(Distance);
%----------------------3D------------
% row=size(pareto_fun,1);
% for i=1:row
%     D1=repmat(pareto_fun(i,:),row,1);
%     D2=D1-pareto_fun;
%     for j=1:row
%         D3(j)=norm(D2(j,:));        
%     end
%     D3(i)=[];
%     Distance(i)=min(D3);    
% end

% d_average= (1/row)*sum(Distance);
%--------------------------------------
total_distance=0;
for i=1:row-1
 total_distance=total_distance+abs((Distance(i)-d_average)); 
end
%------------------------------------2D-------------------
[~,index]=sort(Factual(:,1));
BB=Factual(index,:);

dF=norm(B(1,:)-BB(1,:));
dL=norm(B(end,:)-BB(end,:));

DELTA=(dF+dL+total_distance)/(dF+dL+((row-1)*d_average));
%----------------------------3D--------------------
% for j=1:size(Fmin,2)
%    [~,index1]=max(Factual(:,j));FF=Factual(index1,:);
%    [~,index1]=max(Fmin(:,j));ff=Fmin(index1,:);
%     DFL(j)=norm(FF-ff);
% end
% 
% DELTA=(sum(DFL)+total_distance)/(sum(DFL)+((row)*d_average));
end

function o=dominates(f1,g1,f2,g2)
    if numel(f1)>0 && numel(g1)>0 && numel(f2)>0 && numel(g2)>0
        if max(g1)<=0 && max(g2)<=0
            o=all(f1<=f2) && any(f1<f2);
        elseif max(g1)<=0 && max(g2)>0
            o=1;
        elseif max(g1)>0 && max(g2)<=0
            o=0;
        elseif max(g1)>0 && max(g2)>0
            if max(g1)>=max(g2)
                o=0;
            else
                o=1;
            end
        else
            error('unknown error, please check');
        end
    elseif numel(f1)>0 && numel(g1)>0 && numel(f2)==0 && numel(g2)==0
        o=1;
    elseif numel(f1)==0 && numel(g1)==0 && numel(f2)>0 && numel(g2)>0
    else
        o=0;
    end
end

function eps = epsilon_matlab(PF, truePF)

%STEP 1. Obtain the maximum and minimum values of the Pareto front
m1 = size(PF, 1);
m = size(truePF, 1);

%STEP 2. find the epsilon value
for i = 1:m
    temp = PF - repmat(truePF(i,:), m1, 1);
    eps_k = max(temp, [], 2);
    eps_j = min(eps_k);
    
    if i == 1
        eps = eps_j;
    elseif eps < eps_j
        eps = eps_j;
    end
end

end

function GD = GD_matlab(PF, truePF)

q = 2; %the parameter of GD
%STEP 1. Obtain the maximum and minimum values of the Pareto front
m1 = size(PF, 1);
m = size(truePF, 1);
maxVals = max(truePF);
minVals = min(truePF);

%STEP 2. Get the normalized front
normalizedPF = (PF - repmat(minVals, m1, 1)) ./ repmat(maxVals - minVals, m1, 1);
normalizedTruePF = (truePF - repmat(minVals, m, 1)) ./ repmat(maxVals - minVals, m, 1);

%STEP 3. Sum the distances between each point of the front and the nearest point in the true Pareto front
GD = 0;
for i = 1:m1
    diff = repmat(normalizedPF(i,:), m, 1) - normalizedTruePF;
    dist = sqrt(sum(diff.^2, 2));         
    GD = GD + min(dist)^q;
end
GD = GD^(1.0/q)/m1;

end

function generalizedspread = GeneralizedSpread_matlab(PF, truePF)

%STEP 1. Obtain the maximum and minimum values of the Pareto front
m1 = size(PF, 1);
m = size(truePF, 1);
maxVals = max(truePF);
minVals = min(truePF);

%STEP 2. Get the normalized front
normalizedPF = (PF - repmat(minVals, m1, 1)) ./ repmat(maxVals - minVals, m1, 1);
normalizedTruePF = (truePF - repmat(minVals, m, 1)) ./ repmat(maxVals - minVals, m, 1);

%STEP 3. find extrmal values
[~, I] = max(normalizedTruePF);
extremValues = normalizedTruePF(I, :);

%STEP 4. Sort normalizedFront and normalizedParetoFront;
normalizedPF = sortrows(normalizedPF);

%STEP 5. Calculate the metric value. The value is 1.0 by default
diff = normalizedPF(1, :) - normalizedPF(m1-1, :);
dist = sqrt(sum(diff.^2, 2));
if dist == 0.0
    generalizedspread = 1.0;
else
    %STEP 6. Calculate the mean distance between each point and its nearest neighbor
    distance_near = zeros(m1, 1);
    for i = 1:m1
        diff = repmat(normalizedPF(i,:), m1, 1) - normalizedPF;
        dist = sqrt(sum(diff.^2, 2));
        dist(i) = realmax; %maximum value in double
        distance_near(i) = min(dist);
    end
    dmean = mean(distance_near);
    % STEP 7. Calculate the distance to extremal values
    dExtrems = 0;
    for i = 1:size(extremValues, 1)
        diff = repmat(extremValues(i,:), m1, 1) - normalizedPF;
        dist = sqrt(sum(diff.^2, 2));
        dExtrems = dExtrems + min(dist);
    end
    %STEP 8. Computing the value of the metric
    sumVal = sum(abs(distance_near - repmat(dmean, m1, 1)));
    generalizedspread = (dExtrems + sumVal) / (dExtrems + m1 * dmean);
end

end

function GD=Generational_distance(pareto_fun,Factual)

[row,~]=size(Factual);
npf=size(pareto_fun,1);

total_distance=0;
for i=1:npf
    D1=repmat(pareto_fun(i,:),row,1);
    D2=D1-Factual;
    for j=1:row
        D3(j)=norm(D2(j,:));      
    end
    Distance=min(D3);
    total_distance=total_distance+(Distance^2);
    % total_distance=total_distance+(Distance);
end


%GD=sqrt((1/npf)*total_distance);
 GD=(1/npf)*sqrt(total_distance);
% GD=((1/npf)*total_distance);


end

function [Archive_X_Chopped, Archive_F_Chopped, Archive_G_Chopped, Archive_mem_ranks_updated, Archive_member_no]=HandleFullArchive(Archive_X, Archive_F, Archive_G, Archive_member_no, Archive_mem_ranks, ArchiveMaxSize)
    for i=1:size(Archive_F,1)-ArchiveMaxSize
        index=RouletteWheelSelection(Archive_mem_ranks);
        %[value index]=min(Archive_mem_ranks);
        
        if index==-1
        	index=ceil(size(Archive_F,1)*rand(size(Archive_F,1)-ArchiveMaxSize,1));
        end
    
        Archive_X=[Archive_X(1:index-1,:) ; Archive_X(index+1:Archive_member_no,:)];
        Archive_F=[Archive_F(1:index-1,:) ; Archive_F(index+1:Archive_member_no,:)];
        Archive_G=[Archive_G(1:index-1,:) ; Archive_G(index+1:Archive_member_no,:)];
        Archive_mem_ranks=[Archive_mem_ranks(1:index-1) Archive_mem_ranks(index+1:Archive_member_no)];
        Archive_member_no=Archive_member_no-1;
    end

    Archive_X_Chopped=Archive_X;
    Archive_F_Chopped=Archive_F;
    Archive_G_Chopped=Archive_G;
    Archive_mem_ranks_updated=Archive_mem_ranks;
end

function [hv, nrW] = hypeIndicatorExact8(A, bounds, dim)

hv  = 0; nrW = 0;
nrA = size(A, 1);
if nrA == 0
    return;
elseif nrA == 1
    hv = prod(bounds - A(1,:));
    nrW = 0;
elseif nrA == 2
    w = max(A(1,:), A(2,:));       
    hv = prod(bounds - A(1,:)) + prod(bounds - A(2,:)) - prod(bounds - w);
%     nrW = 1;   
elseif nrA == 3
    w1 = max(A(1,:), A(2,:));       
    w2 = max(A(1,:), A(3,:));
    w3 = max(A(2,:), A(3,:));
    w4 = max(w1, A(3,:));
    hv = prod(bounds - A(1,:)) + prod(bounds - A(2,:)) + prod(bounds - A(3,:)) - prod(bounds - w1) - prod(bounds - w2) - prod(bounds - w3) + prod(bounds - w4);
elseif dim == 2
        %the points are slice to 2d, same values in other dimensions
        hv = hypervolume2(A(:,1:2), bounds(1:2));
        for i = 3:size(A,2)
            hv = hv * (bounds(i) - A(1,i));
        end
        nrW = 0;
else    
    A = sortrows(A, dim);
    while A(nrA, dim) == A(1, dim)
%         fprintf('s');
        dim = dim - 1;
        A = sortrows(A, dim);
    end
    if dim == 2
        %the points are slice to 2d, same values in other dimensions
        hv = hypervolume2(A(:,1:2), bounds(1:2));
        for i = 3:size(A,2)
            hv = hv * (bounds(i) - A(1,i));
        end
        nrW = 0;
        %     nrW = nrA - 1;
    else
        hv = hv + prod(bounds - A(1,:));
        for i = nrA:-1:2
            v = prod(bounds - A(i,:));
            if v > 0
                W = worse(A(i,:), A(1:i-1,:));
                [hv1, nrW1] = hypeIndicatorExact8(W, bounds, dim - 1);
                hv = hv + v - hv1;
                nrW = nrW + size(W,1) +nrW1;
            end
        end
    end
end

end


function [W] = worse(A1, A2) 
%how to quickly get the worse set

[nrA1, dim] = size(A1);
nrA2 = size(A2,1);
W = zeros(nrA1 * nrA2, dim);
nrW = 0;
for i = 1:nrA1
    worse = max(repmat(A1(i,:), nrA2, 1), A2);    
    %find the non-dominated points in worse    
    for j = nrA2:-1:1  
        %the index (nrA2:-1:1) is important, because the last one is alwasys better (nondominated) than others
        [W, nrW] = insertPoint(worse(j,:), W, nrW);
    end
end
%find the non-equal, and non-dominated points in W
W = W(1:nrW,:);

end

function [W, nrW] = insertPoint(p, W, nrW)

if nrW == 0
    W(1,:) = p;
    nrW = 1;
else
    dim = size(p,2);
    i = 1;
    while i <= nrW
        diff = p - W(i,:);
        if sum(diff >=0, 2) == dim
            %same point, or p is dominated
%             fprintf('%c', '1');  
            return;
        end
        if sum(diff <=0, 2) == dim
            %some points in W are dominaed by p
            W(i, :) = W(nrW, :);
%             fprintf('%c', '2'); 
            i = i - 1;
            nrW = nrW - 1;
        end
        i = i + 1;
    end
    %non-dominated each other, add this point
    nrW = nrW + 1;
    W(nrW,:) = p;    
end
    
end

function [hv] = hypervolume2(A, bounds)
%%get hypervolume for 2D

hv = 0;
[nrA, dim] = size(A);
if nrA == 0
    return;
end
A = sortrows(A, dim);
for i=1:nrA-1
    hv = hv + (A(i+1, 2) - A(i, 2)) * (bounds(1) - A(i,1));
end
hv = hv + (bounds(2) - A(nrA, 2)) * (bounds(1) - A(nrA,1));

end

function IGD = IGD_matlab(PF, truePF)

q = 2; %the parameter of IGD
%STEP 1. Obtain the maximum and minimum values of the Pareto front
m1 = size(PF, 1);
m = size(truePF, 1);
maxVals = max(truePF);
minVals = min(truePF);

%STEP 2. Get the normalized front
normalizedPF = (PF - repmat(minVals, m1, 1)) ./ repmat(maxVals - minVals, m1, 1);
normalizedTruePF = (truePF - repmat(minVals, m, 1)) ./ repmat(maxVals - minVals, m, 1);

%STEP 3. Sum the distances between each point of the front and the nearest point in the true Pareto front
IGD = 0;
for i = 1:m
    diff = repmat(normalizedTruePF(i,:), m1, 1) - normalizedPF;
    dist = sqrt(sum(diff.^2, 2));         
    IGD = IGD + min(dist)^q;
end
IGD = IGD^(1.0/q)/m;

end

% This function initialize the first population of search agents
function X=initialization(SearchAgents_no,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    ub_new=ones(1,dim)*ub;
    lb_new=ones(1,dim)*lb;
else
    ub_new=ub;
    lb_new=lb;
end

% If each variable has a different lb and ub
for i=1:dim
    ub_i=ub_new(i);
    lb_i=lb_new(i);
    X(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
end

X=X;
end

function MS=metric_of_maximum_spread(pareto_fun,Factual)

[~,col]=size(Factual);

ms=0;
for i=1:col
   ms=ms+((max(pareto_fun(:,i))-min(pareto_fun(:,i)))/(max(Factual(:,i))-min(Factual(:,i))))^2;
   % ms=ms+((max(Fmin(:,i))-min(Fmin(:,i)))^2);
end
 
MS=sqrt(ms/col);
% MS=sqrt(ms);  

end

function S=metric_of_spacing(pareto_fun)

npf=size(pareto_fun,1);

for i=1:npf
    D1=repmat(pareto_fun(i,:),npf,1);
    D2=D1-pareto_fun;
    for j=1:npf
        % D3(j)=norm(D2(j,:)); 
        D3(j)=sum(abs(D2(j,:)));   
    end
    D3(i)=[];
    Distance(i)=min(D3);    
end

d_average=(1/npf)*sum(Distance);

total_distance=0;
for i=1:npf
 total_distance=total_distance+((Distance(i)-d_average)^2); 
end

%S=(sqrt(total_distance/npf))/d_average;
 S=(sqrt(total_distance/(npf-1)));
% S=(sqrt(total_distance/(npf)));

end

function ranks=RankingProcess(Archive_F, Archive_G, ArchiveMaxSize, obj_no)

my_min=min(Archive_F);
my_max=max(Archive_F);

if size(Archive_F,1)==1
    my_min=Archive_F;
    my_max=Archive_F;
end

r=(my_max-my_min)/(20);

ranks=zeros(1,size(Archive_F,1));

for i=1:size(Archive_F,1)
    ranks(i)=0;
    for j=1:size(Archive_F,1)
        flag=0; % a flag to see if the point is in the neoghbourhood in all dimensions.
        for k=1:obj_no
            %
            %
            F1=fpenal_1(Archive_F(j,k),Archive_G(j,:));
            F2=fpenal_1(Archive_F(i,k),Archive_G(i,:));
            if (abs(F1-F2)<r(k))
                flag=flag+1;
            end
            %
            %
        end
        if flag==obj_no
            ranks(i)=ranks(i)+1;
        end
    end
end
end

function o = RouletteWheelSelection(weights)
accumulation = cumsum(weights);
p = rand() * accumulation(end);
chosen_index = -1;
for index = 1 : length(accumulation)
    if (accumulation(index) > p)
        chosen_index = index;
        break;
    end
end
o = chosen_index;
end

function spread = Spread_matlab(PF, truePF)

%STEP 1. Obtain the maximum and minimum values of the Pareto front
m1 = size(PF, 1);
m = size(truePF, 1);
maxVals = max(truePF);
minVals = min(truePF);

%STEP 2. Get the normalized front
normalizedPF = (PF - repmat(minVals, m1, 1)) ./ repmat(maxVals - minVals, m1, 1);
normalizedTruePF = (truePF - repmat(minVals, m, 1)) ./ repmat(maxVals - minVals, m, 1);

%STEP 3. Sort normalizedFront and normalizedParetoFront;
normalizedPF = sortrows(normalizedPF);
normalizedTruePF = sortrows(normalizedTruePF);

%STEP 4. Compute df and dl (See specifications in Deb's description of the metric)
diff = normalizedPF(1, :) - normalizedTruePF(1, :);
df = sqrt(sum(diff.^2, 2));
diff = normalizedPF(m1, :) - normalizedTruePF(m, :);
dl = sqrt(sum(diff.^2, 2));

%STEP 5. Calculate the mean of distances between points i and (i - 1). (the poins are in lexicografical order)
diff = normalizedPF(1:m1-1,:) - normalizedPF(2:m1,:);
distNears = sqrt(sum(diff.^2, 2));
meanVal = mean(distNears);

%STEP 6. If there are more than a single point, continue computing the metric. In other case, return the worse value (1.0, see metric's description).
if m1 > 1
    diversitySum = sum(abs(distNears - repmat(meanVal, m1-1, 1)));
    spread = (df +dl + diversitySum) / (df + dl + (m1-1) * meanVal);
else
    spread = 1.0;
end

end

function [Archive_X_updated, Archive_F_updated, Archive_G_updated, Archive_member_no]=UpdateArchive(Archive_X, Archive_F, Archive_G, Particles_X, Particles_F, Particles_G, Archive_member_no)
    [Archive_X_temp,uind]=unique([Archive_X ; Particles_X],'rows');
    Archive_F_temp=[Archive_F ; Particles_F];
    Archive_F_temp=Archive_F_temp(uind,:);
    Archive_G_temp=[Archive_G ; Particles_G];
    Archive_G_temp=Archive_G_temp(uind,:);
    
    o=zeros(1,size(Archive_F_temp,1));

    for i=1:size(Archive_F_temp,1)
        o(i)=0;
        for j=1:i-1
            if any(Archive_F_temp(i,:) ~= Archive_F_temp(j,:))
                if dominates(Archive_F_temp(i,:),Archive_G_temp(i,:),Archive_F_temp(j,:),Archive_G_temp(j,:))
                    o(j)=1;
                elseif dominates(Archive_F_temp(j,:),Archive_G_temp(j,:),Archive_F_temp(i,:),Archive_G_temp(i,:))
                    o(i)=1;
                    break;
                end
            else
                o(j)=0;
                o(i)=0;
            end
        end
    end

    Archive_member_no=0;
    index=0;
    for i=1:size(Archive_X_temp,1)
        if o(i)==0
            Archive_member_no=Archive_member_no+1;
            Archive_X_updated(Archive_member_no,:)=Archive_X_temp(i,:);
            Archive_F_updated(Archive_member_no,:)=Archive_F_temp(i,:);
            Archive_G_updated(Archive_member_no,:)=Archive_G_temp(i,:);
        else
            index=index+1;
            %         dominated_X(index,:)=Archive_X_temp(i,:);
            %         dominated_F(index,:)=Archive_F_temp(i,:);
        end
    end
end


function [pareto1,fpareto1,gpareto1,A]=pbil_selection0(x1,f1,g1,pareto,fpareto,gpareto,narchive)
    %
    % GA Elite strategy
    % keep 1 elite from the old generation and another from the
    % new generation
    x=[pareto x1];
    f=[fpareto f1];
    g=[gpareto g1];

    [m0,n0]=size(fpareto);
    [m1,n1]=size(x);
    
    for i=1:n1
        xi=x(:,i);
        fi=f(:,i);
        gi=g(:,i);
        A(i,i)=0;
        for j=(i+1):n1
            xj=x(:,j);
            fj=f(:,j);
            gj=g(:,j);
            %%%%%%%%%%%%%%%%%%%%%%%%%
            [p_count1,p_count2]=fdominated(fi,gi,fj,gj);
            A(i,j)=p_count1;
            A(j,i)=p_count2;
            %%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    B=sum(A,1);
    Indm=[];
    for i=1:n1
        if B(i)==0
            Indm=[Indm i];
        end
    end
    nndm=length(Indm);

    pareto1=x(:,Indm);
    fpareto1=f(:,Indm);
    gpareto1=g(:,Indm);
    A=A(Indm,Indm);
end

function [p1,p2]=fdominated(f1,g1,f2,g2)
    n=length(f1);
    mg1=max(g1);
    mg2=max(g2);

    icount11=0;
    icount12=0;
    icount21=0;
    icount22=0;

    if mg1<=0 && mg2<=0
        for i=1:n
            if f1(i) <= f2(i)
                icount11=icount11+1;
            end
            if f1(i) < f2(i)
                icount12=icount12+1;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%
            if f2(i) <= f1(i)
                icount21=icount21+1;
            end
            if f2(i) < f1(i)
                icount22=icount22+1;
            end
        end
        if icount11 == n && icount12 > 0
            p1=1;
        else
            p1=0;
        end
        if icount21 == n && icount22 > 0
            p2=1;
        else
            p2=0;
        end
    elseif mg1 <=0 && mg2 > 0
        p1=1;p2=0;
    elseif mg2 <=0 && mg1 > 0
        p1=0;p2=1;
    else
        if mg1 <= mg2
            p1=1;p2=0;
        else
            p1=0;p2=1;
        end
    end
end

function [fp] = fpenal_1(f,g)
    % penalty function
    if max(g)>0
        fp=f+100*max(g);
    else
        fp=f;
    end
end