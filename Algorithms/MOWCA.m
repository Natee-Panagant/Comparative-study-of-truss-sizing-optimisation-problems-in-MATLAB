function rst = MOWCA(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)
% function [Non_Dominated_Solutions,Pareto_Front,NFEs,Elapsed_Time]=MOWCA_Unconstrained(objective_function,LB,UB,nvars,Npop,Nsr,dmax,max_it)
% %% Information
% Last Revised: 10th April 2016

% Multi-objective Water Cycle Algorithm (MOWCA) (Standard Version)

% This code is prepared for multiple objective functions, unconstrained, and continuous problems.
% Note that in order to obey the copy-right rules, please kindly cite our published paper properly. Thank you.

% A. Sadollah, H. Eskandar, J.H. Kim, A. Bahreininejad, “Water cycle algorithm for solving multi-objective optimization problems”, Soft Computing, 19(9) (2015) 2587-2603.
% A. Sadollah, H. Eskandar, J.H. Kim, “Water cycle algorithm for solving constrained multi-objective problems”, Applied Soft Computing, 27 (2015) 279-298.

% INPUTS:

% objective_function:           Objective functions which you wish to minimize or maximize
% LB:                           Lower bound of a problem
% UB:                           Upper bound of a problem
% nvars:                        Number of design variables
% Npop                          Population size
% Nsr                           Number of rivers + sea
% dmax                          Evporation condition constant
% max_it:                       Maximum number of iterations

% OUTPUTS:

% Non_Dominated_Solutions:      Optimum non-dominated solutions
% Pareto_Front:                 Obtained optimum objective function values
% NFEs:                         Number of function evaluations
% Elapsed_Time                  Elasped time for solving an optimization problem
% %% Default Values for MOWCA
% format long g
% if (nargin <5 || isempty(Npop)), Npop=50; end
% if (nargin <6 || isempty(Nsr)), Nsr=4; end
% if (nargin <7 || isempty(dmax)), dmax=1e-16; end
% if (nargin <8 || isempty(max_it)), max_it=100; end
% %% --------------------------------------------------------------------------
% % Create initial population
% % tic

%
%
rand('state',sum(100*clock));
tic
%
Npop=nsol;
Nsr=4;
dmax=1e-16;
max_it=nloop;
nvars=nvar;
LB=a';
UB=b';
%
%
N_stream=Npop-Nsr;
NPF=Npop;    % Pareto Front Archive Size

ind.Position=[];
ind.Cost=[];
ind.Rank=[];
ind.DominationSet=[];
ind.DominatedCount=[];
ind.CrowdingDistance=[];

pop=repmat(ind,Npop,1);

for i=1:Npop
    pop(i).Position=LB+(UB-LB).*rand(1,nvar);
%     pop(i).Cost=objective_function(pop(i).Position);
    
    % Objective Function
    [fi,gi]=feval(fun,pop(i).Position');
    % Penalty Function
    if max(gi)>0
        pop(i).Cost=fpenal_1(fi',gi');
    else
        pop(i).Cost=fi';
    end
    pop(i).G=gi';
    %
    %
end

[pop, F]=NonDominatedSorting(pop);  % Non-dominated sorting
pop=CalcCrowdingDistance(pop,F);   % Calculate crowding distance
pop=SortPopulation(pop);     % Sort population
%------------- Forming Sea, Rivers, and Streams  --------------------------
sea=pop(1);
river=pop(2:Nsr);
stream=pop(Nsr+1:end);

cs=[sea.CrowdingDistance';[river.CrowdingDistance]';stream(1).CrowdingDistance];

f=0;
if length(unique(cs))~=1
    CN=cs-max(cs);
else
    CN=cs;
    f=1;
end

NS=round(abs(CN/(sum(CN)+eps))*N_stream);

if f~=1
    NS(end)=[];
end
NS=sort(NS,'descend');
% ------------------------- Modification on NS -----------------------
i=Nsr;
while sum(NS)>N_stream
    if NS(i)>1
        NS(i)=NS(i)-1;
    else
        i=i-1;
    end
end

i=1;
while sum(NS)<N_stream
    NS(i)=NS(i)+1;
end

if find(NS==0)
    index=find(NS==0);
    for i=1:size(index,1)
        count=0;
        while NS(index(i))==0
            NS(index(i))=NS(index(i))+round(NS(i)/6);
            NS(i)=NS(i)-round(NS(i)/6);
            
            count=count+1;
            if count>100
                break;
            end
        end
    end
end

NS=sort(NS,'descend');
NB=NS(2:end);
% %%
%----------- Main Loop for MOWCA --------------------------------------------
% disp('********** Multi-objective Water Cycle Algorithm (MOWCA)************');
% disp('*Iterations         Number of Pareto Front Members *');
% disp('********************************************************************');
FF=zeros(max_it,numel(sea.Cost));
for i=1:max_it
    %---------- Moving stream to sea---------------------------------------
    for j=1:NS(1)
        stream(j).Position=stream(j).Position+2.*rand(1).*(sea.Position-stream(j).Position);
        
        stream(j).Position=min(stream(j).Position,UB);
        stream(j).Position=max(stream(j).Position,LB);
        
%         stream(j).Cost=objective_function(stream(j).Position);
        
        % Objective Function
        [fi,gi]=feval(fun,stream(j).Position');
        % Penalty Function
        if max(gi)>0
            stream(j).Cost=fpenal_1(fi',gi');
        else
            stream(j).Cost=fi';
        end
        stream(j).G=gi';
        %
        %
    end
    %---------- Moving Streams to rivers-----------------------------------
    for k=1:Nsr-1
        for j=1:NB(k)
            stream(j+sum(NS(1:k))).Position=stream(j+sum(NS(1:k))).Position+2.*rand(1,nvars).*(river(k).Position-stream(j+sum(NS(1:k))).Position);
            
            stream(j+sum(NS(1:k))).Position=min(stream(j+sum(NS(1:k))).Position,UB);
            stream(j+sum(NS(1:k))).Position=max(stream(j+sum(NS(1:k))).Position,LB);
            
%             stream(j+sum(NS(1:k))).Cost=objective_function(stream(j+sum(NS(1:k))).Position);
            
            % Objective Function
            [fi,gi]=feval(fun,stream(j+sum(NS(1:k))).Position');
            % Penalty Function
            if max(gi)>0
                stream(j+sum(NS(1:k))).Cost=fpenal_1(fi',gi');
            else
                stream(j+sum(NS(1:k))).Cost=fi';
            end
            stream(j+sum(NS(1:k))).G=gi';
            %
            %
        end
    end
    %---------- Moving rivers to Sea --------------------------------------
    for j=1:Nsr-1
        river(j).Position=river(j).Position+2.*rand(1,nvars).*(sea.Position-river(j).Position);
        
        river(j).Position=min(river(j).Position,UB);
        river(j).Position=max(river(j).Position,LB);
        
%         river(j).Cost=objective_function(river(j).Position);
            
        % Objective Function
        [fi,gi]=feval(fun,river(j).Position');
        % Penalty Function
        if max(gi)>0
            river(j).Cost=fpenal_1(fi',gi');
        else
            river(j).Cost=fi';
        end
         river(j).G=gi';
        %
        %
    end
    %-------------- Evaporation condition and raining process--------------
    % Check the evaporation condition for rivers and sea
    for k=1:Nsr-1
        if ((norm(river(k).Position-sea.Position)<dmax) || rand<0.1)
            for j=1:NB(k)
                stream(j+sum(NS(1:k))).Position=LB+rand(1,nvars).*(UB-LB);
            end
        end
    end
    % Check the evaporation condition for streams and sea
    for j=1:NS(1)
        if ((norm(stream(j).Position-sea.Position)<dmax))
            stream(j).Position=LB+rand(1,nvars).*(UB-LB);
        end
    end
    %----------------------------------------------------------------------
    pop=[pop;sea;river;stream];
    
    [pop, F]=NonDominatedSorting(pop);  % Non-dominated sorting
    pop=CalcCrowdingDistance(pop,F);   % Calculate crowding distance
    pop=SortPopulation(pop);     % Sort population
    
    pop=pop(1:NPF);
    
    [pop, F]=NonDominatedSorting(pop);  % Non-dominated sorting
    pop=CalcCrowdingDistance(pop,F);   % Calculate crowding distance
    [pop, F]=SortPopulation(pop);     % Sort population
    
    sea=pop(1);river=pop(2:Nsr);stream=pop(Nsr+1:end);
    F1=pop(F{1});
    
    dmax=dmax-(dmax/max_it);
%     disp(['Iteration: ',num2str(i),'      ((Number of Members in the Pareto Front)):  ',num2str(numel(F{1}))]);
    FF(i,:)=sea.Cost;
    
    % Plot F1 costs
%     figure(1);
%     PlotCosts(F1);
%     pause(0.01)

    %save results
    ppareto=reshape([pop(F{1}).Position],nvar,[]);
    fpareto=reshape([pop(F{1}).Cost],numel(sea(1).Cost),[]);
    gpareto=reshape([pop(F{1}).G],numel(sea(1).G),[]);
    [ppareto,fpareto,gpareto,~]=pbil_selection0([],[],[],ppareto,fpareto,gpareto,narchive);
    
    rst.ppareto{i}=ppareto;
    rst.fpareto{i}=fpareto;
    rst.gpareto{i}=gpareto;
    rst.timestamp=datetime('now');
end
% %% Optimization Results
% toc;
% Elapsed_Time=toc;
% NFEs=Npop*max_it;

% k=1;
% a=[river.Position stream.Position];
% A=zeros(N_stream+Nsr-1,nvars);
% for i=1:N_stream+Nsr-1
%     A(i,:)=a(1,k:i*nvars);
%     k=k+nvars;
% end
% Non_Dominated_Solutions=[sea.Position;A];
% 
% b=[[river.Cost]';[stream.Cost]'];
% Pareto_Front=[sea.Cost';b];

%     ppareto=Non_Dominated_Solutions';
%     for i=1:size(Non_Dominated_Solutions,1)
%         [fpareto(:,i),gpareto(:,i)]=feval(fun,ppareto(:,i));
%     end
    
%     [ppareto,fpareto,gpareto,~]=pbil_selection0([],[],[],ppareto,fpareto,gpareto,narchive);
%     save(fout,'ppareto','fpareto','gpareto')
end

%%%%%%%%%%%%%%%%
% Sub-Function %
%%%%%%%%%%%%%%%%
function pop=CalcCrowdingDistance(pop,F)

nF=numel(F);

for k=1:nF
    
    Costs=[pop(F{k}).Cost];
    nObj=size(Costs,1);
    n=numel(F{k});
    d=zeros(n,nObj);
    
    for j=1:nObj
        [cj, so]=sort(Costs(j,:));
        d(so(1),j)=1e3;
        
        for i=2:n-1
            d(so(i),j)=abs(cj(i+1)-cj(i-1))/abs(cj(1)-cj(end));
        end
        d(so(end),j)=1e3;
    end
    
    for i=1:n
        pop(F{k}(i)).CrowdingDistance=sum(d(i,:));
    end
end
end

function DELTA=Diversity_metric_delta(Pareto_Front,Factual)

[~,index]=sort(Pareto_Front(:,1));

B=Pareto_Front(index,:);
row=size(B,1);

for i=1:row-1
    Distance(i)=norm(B(i,:)-B(i+1,:));
end

d_average= (1/row)*sum(Distance);

total_distance=0;
for i=1:row-1
    total_distance=total_distance+abs((Distance(i)-d_average));
end

[~,index]=sort(Factual(:,1));
BB=Factual(index,:);

dF=norm(B(1,:)-BB(1,:));
dL=norm(B(end,:)-BB(end,:));

DELTA=(dF+dL+total_distance)/(dF+dL+((row-1)*d_average));
end

function o = Dominates(f1,g1,f2,g2)
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

function Feval=FON(x)

global nvars

a=0;
b=0;
for i=1:nvars
    a=a+(x(i)-(1/8^0.5))^2;
    b=b+(x(i)+(1/8^0.5))^2;
end
f1=1-exp(-a);
f2=1-exp(-b);

Feval=[f1;f2];

end

function GD=Generational_distance(Pareto_Front,Factual)

row=size(Factual,1);

total_distance=0;
for i=1:size(Pareto_Front,1)
    D1=repmat(Pareto_Front(i,:),row,1);
    D2=D1-Factual;
    for j=1:row
        D3(j)=norm(D2(j,:));        
    end
    Distance=min(D3);
    total_distance = total_distance + (Distance^2);
end

GD=sqrt(((1/size(Pareto_Front,1))*total_distance));
end

function MS=metric_of_maximum_spread(Pareto_Front,Factual)

[~,cul]=size(Factual);

ms=0;
for i=1:cul
    ms=ms+((max(Pareto_Front(:,i))-min(Pareto_Front(:,i)))^2);
end

% MS=sqrt(ms/cul);
MS=sqrt(ms);
end

function S=metric_of_spacing(Pareto_Front)

npf=size(Pareto_Front,1);
for i=1:npf
    D1=repmat(Pareto_Front(i,:),npf,1);
    D2=D1-Pareto_Front;
    for j=1:npf
        D3(j)=norm(D2(j,:));
    end
    D3(i)=[];
    
    [~ ,index]=min(D3);
    Distance(i)=sum(D2(index,:));
end

d_average=(1/npf)*sum(Distance);

total_distance=0;
for i=1:npf
    total_distance=total_distance+((Distance(i)-d_average)^2);
end

S=(sqrt(total_distance/(npf-1)));
end

function [pop, F]=NonDominatedSorting(pop)

nPop=numel(pop);
for i=1:nPop
    pop(i).DominationSet=[];
    pop(i).DominatedCount=0;
end

F{1}=[];
for i=1:nPop
    for j=i+1:nPop
        p=pop(i);
        q=pop(j);
        
%         if Dominates(p,q)
%             p.DominationSet=[p.DominationSet j];
%             q.DominatedCount=q.DominatedCount+1;
%         end
        
        if Dominates(q.Cost,q.G,p.Cost,p.G)
            q.DominationSet=[q.DominationSet i];
            p.DominatedCount=p.DominatedCount+1;
        end
        
        pop(i)=p;
        pop(j)=q;
    end
    
    if pop(i).DominatedCount==0
        F{1}=[F{1} i];
        pop(i).Rank=1;
    end
end

k=1;
while true
    
    Q=[];
    
    for i=F{k}
        p=pop(i);
        
        for j=p.DominationSet
            q=pop(j);
            q.DominatedCount=q.DominatedCount-1;
            
            if q.DominatedCount==0
                Q=[Q j];
                q.Rank=k+1;
            end
            pop(j)=q;
        end
    end
    
    if isempty(Q)
        break;
    end
    
    F{k+1}=Q;
    k=k+1;
end
end

function PlotCosts(pop)

global Factual

Costs=[pop.Cost];

plot(Costs(1,:),Costs(2,:),'r*','MarkerSize',8);
hold on
plot(Factual(:,1),Factual(:,2),'Color','blue','LineWidth',2);
legend('Obtained PF','Optimal PF');
hold off

xlabel('F_1 ');
ylabel('F_2 ');
grid on;

end

function RGD=Reverse_generational_distance(Pareto_Front,Factual)

% row=size(Pareto_Front,1);
%
% total_distance=0;
% for i=1:size(Factual,1)
%     D1=repmat(Factual(i,:),row,1);
%     D2=D1-Pareto_Front;
%     for j=1:row
%         D3(j)=norm(D2(j,:));
%     end
%     Distance=min(D3);
%     total_distance = total_distance + (Distance^2);
% end
%
% RGD=(1/size(Factual,1))*sqrt(total_distance);

% Most common used version

row=size(Pareto_Front,1);

total_distance=0;
for i=1:size(Factual,1)
    D1=repmat(Factual(i,:),row,1);
    D2=D1-Pareto_Front;
    for j=1:row
        D3(j)=norm(D2(j,:));
    end
    Distance=min(D3);
    total_distance = total_distance + Distance;
end

RGD=total_distance/size(Factual,1);

end

function [pop, F]=SortPopulation(pop)

% Sort Based on Crowding Distance
[~, CDSO]=sort([pop.CrowdingDistance],'descend');
pop=pop(CDSO);

% Sort Based on Rank
[~, RSO]=sort([pop.Rank]);
pop=pop(RSO);

% Update Fronts
Ranks=[pop.Rank];
MaxRank=max(Ranks);
F=cell(MaxRank,1);
for r=1:MaxRank
    F{r}=find(Ranks==r);
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