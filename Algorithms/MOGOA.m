function rst = MOGOA(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)
    %_________________________________________________________________________________
    %  Multi-objective Grasshopper Optimization Algorithm (MOGOA) source codes version 1.0
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
    %   S. Z. Mirjalili, S. Mirjalili, S. Saremi, H. Fatis, H. Aljarah, 
    %   Grasshopper optimization algorithm for multi-objective optimization problems, 
    %   Applied Intelligence, 2017, DOI: http://dx.doi.org/10.1007/s10489-017-1019-8
    %____________________________________________________________________________________
    rand('state',sum(100*clock));
    tic
    % clc;
    % clear;
    % close all;

    % Change these details with respect to your problem%%%%%%%%%%%%%%
    dim=nvar;
    lb=a';
    ub=b';
    obj_no=2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    flag=0;
    if (rem(dim,2)~=0)
        dim = dim+1;
        ub = [ub, 1];
        lb = [lb, 0];
        flag=1;
    end


    max_iter=nloop;
    N=nsol;
    ArchiveMaxSize=narchive;

    Archive_X=[];
    Archive_F=[];
    Archive_G=[];
    
    Archive_member_no=0;

    %Initialize the positions of artificial whales
    GrassHopperPositions=initialization(N,dim,ub,lb);

    TargetPosition=[];
    TargetFitness=[];
    TargetG=[];
    
    cMax=1;
    cMin=0.00004;
    %calculate the fitness of initial grasshoppers

    for iter=1:max_iter
        for i=1:N

            Flag4ub=GrassHopperPositions(:,i)>ub';
            Flag4lb=GrassHopperPositions(:,i)<lb';
            GrassHopperPositions(:,i)=(GrassHopperPositions(:,i).*(~(Flag4ub+Flag4lb)))+ub'.*Flag4ub+lb'.*Flag4lb;
    %         GrassHopperFitness(i,:)=ObjectiveFunction(GrassHopperPositions(:,i)');

            % Function Evaluation
            [fi,gi]=feval(fun,GrassHopperPositions(1:nvar,i));
            GrassHopperFitness(i,:)=fi';
            GrassHopperG(i,:)=gi';
            %
            %
            %
            if dominates(GrassHopperFitness(i,:),GrassHopperG(i,:),TargetFitness,TargetG)
                TargetFitness=GrassHopperFitness(i,:);
                TargetG=GrassHopperG(i,:);
                TargetPosition=GrassHopperPositions(:,i);
            end

        end
        
        [Archive_X, Archive_F, Archive_G, Archive_member_no]=UpdateArchive(Archive_X, Archive_F, Archive_G, GrassHopperPositions, GrassHopperFitness, GrassHopperG, Archive_member_no);
        
        if Archive_member_no>ArchiveMaxSize
            Archive_mem_ranks=RankingProcess(Archive_F, Archive_G, ArchiveMaxSize, obj_no);
            [Archive_X, Archive_F, Archive_G, Archive_mem_ranks, Archive_member_no]=HandleFullArchive(Archive_X, Archive_F, Archive_G, Archive_member_no, Archive_mem_ranks, ArchiveMaxSize);
        else
            Archive_mem_ranks=RankingProcess(Archive_F, Archive_G, ArchiveMaxSize, obj_no);
        end

        Archive_mem_ranks=RankingProcess(Archive_F, Archive_G, ArchiveMaxSize, obj_no);
        index=RouletteWheelSelection(1./Archive_mem_ranks);
        if index==-1
            index=1;
        end
        TargetFitness=Archive_F(index,:);
        TargetG=Archive_G(index,:);
        TargetPosition=Archive_X(index,:)';

        c=cMax-iter*((cMax-cMin)/max_iter); % Eq. (3.8) in the paper

        for i=1:N

            temp= GrassHopperPositions;

            for k=1:2:dim
                S_i=zeros(2,1);
                for j=1:N
                    if i~=j
                        Dist=distance(temp(k:k+1,j), temp(k:k+1,i));
                        r_ij_vec=(temp(k:k+1,j)-temp(k:k+1,i))/(Dist+eps);
                        xj_xi=2+rem(Dist,2);

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % Eq. (3.2) in the paper 
                        s_ij=((ub(k:k+1)' - lb(k:k+1)') .*c/2)*S_func(xj_xi).*r_ij_vec;
                        S_i=S_i+s_ij;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    end
                end
                S_i_total(k:k+1, :) = S_i;

            end

            X_new=c*S_i_total'+(TargetPosition)'; % Eq. (3.7) in the paper
            GrassHopperPositions_temp(i,:)=X_new';
        end
        % GrassHopperPositions
        GrassHopperPositions=GrassHopperPositions_temp';

%         display(['At the iteration ', num2str(iter), ' there are ', num2str(Archive_member_no), ' non-dominated solutions in the archive']);
        
        % Save results
        ppareto=Archive_X(:,1:nvar)';
        fpareto=Archive_F';
        gpareto=Archive_G';
        [ppareto,fpareto,gpareto,~]=pbil_selection0([],[],[],ppareto,fpareto,gpareto,narchive);
        
        rst.ppareto{iter}=ppareto;
        rst.fpareto{iter}=fpareto;
        rst.gpareto{iter}=gpareto;
        rst.timestamp=datetime('now');
    end


%     if (flag==1)
%         TargetPosition = TargetPosition(1:dim-1);
%     end

%     figure
% 
%     Draw_ZDT1();
% 
%     hold on
% 
%     plot(Archive_F(:,1),Archive_F(:,2),'ro','MarkerSize',8,'markerfacecolor','k');
% 
%     legend('True PF','Obtained PF');
%     title('MOGOA');
% 
%     set(gcf, 'pos', [403   466   230   200])

    % Save Results
%     ppareto=Archive_X(:,1:nvar)';
%     for i=1:size(ppareto,2)
%         [fpareto(:,i),gpareto(:,i)]=feval(fun,ppareto(:,i));
%     end
%     [ppareto,fpareto,gpareto,A]=pbil_selection0([],[],[],ppareto,fpareto,gpareto,narchive);
%     save(fout,'ppareto','fpareto','gpareto')
end

%%%%%%%%%%%%%%%%
% Sub-Function %
%%%%%%%%%%%%%%%%
function d = distance(a,b)
% DISTANCE - computes Euclidean distance matrix
%
% E = distance(A,B)
%
%    A - (DxM) matrix 
%    B - (DxN) matrix
%
% Returns:
%    E - (MxN) Euclidean distances between vectors in A and B
%
%
% Description : 
%    This fully vectorized (VERY FAST!) m-file computes the 
%    Euclidean distance between two vectors by:
%
%                 ||A-B|| = sqrt ( ||A||^2 + ||B||^2 - 2*A.B )
%
% Example : 
%    A = rand(400,100); B = rand(400,200);
%    d = distance(A,B);

% Author   : Roland Bunschoten
%            University of Amsterdam
%            Intelligent Autonomous Systems (IAS) group
%            Kruislaan 403  1098 SJ Amsterdam
%            tel.(+31)20-5257524
%            bunschot@wins.uva.nl
% Last Rev : Oct 29 16:35:48 MET DST 1999
% Tested   : PC Matlab v5.2 and Solaris Matlab v5.3
% Thanx    : Nikos Vlassis

% Copyright notice: You are free to modify, extend and distribute 
%    this code granted that the author of the original code is 
%    mentioned as the original author of the code.

% if (nargin ~= 2)
%    error('Not enough input arguments');
% end
% 
% if (size(a,1) ~= size(b,1))
%    error('A and B should be of same dimensionality');
% end
% 
% aa=sum(a.*a,1); bb=sum(b.*b,1); ab=a'*b; 
% d = sqrt(abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab));
d=sqrt((a(1)-b(1))^2+(a(2)-b(2))^2);
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

function TPF=Draw_ZDT1()

% TPF is the true Pareto optimal front
addpath('ZDT_set')

ObjectiveFunction=@(x) ZDT1(x);
x=0:0.01:1;
for i=1:size(x,2)
    TPF(i,:)=ObjectiveFunction([x(i) 0 0 0]);
end
line(TPF(:,1),TPF(:,2));
title('ZDT1')

xlabel('f1')
ylabel('f2')
box on

fig=gcf;
set(findall(fig,'-property','FontName'),'FontName','Garamond')
set(findall(fig,'-property','FontAngle'),'FontAngle','italic')
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

function Positions=initialization(SearchAgents_no,dim,ub,lb)

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
        Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end

Positions=Positions';
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

function choice = RouletteWheelSelection(weights)
  accumulation = cumsum(weights);
  p = rand() * accumulation(end);
  chosen_index = -1;
  for index = 1 : length(accumulation)
    if (accumulation(index) > p)
      chosen_index = index;
      break;
    end
  end
  choice = chosen_index;
end

function o=S_func(r)
F=0.5;
L=1.5;
o=F*exp(-r/L)-exp(-r); % Eq. (3.3) in the paper
end

function [Archive_X_updated, Archive_F_updated, Archive_G_updated, Archive_member_no]=UpdateArchive(Archive_X, Archive_F, Archive_G, Particles_X, Particles_F, Particles_G, Archive_member_no)
    [Archive_X_temp,uind]=unique([Archive_X ; Particles_X'],'rows');
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