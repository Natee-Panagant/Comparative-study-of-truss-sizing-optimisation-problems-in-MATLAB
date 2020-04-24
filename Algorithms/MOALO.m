function rst = MOALO(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)
%_________________________________________________________________________%
%  Multi-Objective Ant Lion Optimizer (MALO) source codes demo            %
%                           version 1.0                                   %
%                                                                         %
%  Developed in MATLAB R2011b(7.13)                                       %
%                                                                         %
%  Author and programmer: Seyedali Mirjalili                              %
%                                                                         %
%         e-Mail: ali.mirjalili@gmail.com                                 %
%                 seyedali.mirjalili@griffithuni.edu.au                   %
%                                                                         %
%       Homepage: http://www.alimirjalili.com                             %
% Paper: Mirjalili, Seyedali, Pradeep Jangir, and Shahrzad Saremi.        %
% Multi-objective ant lion optimizer: a multi-objective optimization      %
%  algorithm for solving engineering problems." Applied Intelligence      %
% (2016): 1-17, DOI: http://dx.doi.org/10.1007/s10489-016-0825-8          %
%_________________________________________________________________________%

% clc;
% clear;
% close all;
tic
rand('state',sum(100*clock));
    % Change these details with respect to your problem%%%%%%%%%%%%%%
%     ObjectiveFunction=@ZDT1;
    dim=nvar;
    lb=a';
    ub=b';
    obj_no=2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initial parameters of the MODA algorithm
    max_iter=nloop;
    N=nsol;
    ArchiveMaxSize=narchive;

    Archive_X=[];%
    Archive_F=[];%
    Archive_G=[];%
    
    Archive_member_no=0;

    r=(ub-lb)/2;
    V_max=(ub(1)-lb(1))/10;

    Elite_fitness=[];
    Elite_constraint=[];
    Elite_position=[];

    Ant_Position=initialization(N,dim,ub,lb);
    fitness=zeros(N,2);

    V=initialization(N,dim,ub,lb);
    iter=0;

    position_history=zeros(N,max_iter,dim);

    for iter=1:max_iter

        for i=1:N %Calculate all the objective values first
%             Particles_F(i,:)=ObjectiveFunction(Ant_Position(:,i)');
            %
            %
            % Objective Function
            [fi,gi]=feval(fun,Ant_Position(:,i));
            Particles_F(i,:)=fi';
            Particles_G(i,:)=gi';
            %
            %
            if dominates(Particles_F(i,:),Particles_G(i,:),Elite_fitness,Elite_constraint)
                Elite_fitness=Particles_F(i,:);
                Elite_constraint=Particles_G(i,:);
                Elite_position=Ant_Position(:,i);
            end
        end

        [Archive_X, Archive_F, Archive_G, Archive_member_no]=UpdateArchive(Archive_X, Archive_F, Archive_G, Ant_Position, Particles_F, Particles_G, Archive_member_no);

        if Archive_member_no>ArchiveMaxSize
            Archive_mem_ranks=RankingProcess(Archive_F, Archive_G, ArchiveMaxSize, obj_no);
            [Archive_X, Archive_F, Archive_G, Archive_mem_ranks, Archive_member_no]=HandleFullArchive(Archive_X, Archive_F, Archive_G, Archive_member_no, Archive_mem_ranks, ArchiveMaxSize);
        else
            Archive_mem_ranks=RankingProcess(Archive_F, Archive_G, ArchiveMaxSize, obj_no);
        end
        
        Archive_mem_ranks=RankingProcess(Archive_F, Archive_G, ArchiveMaxSize, obj_no);

        % Chose the archive member in the least population area as arrtactor
        % to improve coverage
        index=RouletteWheelSelection(1./Archive_mem_ranks);
        if index==-1
            index=1;
        end
        Elite_fitness=Archive_F(index,:);
        Elite_constraint=Archive_G(index,:);
        Elite_position=Archive_X(index,:)';

        Random_antlion_fitness=Archive_F(1,:);
        Random_antlion_position=Archive_X(1,:)';

        for i=1:N

            index=0;
            neighbours_no=0;

            RA=Random_walk_around_antlion(dim,max_iter,lb,ub, Random_antlion_position',iter);

            [RE]=Random_walk_around_antlion(dim,max_iter,lb,ub, Elite_position',iter);

            Ant_Position(:,i)=(RE(iter,:)'+RA(iter,:)')/2;



            Flag4ub=Ant_Position(:,i)>ub';
            Flag4lb=Ant_Position(:,i)<lb';
            Ant_Position(:,i)=(Ant_Position(:,i).*(~(Flag4ub+Flag4lb)))+ub'.*Flag4ub+lb'.*Flag4lb;

        end
%         display(['At the iteration ', num2str(iter), ' there are ', num2str(Archive_member_no), ' non-dominated solutions in the archive']);
        
        %save results
        [ppareto,fpareto,gpareto,~]=pbil_selection0([],[],[],Archive_X',Archive_F',Archive_G',narchive);
        
        rst.ppareto{iter}=ppareto;
        rst.fpareto{iter}=fpareto;
        rst.gpareto{iter}=gpareto;
        rst.timestamp=datetime('now');
    end

%     figure
% 
%     Draw_ZDT1();
% 
%     hold on
% 
%     plot(Archive_F(:,1),Archive_F(:,2),'ko','MarkerSize',8,'markerfacecolor','k');
% 
%     legend('True PF','Obtained PF');
%     title('MALO');
% 
%     set(gcf, 'pos', [403   466   230   200])
%     ppareto=Archive_X';
%     for i=1:size(Archive_X,1)
%         [fpareto(:,i),gpareto(:,i)]=feval(fun,ppareto(:,i));
%     end
    
%     [ppareto,fpareto,gpareto,~]=pbil_selection0([],[],[],ppareto,fpareto,gpareto,narchive);
%     save(fout,'ppareto','fpareto','gpareto');
end

%%%%%%%%%%%%%%%%
% Sub-Function %
%%%%%%%%%%%%%%%%
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

function [RWs]=Random_walk_around_antlion(Dim,max_iter,lb, ub,antlion,current_iter)
    if size(lb,1) ==1 && size(lb,2)==1 %Check if the bounds are scalar
        lb=ones(1,Dim)*lb;
        ub=ones(1,Dim)*ub;
    end

    if size(lb,1) > size(lb,2) %Check if boundary vectors are horizontal or vertical
        lb=lb';
        ub=ub';
    end

    I=1; % I is the ratio in Equations (2.10) and (2.11)

    if current_iter>max_iter/10
        I=1+100*(current_iter/max_iter);
    end

    if current_iter>max_iter/2
        I=1+1000*(current_iter/max_iter);
    end

    if current_iter>max_iter*(3/4)
        I=1+10000*(current_iter/max_iter);
    end

    if current_iter>max_iter*(0.9)
        I=1+100000*(current_iter/max_iter);
    end

    if current_iter>max_iter*(0.95)
        I=1+1000000*(current_iter/max_iter);
    end


    % Dicrease boundaries to converge towards antlion
    lb=lb/(I); % Equation (2.10) in the paper
    ub=ub/(I); % Equation (2.11) in the paper

    % Move the interval of [lb ub] around the antlion [lb+anlion ub+antlion]
    if rand<0.5
        lb=lb+antlion; % Equation (2.8) in the paper
    else
        lb=-lb+antlion;
    end

    if rand>=0.5
        ub=ub+antlion; % Equation (2.9) in the paper
    else
        ub=-ub+antlion;
    end

    % This function creates n random walks and normalize accroding to lb and ub
    % vectors 
    for i=1:Dim
        X = [0 cumsum(2*(rand(max_iter,1)>0.5)-1)']; % Equation (2.1) in the paper
        %[a b]--->[c d]
        a=min(X);
        b=max(X);
        c=lb(i);
        d=ub(i);      
        X_norm=((X-a).*(d-c))./(b-a)+c; % Equation (2.7) in the paper
        RWs(:,i)=X_norm;
    end
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