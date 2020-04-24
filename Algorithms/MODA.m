function rst = MODA(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)
    %___________________________________________________________________%
    %  Multi-Objective Dragonfly Algorithm (MODA) source codes demo     %
    %                           version 1.0                             %
    %                                                                   %
    %  Developed in MATLAB R2011b(7.13)                                 %
    %                                                                   %
    %  Author and programmer: Seyedali Mirjalili                        %
    %                                                                   %
    %         e-Mail: ali.mirjalili@gmail.com                           %
    %                 seyedali.mirjalili@griffithuni.edu.au             %
    %                                                                   %
    %       Homepage: http://www.alimirjalili.com                       %
    %                                                                   %
    %   Main paper:                                                     %
    %                                                                   %
    %   S. Mirjalili, Dragonfly algorithm: a new meta-heuristic         %
    %   optimization technique for solving single-objective, discrete,  %
    %   and multi-objective problems, Neural Computing and Applications %
    %   DOI: http://dx.doi.org/10.1007/s00521-015-1920-1                %
    %___________________________________________________________________%
    rand('state',sum(100*clock));
    tic;
%     clc;
%     clear;
%     close all;

    % Change these details with respect to your problem%%%%%%%%%%%%%%
    dim=nvar;
    lb=a';
    ub=b';
    obj_no=2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initial parameters of the MODA algorithm
    max_iter=nloop;
    N=nsol;
    ArchiveMaxSize=narchive;

    Archive_X=[];
    Archive_F=[];
    Archive_G=[];

    Archive_member_no=0;

    r=(ub-lb)/2;
    V_max=(ub(1)-lb(1))/10;

    Food_fitness=[];
    Food_G=[];
    Food_pos=[];

    Enemy_fitness=[];
    Enemy_G=[];
    Enemy_pos=[];

    X=initialization(N,dim,ub,lb);
    fitness=[];
    G=[];

    DeltaX=initialization(N,dim,ub,lb);
    iter=0;

    position_history=zeros(N,max_iter,dim);

    for iter=1:max_iter

        r=(ub-lb)/4+((ub-lb)*(iter/max_iter)*2);

        w=0.9-iter*((0.9-0.2)/max_iter);

        my_c=0.1-iter*((0.1-0)/(max_iter/2));
        if my_c<0
            my_c=0;
        end

        if iter<(3*max_iter/4)
            s=my_c;             % Seperation weight
            a=my_c;             % Alignment weight
            c=my_c;             % Cohesion weight
            f=2*rand;           % Food attraction weight
            e=my_c;             % Enemy distraction weight
        else
            s=my_c/iter;        % Seperation weight
            a=my_c/iter;        % Alignment weight
            c=my_c/iter;        % Cohesion weight
            f=2*rand;           % Food attraction weight
            e=my_c/iter;        % Enemy distraction weight
        end

        for i=1:N %Calculate all the objective values first
%             Particles_F(i,:)=ObjectiveFunction(X(:,i)');
            
            % Function Evaluation
            [fi,gi]=feval(fun,X(:,i));
            Particles_F(i,:)=fi';
            Particles_G(i,:)=gi';
            %
            %
            %
            
            if dominates(Particles_F(i,:),Particles_G(i,:),Food_fitness,Food_G)
                Food_fitness=Particles_F(i,:);
                Food_G=Particles_G(i,:);
                Food_pos=X(:,i);
            end

            if dominates(Enemy_fitness,Enemy_G,Particles_F(i,:),Particles_G(i,:))
                if all(X(:,i)<ub') && all( X(:,i)>lb')
                    Enemy_fitness=Particles_F(i,:);
                    Enemy_G=Particles_G(i,:);
                    Enemy_pos=X(:,i);
                end
            end
        end
        [Archive_X, Archive_F, Archive_G, Archive_member_no]=UpdateArchive(Archive_X, Archive_F, Archive_G, X, Particles_F, Particles_G, Archive_member_no);

        if Archive_member_no>ArchiveMaxSize
            Archive_mem_ranks=RankingProcess(Archive_F, Archive_G, ArchiveMaxSize, obj_no);
            [Archive_X, Archive_F, Archive_G, Archive_mem_ranks, Archive_member_no]=HandleFullArchive(Archive_X, Archive_F, Archive_G, Archive_member_no, Archive_mem_ranks, ArchiveMaxSize);
        else
            Archive_mem_ranks=RankingProcess(Archive_F, Archive_G, ArchiveMaxSize, obj_no);
        end

        Archive_mem_ranks=RankingProcess(Archive_F, Archive_G, ArchiveMaxSize, obj_no);

        % Chose the archive member in the least population area as foods
        % to improve coverage
        index=RouletteWheelSelection(1./Archive_mem_ranks);
        if index==-1
            index=1;
        end
        Food_fitness=Archive_F(index,:);
        Food_G=Archive_G(index,:);
        Food_pos=Archive_X(index,:)';

        % Chose the archive member in the most population area as enemies
        % to improve coverage
        index=RouletteWheelSelection(Archive_mem_ranks);
        if index==-1
            index=1;
        end
        Enemy_fitness=Archive_F(index,:);
        Enemy_G=Archive_G(index,:);
        Enemy_pos=Archive_X(index,:)';

        for i=1:N
            index=0;
            neighbours_no=0;

            clear Neighbours_V
            clear Neighbours_X
            % Find the neighbouring solutions
            for j=1:N
                Dist=distance(X(:,i),X(:,j));
                if (all(Dist<=r) && all(Dist~=0))
                    index=index+1;
                    neighbours_no=neighbours_no+1;
                    Neighbours_V(:,index)=DeltaX(:,j);
                    Neighbours_X(:,index)=X(:,j);
                end
            end

            % Seperation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Eq. (3.1)
            S=zeros(dim,1);
            if neighbours_no>1
                for k=1:neighbours_no
                    S=S+(Neighbours_X(:,k)-X(:,i));
                end
                S=-S;
            else
                S=zeros(dim,1);
            end

            % Alignment%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Eq. (3.2)
            if neighbours_no>1
                A=(sum(Neighbours_V')')/neighbours_no;
            else
                A=DeltaX(:,i);
            end

            % Cohesion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Eq. (3.3)
            if neighbours_no>1
                C_temp=(sum(Neighbours_X')')/neighbours_no;
            else
                C_temp=X(:,i);
            end

            C=C_temp-X(:,i);

            % Attraction to food%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Eq. (3.4)
            Dist2Attraction=distance(X(:,i),Food_pos(:,1));
            if all(Dist2Attraction<=r)
                F=Food_pos-X(:,i);
                iter;
            else
                F=0;
            end

            % Distraction from enemy%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Eq. (3.5)
            Dist=distance(X(:,i),Enemy_pos(:,1));
            if all(Dist<=r)
                E=Enemy_pos+X(:,i);
            else
                E=zeros(dim,1);
            end

            for tt=1:dim
                if X(tt,i)>ub(tt)
                    X(tt,i)=lb(tt);
                    DeltaX(tt,i)=rand;
                end
                if X(tt,i)<lb(tt)
                    X(tt,i)=ub(tt);
                    DeltaX(tt,i)=rand;
                end
            end


            if any(Dist2Attraction>r)
                if neighbours_no>1
                    for j=1:dim
                        DeltaX(j,i)=w*DeltaX(j,i)+rand*A(j,1)+rand*C(j,1)+rand*S(j,1);
                        if DeltaX(j,i)>V_max
                            DeltaX(j,i)=V_max;
                        end
                        if DeltaX(j,i)<-V_max
                            DeltaX(j,i)=-V_max;
                        end
                        X(j,i)=X(j,i)+DeltaX(j,i);
                    end

                else
                    X(:,i)=X(:,i)+Levy(dim)'.*X(:,i);
                    DeltaX(:,i)=0;
                end
            else    
                for j=1:dim
                    DeltaX(j,i)=s*S(j,1)+a*A(j,1)+c*C(j,1)+f*F(j,1)+e*E(j,1) + w*DeltaX(j,i);
                    if DeltaX(j,i)>V_max
                        DeltaX(j,i)=V_max;
                    end
                    if DeltaX(j,i)<-V_max
                        DeltaX(j,i)=-V_max;
                    end
                    X(j,i)=X(j,i)+DeltaX(j,i);
                end
            end

            Flag4ub=X(:,i)>ub';
            Flag4lb=X(:,i)<lb';
            X(:,i)=(X(:,i).*(~(Flag4ub+Flag4lb)))+ub'.*Flag4ub+lb'.*Flag4lb;

        end

%         display(['At the iteration ', num2str(iter), ' there are ', num2str(Archive_member_no), ' non-dominated solutions in the archive']);

        
        % Save Results
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
%     if obj_no==2
%         plot(Archive_F(:,1),Archive_F(:,2),'ko','MarkerSize',8,'markerfacecolor','k');
%     else
%         plot3(Archive_F(:,1),Archive_F(:,2),Archive_F(:,3),'ko','MarkerSize',8,'markerfacecolor','k');
%     end
%     legend('True PF','Obtained PF');
%     title('MODA');
    
    
    % Save Results
%     ppareto=Archive_X';
%     for i=1:size(ppareto,2)
%         [fpareto(:,i),gpareto(:,i)]=feval(fun,ppareto(:,i));
%     end
%     [ppareto,fpareto,gpareto,A]=pbil_selection0([],[],[],ppareto,fpareto,gpareto,narchive);
%     save(fout,'ppareto','fpareto','gpareto')
end

%%%%%%%%%%%%%%%%
% Sub-Function %
%%%%%%%%%%%%%%%%
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
            %
            %
            if (abs(F1-F2)<r(k))
                flag=flag+1;
            end
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

function o=Levy(d)

beta=3/2;
%Eq. (3.10)
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;
v=randn(1,d);
step=u./abs(v).^(1/beta);

% Eq. (3.9)
o=0.01*step;
end

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

X=X';
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
        o=0;
    else
        o=0;
    end
end

function o = distance(a,b)

for i=1:size(a,1)
    o(1,i)=sqrt((a(i)-b(i))^2);
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

function [bin,x,f,g] = moea_initialisation(fun,nsol,nvar,nbit,a,b)
%
% Randomly initiate the population, design variables 
% nvar=no. of design (decision) variables
% nbit is the number of binary cells in each variable
% nsol is a number of individuls
% 
% bin = a binary population
% x = a real number population of bin
% f = objective functions
% g = constraints
%

bin = round(rand(nbit*nvar,nsol));
x=bin2real(bin,a,b);
for i=1:nsol
    [f(:,i),g(:,i)]=feval(fun,x(:,i));
end
end

function x=bin2real(bin,a,b);% convert binary to decimal and then to real numbers
[m,n]=size(bin);
nvar=length(a);
nbit=m/nvar;

for i=1:n
    for j=1:nvar
        x(j,i)=bin2dec(bin((j-1)*nbit+1:j*nbit,i),a(j),b(j));
    end
end
end

function x=bin2dec(bin,a,b)

% 
% Transformation from binary string to real number
% with lowr limit a and upper limit b

n=max(size(bin));
trans=cumprod(2*ones(size(bin)))/2;
real1=sum(bin.*trans);

x=a+(real1*(b-a))/(2^n-1);
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