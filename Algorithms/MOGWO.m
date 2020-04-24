function rst = MOGWO(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)
    %___________________________________________________________________%
    %  Multi-Objective Grey Wolf Optimizer (MOGWO)                      %
    %  Source codes demo version 1.0                                    %
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
    %    S. Mirjalili, S. Saremi, S. M. Mirjalili, L. Coelho,           %
    %    Multi-objective grey wolf optimizer: A novel algorithm for     %
    %    multi-criterion optimization, Expert Systems with Applications,%
    %    in press, DOI: http://dx.doi.org/10.1016/j.eswa.2015.10.039    %       %
    %                                                                   %
    %___________________________________________________________________%

    % I acknowledge that this version of MOGWO has been written using
    % a large portion of the following code:

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  MATLAB Code for                                                  %
    %                                                                   %
    %  Multi-Objective Particle Swarm Optimization (MOPSO)              %
    %  Version 1.0 - Feb. 2011                                          %
    %                                                                   %
    %  According to:                                                    %
    %  Carlos A. Coello Coello et al.,                                  %
    %  "Handling Multiple Objectives with Particle Swarm Optimization," %
    %  IEEE Transactions on Evolutionary Computation, Vol. 8, No. 3,    %
    %  pp. 256-279, June 2004.                                          %
    %                                                                   %
    %  Developed Using MATLAB R2009b (Version 7.9)                      %
    %                                                                   %
    %  Programmed By: S. Mostapha Kalami Heris                          %
    %                                                                   %
    %         e-Mail: sm.kalami@gmail.com                               %
    %                 kalami@ee.kntu.ac.ir                              %
    %                                                                   %
    %       Homepage: http://www.kalami.ir                              %
    %                                                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rand('state',sum(100*clock));
tic
%     clear all
%     clc

    drawing_flag = 0;

    nVar=nvar;

    fobj = fun;

    % Lower bound and upper bound
    lb=a';
    ub=b';

    VarSize=[1 nVar];

    GreyWolves_num=nsol;
    MaxIt=nloop;  % Maximum Number of Iterations
    Archive_size=narchive;   % Repository Size

    alpha=0.1;  % Grid Inflation Parameter
    nGrid=10;   % Number of Grids per each Dimension
    beta=4; %=4;    % Leader Selection Pressure Parameter
    gamma=2;    % Extra (to be deleted) Repository Member Selection Pressure

    % Initialization

    GreyWolves=CreateEmptyParticle(GreyWolves_num);


    for i=1:GreyWolves_num
        GreyWolves(i).Velocity=0;
        GreyWolves(i).Position=zeros(1,nVar);
        for j=1:nVar
            GreyWolves(i).Position(1,j)=unifrnd(lb(j),ub(j),1);
        end
%         GreyWolves(i).Cost=fobj(GreyWolves(i).Position')';

        % Objective Function
        [fi,gi]=feval(fobj,GreyWolves(i).Position');
        GreyWolves(i).Cost=fi';
        GreyWolves(i).G=gi';
        %
        %
        GreyWolves(i).Best.Position=GreyWolves(i).Position;
        GreyWolves(i).Best.Cost=GreyWolves(i).Cost;
        GreyWolves(i).Best.G=GreyWolves(i).G;
    end

    GreyWolves=DetermineDomination(GreyWolves);

    Archive=GetNonDominatedParticles(GreyWolves);

    Archive_costs=GetCosts(Archive);
    G=CreateHypercubes(Archive_costs,nGrid,alpha);

    for i=1:numel(Archive)
        [Archive(i).GridIndex Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
    end

    % MOGWO main loop

    for iter=1:MaxIt
        a=2-iter*((2)/MaxIt);
        for i=1:GreyWolves_num

            clear rep2
            clear rep3

            % Choose the alpha, beta, and delta grey wolves
            Delta=SelectLeader(Archive,beta);
            Beta=SelectLeader(Archive,beta);
            Alpha=SelectLeader(Archive,beta);

            % If there are less than three solutions in the least crowded
            % hypercube, the second least crowded hypercube is also found
            % to choose other leaders from.
            if size(Archive,1)>1
                counter=0;
                for newi=1:size(Archive,1)
                    if sum(Delta.Position~=Archive(newi).Position)~=0
                        counter=counter+1;
                        rep2(counter,1)=Archive(newi);
                    end
                end
                if exist('rep2')~=0
                    Beta=SelectLeader(rep2,beta);
                end
            end

            % This scenario is the same if the second least crowded hypercube
            % has one solution, so the delta leader should be chosen from the
            % third least crowded hypercube.
            if size(Archive,1)>2 && exist('rep2')~=0
                counter=0;
                for newi=1:size(rep2,1)
                    if sum(Beta.Position~=rep2(newi).Position)~=0
                        counter=counter+1;
                        rep3(counter,1)=rep2(newi);
                    end
                end
                if exist('rep3')~=0
                    Alpha=SelectLeader(rep3,beta);
                end
            end

            % Eq.(3.4) in the paper
            c=2.*rand(1, nVar);
            % Eq.(3.1) in the paper
            D=abs(c.*Delta.Position-GreyWolves(i).Position);
            % Eq.(3.3) in the paper
            A=2.*a.*rand(1, nVar)-a;
            % Eq.(3.8) in the paper
            X1=Delta.Position-A.*abs(D);


            % Eq.(3.4) in the paper
            c=2.*rand(1, nVar);
            % Eq.(3.1) in the paper
            D=abs(c.*Beta.Position-GreyWolves(i).Position);
            % Eq.(3.3) in the paper
            A=2.*a.*rand()-a;
            % Eq.(3.9) in the paper
            X2=Beta.Position-A.*abs(D);


            % Eq.(3.4) in the paper
            c=2.*rand(1, nVar);
            % Eq.(3.1) in the paper
            D=abs(c.*Alpha.Position-GreyWolves(i).Position);
            % Eq.(3.3) in the paper
            A=2.*a.*rand()-a;
            % Eq.(3.10) in the paper
            X3=Alpha.Position-A.*abs(D);

            % Eq.(3.11) in the paper
            GreyWolves(i).Position=(X1+X2+X3)./3;

            % Boundary checking
            GreyWolves(i).Position=min(max(GreyWolves(i).Position,lb),ub);

%             GreyWolves(i).Cost=fobj(GreyWolves(i).Position')';
            
            % Objective Function
            [fi,gi]=feval(fobj,GreyWolves(i).Position');
            GreyWolves(i).Cost=fi';
            GreyWolves(i).G=gi';
            %
            %
        end

        GreyWolves=DetermineDomination(GreyWolves);
        non_dominated_wolves=GetNonDominatedParticles(GreyWolves);

        Archive=[Archive
            non_dominated_wolves];

        Archive=DetermineDomination(Archive);
        Archive=GetNonDominatedParticles(Archive);

        for i=1:numel(Archive)
            [Archive(i).GridIndex Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
        end

        if numel(Archive)>Archive_size
            EXTRA=numel(Archive)-Archive_size;
            Archive=DeleteFromRep(Archive,EXTRA,gamma);

            Archive_costs=GetCosts(Archive);
            G=CreateHypercubes(Archive_costs,nGrid,alpha);

        end

%         disp(['In iteration ' num2str(it) ': Number of solutions in the archive = ' num2str(numel(Archive))]);
%         save results

        % Results

        costs=GetCosts(GreyWolves);
        Archive_costs=GetCosts(Archive);

        if drawing_flag==1
            hold off
            plot(costs(1,:),costs(2,:),'k.');
            hold on
            plot(Archive_costs(1,:),Archive_costs(2,:),'rd');
            legend('Grey wolves','Non-dominated solutions');
            drawnow
        end
        
        % Save results
        ppareto=zeros(nvar,numel(Archive));
        fpareto=zeros(numel(Archive(1).Cost),numel(Archive));
        gpareto=zeros(numel(Archive(1).G),numel(Archive));
        for i=1:numel(Archive)
            ppareto(:,i)=Archive(i).Position';
            fpareto(:,i)=Archive(i).Cost';
            gpareto(:,i)=Archive(i).G';
        end
        [ppareto,fpareto,gpareto,~]=pbil_selection0([],[],[],ppareto,fpareto,gpareto,narchive);
        
        rst.ppareto{iter}=ppareto;
        rst.fpareto{iter}=fpareto;
        rst.gpareto{iter}=gpareto;
        rst.timestamp=datetime('now');
    end
    
%     for i=1:numel(Archive)
%         ppareto(:,i)=Archive(i).Position';
%     end
%     for i=1:numel(Archive)
%         [fpareto(:,i),gpareto(:,i)]=feval(fun,ppareto(:,i));
%     end
    
%     [ppareto,fpareto,gpareto,~]=pbil_selection0([],[],[],ppareto,fpareto,gpareto,narchive);
%     save(fout,'ppareto','fpareto','gpareto')
end

%%%%%%%%%%%%%%%%
% Sub-Function %
%%%%%%%%%%%%%%%%
function particle=CreateEmptyParticle(n)
    
    if nargin<1
        n=1;
    end

    empty_particle.Position=[];
    empty_particle.Velocity=[];
    empty_particle.Cost=[];
    empty_particle.Dominated=false;
    empty_particle.Best.Position=[];
    empty_particle.Best.Cost=[];
    empty_particle.GridIndex=[];
    empty_particle.GridSubIndex=[];
    
    particle=repmat(empty_particle,n,1);
    
end

function G=CreateHypercubes(costs,ngrid,alpha)

    nobj=size(costs,1);
    
    empty_grid.Lower=[];
    empty_grid.Upper=[];
    G=repmat(empty_grid,nobj,1);
    
    for j=1:nobj
        
        min_cj=min(costs(j,:));
        max_cj=max(costs(j,:));
        
        dcj=alpha*(max_cj-min_cj);
        
        min_cj=min_cj-dcj;
        max_cj=max_cj+dcj;
        
        gx=linspace(min_cj,max_cj,ngrid-1);
        
        G(j).Lower=[-inf gx];
        G(j).Upper=[gx inf];
        
    end

end

function rep=DeleteFromRep(rep,EXTRA,gamma)

    if nargin<3
        gamma=1;
    end

    for k=1:EXTRA
        [occ_cell_index occ_cell_member_count]=GetOccupiedCells(rep);

        p=occ_cell_member_count.^gamma;
        p=p/sum(p);

        selected_cell_index=occ_cell_index(RouletteWheelSelection(p));

        GridIndices=[rep.GridIndex];

        selected_cell_members=find(GridIndices==selected_cell_index);

        n=numel(selected_cell_members);

        selected_memebr_index=randi([1 n]);

        j=selected_cell_members(selected_memebr_index);
        
        rep=[rep(1:j-1); rep(j+1:end)];
    end
    
end

function pop=DetermineDomination(pop)

    npop=numel(pop);
    
    for i=1:npop
        pop(i).Dominated=false;
        for j=1:i-1
            if ~pop(j).Dominated
                if Dominates(pop(i).Cost,pop(i).G,pop(j).Cost,pop(j).G)
                    pop(j).Dominated=true;
                elseif Dominates(pop(j).Cost,pop(j).G,pop(i).Cost,pop(i).G)
                    pop(i).Dominated=true;
                    break;
                end
            end
        end
    end

end

function o=Dominates(f1,g1,f2,g2)

    if isstruct(f1)
        f1=f1.Cost;
    end
    if isstruct(g1)
        g1=g1.G;
    end
    if isstruct(f2)
        f2=f2.Cost;
    end
    if isstruct(g2)
        g2=g2.G;
    end
    
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

function costs=GetCosts(pop)

    nobj=numel(pop(1).Cost);
    costs=reshape([pop.Cost],nobj,[]);

end

function [Index SubIndex]=GetGridIndex(particle,G)

    c=particle.Cost;
    
    nobj=numel(c);
    ngrid=numel(G(1).Upper);
    
    str=['sub2ind(' mat2str(ones(1,nobj)*ngrid)];

    SubIndex=zeros(1,nobj);
    for j=1:nobj
        
        U=G(j).Upper;
        
        i=find(c(j)<U,1,'first');
        
        SubIndex(j)=i;
        
        str=[str ',' num2str(i)];
    end
    
    str=[str ');'];
    
    Index=eval(str);
    
end

function nd_pop=GetNonDominatedParticles(pop)

    ND=~[pop.Dominated];
    
    nd_pop=pop(ND);

end

function [occ_cell_index occ_cell_member_count]=GetOccupiedCells(pop)

    GridIndices=[pop.GridIndex];
    
    occ_cell_index=unique(GridIndices);
    
    occ_cell_member_count=zeros(size(occ_cell_index));

    m=numel(occ_cell_index);
    for k=1:m
        occ_cell_member_count(k)=sum(GridIndices==occ_cell_index(k));
    end
    
end

function i=RouletteWheelSelection(p)

    r=rand;
    c=cumsum(p);
    i=find(r<=c,1,'first');

end

function rep_h=SelectLeader(rep,beta)
    if nargin<2
        beta=1;
    end

    [occ_cell_index occ_cell_member_count]=GetOccupiedCells(rep);
    
    p=occ_cell_member_count.^(-beta);
    p=p/sum(p);
    
    selected_cell_index=occ_cell_index(RouletteWheelSelection(p));
    
    GridIndices=[rep.GridIndex];
    
    selected_cell_members=find(GridIndices==selected_cell_index);
    
    n=numel(selected_cell_members);
    
    selected_memebr_index=randi([1 n]);
    
    h=selected_cell_members(selected_memebr_index);
    
    rep_h=rep(h);
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