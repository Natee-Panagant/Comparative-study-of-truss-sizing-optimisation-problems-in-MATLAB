function rst = SHAMODE_WO(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)
    rand('state',sum(100*clock));
    tic
    
    min_nsol=4;%
    max_nsol=nsol;
    current_nsol=nsol;


    % Generate initial population
    [bin0,x0,f0,g0] = moea_initialisation(fun,nsol,nvar,nbit,a,b);
    A=zeros(nsol);
    A0=A;
    [ppareto,fpareto,gpareto,A]=pareto_sorting(x0,f0,g0,[],[],[],A,nsol);
    
    %Adaptive Processor (Initialize)
    adapt_dat=struct();
    adapt_dat.achieve_ratio=1.4;
    adapt_dat.achieve_size=adapt_dat.achieve_ratio*current_nsol;
    adapt_dat.x0a=[];
    adapt_dat.f0a=[];
    adapt_dat.g0a=[];
%     adapt_dat.xpbest_ratio=0.11;
    adapt_dat.mem_size=5;
    adapt_dat.mem_pos=1;
    adapt_dat.mem_sf=0.5*ones(1,adapt_dat.mem_size);
    adapt_dat.mem_cr=0.5*ones(1,adapt_dat.mem_size);
    
    % MAIN LOOP
    for iter=1:nloop
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    DE operatos & FSD resizing   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update population
        [x1,f1,g1,adapt_dat]=de_reproduct(fun,x0,f0,g0,ppareto,fpareto,gpareto,a,b,iter,nloop,adapt_dat);
        % Selection
        [x2,f2,g2,A2,nsort]=DEselection(x1,f1,g1,x0,f0,g0,A0);
        [ppareto,fpareto,gpareto,A]=pareto_sorting(x1,f1,g1,ppareto,fpareto,gpareto,A,narchive);

%         figure(1),clf,hold on,
%         plot(fpareto(1,:),fpareto(2,:),'s')
%         plot(f0(1,:),f0(2,:),'or')
%         pause(1)
        %
        %
        %
        
        x0=x2;f0=f2;g0=g2;A0=A2;
    
        %DE Adaptive Processor
        [adapt_dat]=DE_adaptive_processor(adapt_dat,x0,f0,g0,nsort,nsol);
        
        %Population sizing control
%         [x0,f0,g0,A0,adapt_dat]=popControl(x0,f0,g0,A0,adapt_dat,min_nsol,max_nsol,iter,nloop);
        
        % Save Results
        rst.ppareto{iter}=ppareto;
        rst.fpareto{iter}=fpareto;
        rst.gpareto{iter}=gpareto;
        rst.timestamp=datetime('now');
    end

    %Save Result
%     [ppareto,fpareto,gpareto,A]=pareto_sorting(ppareto,fpareto,gpareto,[],[],[],A0,nsol);
%     save(fout,'ppareto','fpareto','gpareto')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%
% SUB-Program %
%%%%%%%%%%%%%%%
function [x1,f1,g1,adapt_dat] = de_reproduct(fun,x0,f0,g0,ppareto,fpareto,gpareto,a,b,iter,nloop,adapt_dat)
    [nvar,nsol]=size(x0);
    f1=zeros(size(f0,1),nsol);
    g1=zeros(size(g0,1),nsol);
    
    % DE point generation for n1 parents (x1)
    % x0 = population from the previous generation
    % f0 = corresponding objective values
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % DE Mutation Operator %
    %%%%%%%%%%%%%%%%%%%%%%%%
    %Generate scaling factror(sf)
    sf=zeros(1,nsol);
    sf_ind=ceil(adapt_dat.mem_size*rand(1,nsol));
    sf_mid=adapt_dat.mem_sf(sf_ind);
    pos=true(size(sf));
    while sum(pos)>0
        sf(pos)=sf_mid(pos)+0.1*tan(pi*(rand(1,sum(pos))-0.5));
        pos=sf<=0;
    end
    sf=min(1,sf);
    adapt_dat.sf=sf;
    x0a=[x0,adapt_dat.x0a];
    
    Npbest=size(ppareto,2);%max(2,round(adapt_dat.xpbest_ratio*nsol)); %Select at least 2 pbest
    pbi=ceil(Npbest*rand(1,nsol));
    xpbest=ppareto(:,pbi);
    
    [i1,i2]=ind_gen(nsol,size(x0a,2));
    xm=x0+sf(ones(nvar,1),:).*(xpbest-x0+x0(:,i1)-x0a(:,i2));

    xm=LUbound(xm,x0,[a,b]);
    % Updatae with Bubble-net attacking method
    % of Whale Optimization Algorithm (WOA)
    bb=1;
    pbi2=ceil(Npbest*rand(1,nsol));
    xpbest2=ppareto(:,pbi2);
    for wi=1:nsol
        if rand<0.5
            l=2*rand-1;
            xm(:,wi)=abs(xpbest2(:,wi)-xm(:,wi))*exp(bb*l)*cos(2*pi*l)+xpbest2(:,wi);
        end
    end
    xm=LUbound(xm,x0,[a,b]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % DE Crossover Operation %
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %Generate crossover ratio
    cr_ind=ceil(adapt_dat.mem_size*rand(1,nsol));
    cr_mid=adapt_dat.mem_cr(cr_ind);
    cr=normrnd(cr_mid, 0.1);
    cr=max(0,min(1,cr));
    adapt_dat.cr=cr;
    
    %Binomial Crossover
    x1=x0;
    mask=rand(nvar,nsol)<=cr(ones(nvar,1),:);%mask is used to indicate which elements of xr comes from xm
    rows=floor(rand(1,nsol)*nvar)+1; % choose one position where the element of xr always come from xm
    cols=(1:nsol);
    jrand=sub2ind([nvar,nsol],rows,cols);
    mask(jrand)=true;
    x1(mask)=xm(mask);
    
    for i=1:size(x1,2)
        [f1(:,i),g1(:,i)]=feval(fun,x1(:,i));
    end
end

function [adapt_dat] = DE_adaptive_processor(adapt_dat,x0,f0,g0,nsort,nsol)
    % Get success update index (Parent of success offspring)
    % index 1:nsol = Parent, (nsol+1):(2*nsol) = Offspring
    sind=nsort(1:nsol); % Success Offspring = Offspring in first 1:nsol rank
    sind=sind(sind>nsol)-nsol; % Offspring index - nsol = Parent index
    %
    adapt_dat.achieve_size=round(adapt_dat.achieve_ratio*nsol);

    % Update adaptive population
    x0a=[adapt_dat.x0a,x0(:,sind)];
    f0a=[adapt_dat.f0a,f0(:,sind)];
    g0a=[adapt_dat.g0a,g0(:,sind)];
    [~,Uind]=unique(x0a','rows');
    x0a=x0a(:,Uind);
    f0a=f0a(:,Uind);
    g0a=g0a(:,Uind);
    
    if size(adapt_dat.f0a,2)<=adapt_dat.achieve_size
        adapt_dat.x0a=x0a;
        adapt_dat.f0a=f0a;
        adapt_dat.g0a=g0a;
    else
        ind=randperm(size(adapt_dat.f0a,2));
        ind=ind(1:adapt_dat.achieve_size);
        adapt_dat.x0a=x0a(:,ind);
        adapt_dat.f0a=f0a(:,ind);
        adapt_dat.g0a=g0a(:,ind);
    end

    %Update scaling factror and crossover memory
    if sum(sind)>0 && adapt_dat.mem_cr(adapt_dat.mem_pos)~=-1
        %Update sf
        sssf=adapt_dat.sf(sind);
        lmeanSF=sum(sssf.^2)/sum(sssf); % Lehmer mean
        adapt_dat.mem_sf(adapt_dat.mem_pos,:)=lmeanSF;

        %Update cr
        sscr=adapt_dat.cr(sind);
        if max(sscr)==0 || adapt_dat.mem_cr(adapt_dat.mem_pos)==-1
            adapt_dat.mem_cr(adapt_dat.mem_pos)=-1;
        else
            lmeanCR=sum(sscr.^2)/sum(sscr); % Lehmer mean
            adapt_dat.mem_cr(adapt_dat.mem_pos)=lmeanCR;
        end
    end
    adapt_dat.mem_sf=max(0,min(1,adapt_dat.mem_sf));
    adapt_dat.mem_cr=max(0,min(1,adapt_dat.mem_cr));
    
    %Shift memory position
    adapt_dat.mem_pos=adapt_dat.mem_pos+1;
    if adapt_dat.mem_pos>adapt_dat.mem_size
        adapt_dat.mem_pos=1;
    end
end

% function [x0,f0,g0,A0,adapt_dat]=popControl(x0,f0,g0,A0,adapt_dat,min_nsol,max_nsol,iter,nloop)
%     current_nsol=round(max_nsol-(max_nsol-min_nsol)*iter/nloop);
%     nsol=size(f0,2);
%     K=nsol;
%     if current_nsol<K
%         Rind=[];
%         B=ndlevel_sort(A0);
%         NL=max(B);
%         for l=NL:-1:1
%             Lind=find(B==l);
%             N=numel(Lind);
%             if K-N>=current_nsol
%                 Rind=[Rind,Lind];
%                 K=K-N;
%             else
%                 ind=randperm(N,K-current_nsol);
%                 Lind=Lind(ind);
%                 Rind=[Rind,Lind];
%             end
%             if nsol-numel(Rind)==current_nsol
%                 break;
%             elseif nsol-numel(Rind)<current_nsol
%                 error('please check');
%             end
%         end
%         x0(:,Rind)=[];
%         f0(:,Rind)=[];
%         g0(:,Rind)=[];
%         A0(:,Rind)=[];
%         A0(Rind,:)=[];
%         adapt_dat.sf(:,Rind)=[];
%         adapt_dat.cr(:,Rind)=[];
%     end
% end

function [i1,i2] = ind_gen(nsol,nsolA)
    i0=1:nsol;

    i1=max(1,ceil(rand(1,nsol)*nsol));
    dpos=i1==i0;
    iter=0;
    while sum(dpos)>0
        iter=iter+1;
        i1(dpos)=max(1,ceil(rand(1,sum(dpos))*nsol));
        dpos=i1==i0;
        if iter>999
            error('check');
        end
    end
    
    i2=max(1,ceil(rand(1,nsol)*nsolA));
    dpos=i2==i0 | i2==i1;
    iter=0;
    while sum(dpos)>0
        iter=iter+1;
        i2(dpos)=max(1,ceil(rand(1,sum(dpos))*nsolA));
        dpos=i2==i0 | i2==i1;
        if iter>999
            error('check');
        end
    end
end

function [vi] = LUbound(vi,pop,lu)
    bound_ratio=2;
    lb=lu(:,ones(1,size(vi,2)));
    ub=lu(:,2*ones(1,size(vi,2)));
    lind=vi<lb;
    vi(lind)=(1/bound_ratio)*((bound_ratio-1)*lb(lind)+pop(lind));
    uind=vi>ub;
    vi(uind)=(1/bound_ratio)*((bound_ratio-1)*ub(uind)+pop(uind));
end

%%%%%%%%%%%%%%%%%%%%%%%
% MODEMO sub function %
%%%%%%%%%%%%%%%%%%%%%%%
function [pop3,f3,g3,A3,nsort]=DEselection(pop2,f2,g2,pop1,f1,g1,A1)
    %
    % Selection procedure of NSGAII
    % 
    pop=[pop1 pop2];
    f=[f1 f2];
    g=[g1 g2];
    
    [m0,n0]=size(A1);
    [m1,n1]=size(pop);
    [m2,n2]=size(f);
    [m3,n3]=size(f2);
    %%%%%%%%%%%%%%%%%%%%%%
    A=zeros(size(pop,2));
    for i=1:n0
        fi=f(:,i);
        gi=g(:,i);
        A(i,i)=0;
        for j=(i+1):n1
            fj=f(:,j);
            gj=g(:,j);
            %%%%%%%%%%%%%%%%%%%%%%%%%
            [p_count1,p_count2]=fdominated(fi,gi,fj,gj);
            A(i,j)=p_count1;
            A(j,i)=p_count2;
        end
    end
    
    B=ndlevel_sort(A);% level of being dominated, 1 = non-dominated
    
    [B,nsort]=sort(B);

    pop=pop(:,nsort);
    f=f(:,nsort);
    g=g(:,nsort);
    A=A(nsort,nsort);

    nlevel=zeros(1,max(B));
    for i=1:n1;nlevel(B(i))=nlevel(B(i))+1;end
    nlevel2=cumsum(nlevel);

    ncheck=nlevel(1);icheck=1;
    while ncheck < n3
        icheck=icheck+1;
        ncheck=ncheck+nlevel(icheck);
    end

    if ncheck == n3
        pop3=pop(:,1:n3);
        f3=f(:,1:n3);
        g3=g(:,1:n3);
        A3=A(1:n3,1:n3);
    elseif ncheck > n3 && icheck==1
        nc2=1:nlevel2(1);

        popc=pop(:,nc2);
        fc=f(:,nc2);
        gc=g(:,nc2);
        Ac=A(nc2,nc2);

        iselect=farchive(fc,n3);

        pop3=popc(:,iselect);
        f3=fc(:,iselect);
        g3=gc(:,iselect);
        A3=Ac(iselect,iselect);
    else
        nn1=nlevel2(icheck-1);
        nn2=n3-nn1;

        nc1=1:nn1;
        nc2=(nn1+1):nlevel2(icheck);

        popc=pop(:,nc2);
        fc=f(:,nc2);
        gc=g(:,nc2);
        Ac=A(nc2,nc2);

        iselect=farchive(fc,nn2);

        pop3=pop(:,iselect);
        f3=f(:,iselect);
        g3=g(:,iselect);
        A3=A(iselect,iselect);

        nc3=[nc1,nn1+iselect];

        pop3=pop(:,nc3);
        f3=f(:,nc3);
        g3=g(:,nc3);
        A3=A(nc3,nc3);
    end
end
function [pareto1,fpareto1,gpareto1,A]=pareto_sorting(x1,f1,g1,pareto,fpareto,gpareto,A,narchive)
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
        for j=(n0+1):n1
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
    for i=(n0+1):n1
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    if nndm > narchive
        nsl=farchive(fpareto1,narchive);
        pareto1=pareto1(:,nsl);
        fpareto1=fpareto1(:,nsl);
        gpareto1=gpareto1(:,nsl);
        A=A(nsl,nsl);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seletion used when the number of non-dominated solutions
% exceeds the archive size
function iselect=farchive(f,narchive)
    iselect=farchive72(f,narchive);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B=ndlevel_sort(A)
    [m,n]=size(A);
    b=sum(A,1);
    [b,ns]=sort(b);

    ilevel=1;
    B(1)=ilevel;

    for i=2:n
        if b(i)==b(i-1)
            B(i)=ilevel;
        else
            ilevel=ilevel+1;
            B(i)=ilevel;
        end
    end
    B(ns)=B;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fc=fcrownding(f)
    [m,n]=size(f);
    delf=max(f,[],2)-min(f,[],2);

    I=zeros(m,n);
    for i=1:m
        [fs,nsort]=sort(f(i,:));
        I(i,nsort(1))=1e10;I(i,nsort(n))=1e10;
        for j=2:(n-1)
            I(i,nsort(j))=abs(f(i,j+1)-f(i,j-1))/delf(i);
        end
    end     
    fc=sum(I,1);
end

function iselect=farchive72(f,narchive)
[m,n]=size(f);

fmax=max(f,[],2);[fmin nmin]=min(f,[],2);
fdel=max(fmax-fmin,1e-5);
iselect0=[nmin'];

for i=1:n
    f(:,i)=(f(:,i)-fmin)./fdel;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n
    D(i,i)=0;
    fi=f(:,i);
    for j=i:n
        fj=f(:,j);
        dij=sqrt(sum((fi-fj).*(fi-fj)));
        D(i,j)=dij;D(j,i)=dij;       
    end
end

A=[(1:n)' zeros(n,n) ones(n,1)];
ngroup=size(A,1);
k=0;
while ngroup > narchive
    Dg=[];
    for i=1:size(A,1)
        Idi=A(i,1:A(i,n+2));
        Dg(i,i)=1e20;
        for j=(i+1):size(A,1)
            Idj=A(j,1:A(j,n+2));
            Dgij=D(Idi,Idj);
            Dg(i,j)=min(min(Dgij));
            Dg(j,i)=Dg(i,j);
        end
    end
    [Dgs1,ns1]=min(Dg,[],2);
    [Dgs2,ns2]=min(Dgs1);
    I=ns1(ns2);
    J=ns2;
    nA=A(I,n+2)+A(J,n+2);
    A(I,1:nA)=[A(I,1:A(I,n+2)) A(J,1:A(J,n+2))];
    A(I,n+2)=nA;
    A(J,:)=[];
    ngroup=size(A,1);
end
% figure(1),clf,hold on
%     plot(f(1,:),f(2,:),'s')
%     for k=1:size(A,1)
%         for l=1:A(k,n+2)
%             text(f(1,A(k,l)),f(2,A(k,l)),['\bf' num2str(k)],'fontsize',12);
%         end
%     end
%     pause
for k=1:size(A,1)
    islk=A(k,1:A(k,n+2));
    rG=mean(f(:,islk),2);dd=[];
    sml=0;
    for l=1:length(islk)
        fl=f(:,A(k,l));
        dd(l)=sqrt(sum((fl-rG).*(fl-rG)));
        for ii=1:m
            if iselect0(ii)==A(k,l)
                sml=1;Isml=iselect0(ii);
            end
        end
    end
    if sml==1
        iselect(k)=Isml;
    else
        [dmin,nmin]=min(dd);
        iselect(k)=A(k,nmin);
    end
end
% plot(f(1,iselect),f(2,iselect),'xr')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%