function rst = RPBILDE(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)
%%%%%%%%%%%%%%%%% MAIN PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Initialisation Process %%%%%%%%%%%%
rand('state',sum(100*clock));
tic

LR0=0.25;% Learning rate
Mut_prob=0.05;%mutation probability
Mut_shft=0.20;%mutation shift
ninterval=40;%max([10 nvar]);
ntray=ceil(nsol/5);
pc=0.7;% crossover probability
F=0.8;% Scaling factor for differentail evolution (DE) operator
CR=0.5;% probability of chosing element from offspring in crossover

narchive1=narchive;
narchive2=narchive;
% a=0*ones(nvar,1);
% b=1*ones(nvar,1);

% for minimisation
[bin0,x0,f0,g0]=moea_initialisation(fun,nsol,nvar,nbit,a,b);
ppareto=[];fpareto=[];gpareto=[];
[ppareto,fpareto,gpareto,A]=pbil_selection0(x0,f0,g0,ppareto,fpareto,gpareto,narchive1);

for i=1:ntray
    eval(['pop0.l' num2str(i) '=(1/ninterval)*ones(nvar,ninterval);']);
    % initial probaility matrix
end

nsol0(1)=ceil(nsol/ntray);
for i=2:ntray-1;nsol0(i)=nsol0(1);end
nsol0(ntray)=nsol-(ntray-1)*nsol0(1);
%%%

iter=0;% iteration number

while iter < nloop
	iter=iter+1;
    %%%%%%%%% Create a population according to the probability matrix pop0 and nsol    
    x0=[];
    for i=1:ntray
        eval(['pop0li=pop0.l' num2str(i) ';']);
        xx0=pbilsv(pop0li,nsol0(i),a,b);
        x0=[x0 xx0];
    end
    %%%%%%%%% Convert binary strings to real parameters between a and b
    [x,f,g] = de_reproduct(fun,x0,ppareto,a,b,pc,F,CR);

    [ppareto,fpareto,gpareto,A]=pbil_selection(x,f,g,ppareto,fpareto,gpareto,A,narchive1,narchive2);

%     figure(2),clf,plot(fpareto(1,:),fpareto(2,:),'sb');hold on
%     plot(f(1,:),f(2,:),'o')
%     pause(1)
%%%%%%%%%%%%%%%%%%%%%%%%update the probability vector %%%%%%%%%
    npareto=size(ppareto,2);    
    if npareto < ntray
        for i=1:ntray
            Isl=ceil(rand*npareto);
            bestsol(:,i)=mean(ppareto(:,Isl),2);
        end
    else
        bestsol=fgroup01(fpareto,ppareto,ntray);
%         islbest=farchive22(fpareto,ntray);
%         bestsol=ppareto(:,islbest);
    end
%%%%%%%%%%%%%%%%%%%%%%%% update the probability vector %%%%%%%%%
    LR=LR0+0.1*rand*(-1)^round(rand);
    for i=1:ntray
        eval(['pop0li=pop0.l' num2str(i) ';']);
        pop0li=pop_update(pop0li,LR,bestsol(:,i),Mut_prob,Mut_shft,a,b);
        eval(['pop0.l' num2str(i) '=pop0li;']);
    end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     hv=HVcal(ppareto,fpareto,gpareto,ifun);
%     clc
%     disp(['iter = ' num2str(iter)]);
%     disp(['hv = ' num2str(hv)]);

    % Save Results
    rst.ppareto{iter}=ppareto;
    rst.fpareto{iter}=fpareto;
    rst.gpareto{iter}=gpareto;
    rst.timestamp=datetime('now');
end
% [ppareto,fpareto,gpareto,A]=pbil_selection0([],[],[],ppareto,fpareto,gpareto,narchive1);
% save(fout,'ppareto','fpareto','gpareto')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%sub programs%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x2,f2,g2] = de_reproduct(fun,x0,x1,a,b,pc,F,CR)
[m0,n0]=size(x0);[m1,n1]=size(x1);
if n1<n0
    ndiff=n0-n1;
    isl=randperm(n0);
    isl=isl(1:ndiff);
    x1=[x1 x0(:,isl)];
elseif n1>n0
    isl=randperm(n1);
    isl=isl(1:n0);
    x1=x1(:,isl);
end
n1=n0;
% DE point generation for n1 parents (x1)
% x0 = population from the previous generation
% x1 = ramdomly selected parents

for i=1:n1
    Pi=x1(:,i);% current parent
    nr=randperm(n0);
    i1=nr(1);i2=nr(2);
    Pi1=x0(:,i1);% Randomly seletced individual 1
    Pi2=x0(:,i2);% Randomly seletced individual 2
    Ci=Pi+unifrnd(0.5,0.8)*(Pi1-Pi2);
    for j=1:m1
        if Ci(j) < a(j)
            Ci(j)=a(j);
        elseif Ci(j) > b(j)
            Ci(j)=b(j);
        end
    end
    
    if rand < pc
        for j=1:m1 % binary crossover
            if rand < CR
                Vi(j,1)=Ci(j);
            else
                Vi(j,1)=Pi(j);
            end
        end
        Ci=Vi;
    end
    [f2(:,i),g2(:,i)]=feval(fun,Ci);
    x2(:,i)=Ci;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop1=pop_update(pop0,LR,bestsol,Mut_prob,Mut_shft,a,b)
pop1=pop0;
[nvar,nintv]=size(pop0);
for i=1:nvar
    besti=bestsol(i);
    k=0;sw=1;
    dxi=linspace(a(i),b(i),nintv+1);
    while sw==1
        k=k+1;
        if besti>=dxi(k)&besti<dxi(k+1)
            isl=k;
            sw=0;
        end
        if k==nintv-1
            sw=0;
            isl=nintv;
        end
    end
    pop1(i,isl)=(1-LR)*pop0(i,isl) + LR*(1);
	if rand < Mut_prob
        islm=ceil(rand*nintv);
	    pop1(i,islm)=(1-Mut_shft)*pop1(i,islm) + round(rand)*Mut_shft;
    end
    pop1(i,:)=pop1(i,:)/sum(pop1(i,:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x=pbilsv(pop0,nsol,a,b)

% Create solution vectors corresponding to pop0 (in PBIL)
% 
[m,n]=size(pop0);
for i=1:m;
    pop0i=pop0(i,:);
    n1=round(nsol*pop0i);
    if sum(n1) < nsol
        nxc=nsol-sum(n1);
        [n1max,nmax]=max(n1);
        n1(nmax)=n1(nmax)+nxc;
    end
    dx=linspace(a(i),b(i),n+1);
    x0j=[];
    for j=1:n
        for k=1:n1(j)
            x0jk=max(dx(j),min(dx(j+1),dx(j)+rand*(dx(j+1)-dx(j))));
            x0j=[x0j x0jk];
        end
    end
    npst=randperm(nsol);
    x(i,:)=x0j(1,npst(1:nsol));
end	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pareto1,fpareto1,gpareto1,A]=pbil_selection0(x1,f1,g1,pareto,fpareto,gpareto,narchive)
%
% GA Elite strategy
% keep 1 elite from the old generation and another from the
% new generation
x=[pareto x1];
f=[fpareto f1];
g=[gpareto g1];

[m0,n0]=size(fpareto);
% [m1,n1]=size(x);
n1=narchive;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update the external Pareto archive
function [pareto1,fpareto1,gpareto1,A]=pbil_selection(x1,f1,g1,pareto,fpareto,gpareto,A,narchive1,narchive2)
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

if nndm > narchive2
    nsl=farchive(fpareto1,narchive1);
    pareto1=pareto1(:,nsl);
    fpareto1=fpareto1(:,nsl);
    gpareto1=gpareto1(:,nsl);
    A=A(nsl,nsl);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seletion used when the number of non-dominated solutions
% exceeds the archive size
function iselect=farchive(f,narchive)
if size(f,1)==2
    iselect=farchive22(f,narchive);
else
    iselect=farchive72(f,narchive);
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

if mg1<=0&mg2<=0
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
    if icount11 == n & icount12 > 0
        p1=1;
    else
        p1=0;
    end
    if icount21 == n & icount22 > 0
        p2=1;
    else
        p2=0;
    end
elseif mg1 <=0 & mg2 > 0
    p1=1;p2=0;
elseif mg2 <=0 & mg1 > 0
    p1=0;p2=1;
else
    if mg1 <= mg2
        p1=1;p2=0;
    else
        p1=0;p2=1;
    end
end
function iselect=farchive22(f,narchive)

[m,n]=size(f);

fmax=max(f,[],2);fmin=min(f,[],2);
fdel=max(fmax-fmin,1e-5);
for i=1:n
    f(:,i)=(f(:,i)-fmin)./fdel;
end
nn=narchive;
x1=linspace(0,1,nn);
x2=1-x1;
iselect0=1:n;
for i=1:narchive
    df=[f(1,:)-x1(i);f(2,:)-x2(i)];
    d2=std(df);
    [d2min,imin]=min(d2);
    iselect(i)=iselect0(imin);
    iselect0(imin)=[];
    f(:,imin)=[];
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
function x=bin2real(bin,a,b);% convert binary to decimal and then to real numbers
    [m,n]=size(bin);
    nvar=length(a);
    nbit=m/nvar;

    for i=1:n
        for j=1:nvar
            x(j,i)=bin2dec(bin((j-1)*nbit+1:j*nbit,i),a(j),b(j));
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
    
function rG=fgroup01(f,x,narchive)
[m,n]=size(f);
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
for k=1:size(A,1)
    islk=A(k,1:A(k,n+2));
    rG(:,i)=mean(x(:,islk),2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% End of file %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%