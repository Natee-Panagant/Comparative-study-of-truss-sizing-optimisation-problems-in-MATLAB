function rst = UPSEMOA(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unrestricted Population Size Evolutionary Multiobjective      %
% Optimization Algorithm (UPS-EMOA)                             %
%                    Using real codes                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T. Attokoski and K. Miettinen, 2010, Efficient evolutionary 
% approach to approximating the Parto-optimal set in multiobjective 
% optimization, UPS-EMOA, Optimization Methods & Software, 25:841-858 
% VARIABLE DEFINITIONS
% a = lower limit
% b = upper limit
% nvar = the number of design variables
% nloop = the number of iteration
% nsol = the number of individuals(genes) in a population
% narchive = the maximum external Pareto archive size
% fun = m-file for evaluating objective functions and constraints
% fundat = mat-file containing initial population, objectives, and
% constraints
% fout = output mat-file containing non-dominated solutions
%


%%%%%%%%%%Initialisation Process %%%%%%%%%%%%
rand('state',sum(100*clock));% reset randomisation
tic

pc=0.7;% crossover probability
F=0.8;% Scaling factor for differentail evolution (DE) operator
CR=0.5;% probability of chosing element from offspring in crossover

minsize=10;% minimum population size
burstsize=nsol;%25;% no. of solution to be used as parents

% a=0*ones(nvar,1);% lower bounds
% b=1*ones(nvar,1);% upper bounds

% Initial population
[bin0,x0,f0,g0]=moea_initialisation(fun,nsol,nvar,nbit,a,b);
% This optimiser will not use bin0.
A0=fdominance(f0,g0);% dominance matrix for f0

%%%%%%%%% Start UPS-EMOA search %%%%%%%%
tic %
iter = 0;% iteration number
maxeval=nloop*nsol;% maximum number of function evaluations
neval=0;% total number of function evaluations so far

while neval < maxeval
	
	iter = iter+1;
    
    x1=upsemoa_select01(x0,burstsize);% select burstsize parents
	[x2,f2,g2]=upsemoa_reproduct(fun,x0,x1,a,b,pc,F,CR);%DE Point generation
    
    neval=neval+burstsize;% total number of function evaluations so far
    
    [x3,f3,g3,A3]=upsemoa_select02(x2,f2,g2,x0,f0,g0,A0,minsize);% new population
    
%     figure(1),clf,plot(f3(1,:),f3(2,:),'o'),pause(1)
    % the following operator is created just for the cases that the
    % population size (x3) might be too large (> 500 individuals).
    if size(f3,2) > 500
        Iselect=farchive(f3,500);
        x3=x3(:,Iselect);
        f3=f3(:,Iselect);
        g3=g3(:,Iselect);
        A3=A3(Iselect,Iselect);
    end
    x0=x3;f0=f3;g0=g3;A0=A3;
    
	% Save results
    [ppareto,fpareto,gpareto,~]=pareto_sorting(x0,f0,g0,[],[],[],[],narchive);
    rst.ppareto{iter}=ppareto;
    rst.fpareto{iter}=fpareto;
    rst.gpareto{iter}=gpareto;
    rst.timestamp=datetime('now');
end

% ppareto,fpareto,gpareto = Pareto archives for design variables,
% objectives, and constraints
% Time=toc;   %  in seconds 
% save(fout,'ppareto','fpareto','gpareto','Time')
% Save the final solutions to fout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A=fdominance(f,g);
[m,n]=size(f);
A=zeros(n,n);
for i=1:n
    fi=f(:,i);
    gi=g(:,i);
    A(i,i)=0;
    for j=(i+1):n
        fj=f(:,j);
        gj=g(:,j);
        %%%%%%%%%%%%%%%%%%%%%%%%%
        [p_count1,p_count2]=fdominated(fi,gi,fj,gj);
        A(i,j)=p_count1;
        A(j,i)=p_count2;
        %%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x1=upsemoa_select01(x0,burstsize)
% select burstsize parents randomly
[m,n]=size(x0);
if n < burstsize
    x1=x0;
    for i=1:(burstsize-n)
        x1=[x1 x0(:,ceil(rand*n))];
    end
elseif n== burstsize
    x1=x0;
else
    nsh=randperm(n);
    x1=x0(:,nsh(1:burstsize));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pareto1,fpareto1,gpareto1,A]=upsemoa_select02(x1,f1,g1,pareto,fpareto,gpareto,A,minsize)
%
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
if nndm < minsize
    [Bsort,nsort]=sort(B);
    Iselect=nsort(1:minsize);
    pareto1=x(:,Iselect);
    fpareto1=f(:,Iselect);
    gpareto1=g(:,Iselect);
    A=A(Iselect,Iselect);
else
    pareto1=x(:,Indm);
    fpareto1=f(:,Indm);
    gpareto1=g(:,Indm);
    A=A(Indm,Indm);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x2,f2,g2] = upsemoa_reproduct(fun,x0,x1,a,b,pc,F,CR)
[m0,n0]=size(x0);[m1,n1]=size(x1);

% DE point generation for n1 parents (x1)
% x0 = population from the previous generation
% x1 = ramdomly selected parents

for i=1:n1
    Pi=x1(:,i);% current parent
    nr=randperm(n0);
    i1=nr(1);i2=nr(2);
    Pi1=x0(:,i1);% Randomly seletced individual 1
    Pi2=x0(:,i2);% Randomly seletced individual 2
    Ci=Pi+F*(Pi1-Pi2);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Crowding distance algorithm                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iselect=farchive(f,narchive)

[m,n]=size(f);

fmax=max(f,[],2);[fmin nmin]=min(f,[],2);
fdel=max(fmax-fmin,1e-5);

for i=1:n
    f(:,i)=(f(:,i)-fmin)./fdel;
end

ii=1:n;
isl0=[];
for i=1:(n-narchive)
    fcrownd=fcrownding(f(:,ii));
    [fmin,nmin]=min(fcrownd);
    isl0=[isl0 ii(nmin)];
    ii(nmin)=[];
end
iselect=setdiff(1:n,isl0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fc=fcrownding(f);
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
%%%%%%%%%%%%%%%%%%%%%%%%End of file%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
