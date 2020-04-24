function rst = MOEA_D(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)
%%%%%%%%%%%%%%%%% MAIN PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Initialisation Process %%%%%%%%%%%%
rand('state',sum(100*clock));
tic
% for minimisation
[bin0,x0,f0,g0]=moea_initialisation(fun,nsol,nvar,nbit,a,b);

% wieghting factors
nobj=size(f0,1);
lamb0=lhsamp(nsol,nobj,zeros(nobj,1),ones(nobj,1));

for i=1:nsol
    lamb(:,i)=lamb0(:,i)./sum(lamb0(:,i));
end

T=6;% no. of neighboring weight vector

% initial Pareto archive
[ppareto,fpareto,gpareto,A]=pareto_selection(x0,f0,g0,[],[],[],[],narchive);
fp0=fpenal(f0,g0);
[x0,f0,g0,fp0]=findbestsol(x0,f0,g0,fp0,lamb);
%z=min(fp0,[],2);%best values of fpenal
for i=1:nsol
    li=lamb(:,i);
    Dl(i,i)=0;
    for j=(i+1):nsol
        lj=lamb(:,j);
        Dl(i,j)=sum((li-lj).*(li-lj));
        Dl(j,i)=Dl(i,j);
    end
    [sDl,ns]=sort(Dl(i,:));
    B(i,1:T)=ns(1:T);% neighbours of a weighting vector
end

iter=0;% iteration number
pm=0.1;% mutation rate
while iter < nloop
	iter=iter+1;
    
    % reproduce offspring
    for i=1:nsol
        nsh=randperm(T);
        x1(:,i)=moead_crossover(x0(:,B(i,nsh(1))),x0(:,B(i,nsh(2))),a,b);
        x1(:,i)=moead_mutate(x1(:,i),a,b,pm);

        [f1(:,i),g1(:,i)]=feval(fun,x1(:,i));
        fp1(:,i)=fpenal(f1(:,i),g1(:,i));
    end
    
    [x0,f0,g0,fp0]=findbestsol([x0 x1],[f0 f1],[g0 g1],[fp0 fp1],lamb);
    
    [ppareto,fpareto,gpareto,A]=pareto_selection(x1,f1,g1,ppareto,fpareto,gpareto,A,narchive);
    
    % Save results
    rst.ppareto{iter}=ppareto;
    rst.fpareto{iter}=fpareto;
    rst.gpareto{iter}=gpareto;
    rst.timestamp=datetime('now');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%sub programs%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x1,f1,g1,fp1]=findbestsol(x0,f0,g0,fp0,lamb)
[m,n]=size(lamb);
for i=1:n
    fpi=lamb(:,i)'*fp0;
    [fpimin,nmin]=min(fpi);
    x1(:,i)=x0(:,nmin);
    f1(:,i)=f0(:,nmin);
    g1(:,i)=g0(:,nmin);
    fp1(:,i)=fp0(:,nmin);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fp=fpenal(f,g)
% Exterior penalty function technique
c=100;
for i=1:size(f,2);
    fp(:,i)=f(:,i)+c*sum(max(0,g(:,i)).^2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R=ParetoRank(x,f,g)
[m,n]=size(x);

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
R=sum(A,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x3=moead_crossover(x1,x2,a,b)
m=length(x1);
PR=rand;
stl=-0.25+1.5*rand;
for j=1:m
   if PR < 0.33
      if rand < 0.5
         x3(j,1)=x1(j,1);
      else
         x3(j,1)=x2(j,1);
      end
   elseif  PR >= 0.33 & PR < 0.66
      x3(j,1)=x1(j,1)+stl*(x2(j,1)-x1(j,1));
   else
      x3(j,1)=x1(j,1)+(-0.25+1.5*rand)*(x2(j,1)-x1(j,1));
   end
      
   if x3(j,1) < a(j)
       x3(j,1)=a(j);
   elseif x3(j,1) > b(j)
       x3(j,1)=b(j);
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x1] = moead_mutate(x,a,b,pm)
m=length(x);
x1=x;
%%%%%%%%%%%%%%%% main steps %%%%%%%%%%%%%%%%%%%%%%%%%%
Pr=rand;
if Pr < pm
    irow=max([ceil((rand)*m) 1]);
    PR=rand;
    if PR < 0.5
        x1(irow,1)=max([a(irow),a(irow)+(-0.25+1.25*rand)*(x(irow,1)-a(irow))]);
    else
        x1(irow,1)=min([b(irow),x(irow,1)+(1.25*rand)*(b(irow)-x(irow,1))]);
    end
    x1=min(1,max(0,x1));
else
    x1=x;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pareto1,fpareto1,gpareto1,A]=pareto_selection(x1,f1,g1,pareto,fpareto,gpareto,A,narchive)
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
% Seletion used when the number of non-dominated solutions
% exceeds the archive size
function iselect=farchive(f,narchive)
if size(f,1)==2
    iselect=farchive22(f,narchive);
else
    iselect=farchive53(f,narchive);
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
%
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

function S = lhsamp(m, n, a, b)
%LHSAMP  Latin hypercube distributed random numbers
%
% Call:    S = lhsamp
%          S = lhsamp(m)
%          S = lhsamp(m, n)
%
% m : number of sample points to generate, if unspecified m = 1
% n : number of dimensions, if unspecified n = m
%
% a : lower bounds
% b : upper bounds
%
% S : the generated n dimensional m sample points chosen from
%     uniform distributions on m subdivions of the interval ({a}, {b})

% hbn@imm.dtu.dk  
% Last update April 12, 2002

if nargin < 1, m = 1; end
if nargin < 2, n = m; end
if nargin > 2 & nargin < 3, a = 0*ones(n,1),b = 0*ones(n,1); end

S = zeros(n,m);
for i = 1 : n
  S0 = (rand(1, m) + (randperm(m) - 1)) / m;
  S(i, :) = a(i)*ones(1,m)+(b(i)-a(i))*S0;
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

function iselect=farchive53(f,narchive)

[m,n]=size(f);

fmax=max(f,[],2);[fmin nmin]=min(f,[],2);
fdel=max(fmax-fmin,1e-5);

for i=1:n
    f(:,i)=(f(:,i)-fmin)./fdel;
end

ii=1:n;
isl0=[];
for i=1:(n-narchive)
    fcrownd=fcrownding2(f(:,ii));
    [fmin,nmin]=min(fcrownd);
    isl0=[isl0 ii(nmin)];
    ii(nmin)=[];
end
iselect=setdiff(1:n,isl0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fc=fcrownding2(f);
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
% figure(2),clf,plot(f(1,:),f(2,:),'o'),hold on
% for i=1:n
%     text(f(1,i),f(2,i),num2str(fc(i)))
% end
% figure(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% End of file %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%