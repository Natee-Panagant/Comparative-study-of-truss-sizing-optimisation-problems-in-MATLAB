function rst = NSGA_II(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-Dominated Sorting Genetic Algorithm II                    %
%                   using real codes                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deb, K., Pratap, A., Agarwal, S., and Meyarivan, T. (2002). “A fast 
% and elitist multiobjective genetic algorithm: NSGAII.” IEEE Trans. 
% on Evolutionary Computation, 6(2),182-197.
%
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

%%%%%%%%%%Initialisation Process %%%%%%%%%%%%
rand('state',sum(100*clock));
tic

pm=0.1;% pm = mutation probability

% Initial population
[bin0,x0,f0,g0]=moea_initialisation(fun,nsol,nvar,nbit,a,b);
% This optimiser will not use bin0.

pop0=x0;
[pop0,f0,g0,A0]=nsga_select(pop0,f0,g0,[],[],[],[]);
%%%%%%%%% Start NSGA search %%%%%%%%
% tic %
iter = 0;
while iter < nloop
	
	iter = iter+1;

	pop1=nsga_crossover(pop0,a,b);%crossover
	[pop2,f2,g2]=nsga_mutate(fun,pop1,a,b,pm);%mutation
	
    [pop3,f3,g3,A3]=nsga_select(pop2,f2,g2,pop0,f0,g0,A0);%selection

	pop0=pop3;f0=f3;g0=g3;A0=A3;
    
    % Save results
    [ppareto,fpareto,gpareto,~]=pbil_selection0([],[],[],pop0,f0,g0,narchive);
    rst.ppareto{iter}=ppareto;
    rst.fpareto{iter}=fpareto;
    rst.gpareto{iter}=gpareto;
    rst.timestamp=datetime('now');
end

% [ppareto,fpareto,gpareto]=pareto_sorting(pop0,f0,g0,A0);

% Time=toc   %   second 
% save(fout,'ppareto','fpareto','gpareto')
% Save the final solutions to fout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pop3,f3,g3,A3]=nsga_select(pop2,f2,g2,pop1,f1,g1,A1);
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
for i=1:n1
    fi=f(:,i);
    gi=g(:,i);
    A(i,i)=0;
    for j=(n0+1):n1
        fj=f(:,j);
        gj=g(:,j);
        %%%%%%%%%%%%%%%%%%%%%%%%%
        [p_count1,p_count2]=fdominated(fi,gi,fj,gj);
        A(i,j)=p_count1;
        A(j,i)=p_count2;
    end
end
%%%%%%%%%%%%%%%%%%%%%%
for i=(n0+1):n1
    fi=f(:,i);
    gi=g(:,i);
    A(i,i)=0;
    for j=i+1:n1
        fj=f(:,j);
        gj=g(:,j);
        %%%%%%%%%%%%%%%%%%%%%%%%%
        [p_count1,p_count2]=fdominated(fi,gi,fj,gj);
        
        % p_count1 = 1 if i dominates j
        % p_count2 = 1 if j dominates i
        
        A(i,j)=p_count1;
        A(j,i)=p_count2;
        %%%%%%%%%%%%%%%%%%%%%%%%%        
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
elseif ncheck > n3&icheck==1
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B=ndlevel_sort(A);
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
function [ppareto,fpareto,gpareto]=pareto_sorting(pop,f,g,A);
B=sum(A,1)
col=[];
for i=1:size(pop,2)
    if B(i)==0
        col=[col i];
    end
end
ppareto=pop(:,col);
fpareto=f(:,col);
gpareto=g(:,col);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constrained optimisation domination sorting
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
    [p1,p2]=fdominated_uncon(g1,g2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unconstrained optimisation domination sorting
function [p1,p2]=fdominated_uncon(f1,f2)
n=length(f1);

icount11=0;
icount12=0;
icount21=0;
icount22=0;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x1=nsga_crossover(x0,a,b)
% real-code crossover in S. Srisomporn, and S. Bureerat, Geometrical design
% of plate-fin heat sinks using hybridization of MOEA and RSM, IEEE 
% Transaction on components and Packaging Technology, 31 (2008), pp. 351-360. 
[m,n]=size(x0);
x1=x0;
for i=1:n
   S=1+floor(0.999*rand*n);
   T=1+floor(0.999*rand*n);
   PR=rand;
   stl=-0.25+1.5*rand;
   for j=1:m
      if PR < 0.33
         if rand < 0.5
            x1(j,i)=x0(j,S);
         else
            x1(j,i)=x0(j,T);
         end
      elseif  PR >= 0.33 & PR < 0.66
         x1(j,i)=x0(j,S)+stl*(x0(j,T)-x0(j,S));
      else
         x1(j,i)=x0(j,S)+(-0.25+1.5*rand)*(x0(j,T)-x0(j,S));
      end
      
      if x1(j,i) < a(j)
          x1(j,i)=a(j);
      elseif x1(j,i) > b(j)
          x1(j,i)=b(j);
      end
      
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x1,f1,g1] = nsga_mutate(fun,x,a,b,pm)
% real-code mutation in S. Srisomporn, and S. Bureerat, Geometrical design
% of plate-fin heat sinks using hybridization of MOEA and RSM, IEEE 
% Transaction on components and Packaging Technology, 31 (2008), pp.
% 351-360. 
[m,n]=size(x);
alpha=0.5;
n0=ceil(pm*n);
%%%%%%%%%%%%%%%% main steps %%%%%%%%%%%%%%%%%%%%%%%%%%
x1=x;
for i=1:n0
    irow=max([ceil((rand)*m) 1]);
    icol=max([ceil((rand)*n) 1]);
    PR=rand;
	if PR < 0.5
		x1(irow,icol)=max([a(irow),a(irow)+(-0.25+1.25*rand)*(x1(irow,icol)-a(irow))]);
    else
		x1(irow,icol)=min([b(irow),x1(irow,icol)+(1.25*rand)*(b(irow)-x1(irow,icol))]);
    end
end
x2=min(1,max(0,x1));
for i=1:n
    [f1(:,i),g1(:,i)]=feval(fun,x2(:,i));
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
        [p_count1,p_count2]=fdominated2(fi,gi,fj,gj);
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

function [p1,p2]=fdominated2(f1,g1,f2,g2)
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