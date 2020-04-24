function rst = DEMO(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)
%%%%%%%%%%%%%%%%% MAIN PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Initialisation Process %%%%%%%%%%%%
tic
rand('state',sum(100*clock));

pc=0.7;% crossover probability
F=0.8;% Scaling factor for differentail evolution (DE) operator
CR=0.5;% probability of chosing element from offspring in crossover

% for minimisation
% initialisation
[bin0,x0,f0,g0]=moea_initialisation(fun,nsol,nvar,nbit,a,b);

[x0,f0,g0,A0]=pareto_select(x0,f0,g0,[],[],[],[]);

ppareto=[];
fpareto=[];
gpareto=[];

iter=0;% iteration number
while iter < nloop
	iter=iter+1;
    
    % reproduce offspring
    [x1,f1,g1] = de_reproduct(fun,x0,a,b,pc,F,CR);
    [x2,f2,g2,A2]=pareto_select(x1,f1,g1,x0,f0,g0,A0);
%     fpareto=f2
    
%     figure(1),clf,hold on,
%     plot(fpareto(1,:),fpareto(2,:),'s')
%     plot(f0(1,:),f0(2,:),'or')
%     pause(1)
    
    x0=x2;f0=f2;g0=g2;A0=A2;
    
    [ppareto,fpareto,gpareto]=pareto_sorting([ppareto,x0,],[fpareto,f0],[gpareto,g0],narchive);

    
    % Save Results
    rst.ppareto{iter}=ppareto;
    rst.fpareto{iter}=fpareto;
    rst.gpareto{iter}=gpareto;
    rst.timestamp=datetime('now');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%sub programs%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x1,f1,g1] = de_reproduct(fun,x0,a,b,pc,F,CR)
[m,n]=size(x0);
for i=1:n
    Pi=x0(:,i);% current parent
    I=[1:(i-1) (i+1):n];
    nr=randperm(n-1);
    i1=I(nr(1));i2=I(nr(2));
    Pi1=x0(:,i1);% Randomly seletced individual 1
    Pi2=x0(:,i2);% Randomly seletced individual 2
    Ci=Pi+unifrnd(0.5,0.8)*(Pi1-Pi2);
    for j=1:m
        if Ci(j) < a(j)
            Ci(j)=a(j);
        elseif Ci(j) > b(j)
            Ci(j)=b(j);
        end
    end
    
    if rand < pc
        for j=1:m % binary crossover
            if rand < CR
                Vi(j,1)=Ci(j);
            else
                Vi(j,1)=Pi(j);
            end
        end
        Ci=Vi;
    end
    [f1(:,i),g1(:,i)]=feval(fun,Ci);
    x1(:,i)=Ci;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pop3,f3,g3,A3]=pareto_select(pop2,f2,g2,pop1,f1,g1,A1)
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
% [B
% sum(A,1)]
% pause
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
% Seletion used when the number of non-dominated solutions
% exceeds the archive size
function iselect=farchive(f,narchive)
if size(f,1)==2
    iselect=farchive52(f,narchive);
else
    iselect=farchive53(f,narchive);
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
function [ppareto,fpareto,gpareto]=pareto_sorting(pop,f,g,narchive)
[~,uind]=unique(pop','rows');
pop=pop(:,uind);
f=f(:,uind);
g=g(:,uind);

n1=size(f,2);
for i=1:n1
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
B=sum(A,1);
col=[];
for i=1:size(pop,2)
    if B(i)==0
        col=[col i];
    end
end
ppareto=pop(:,col);
fpareto=f(:,col);
gpareto=g(:,col);

if numel(col) > narchive
    nsl=farchive(fpareto,narchive);
    ppareto=ppareto(:,nsl);
    fpareto=fpareto(:,nsl);
    gpareto=gpareto(:,nsl);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fc=fcrownding(f);
[m,n]=size(f);
delf=max(f,[],2)-min(f,[],2);

[fs,nsort]=sort(f(1,:));
ff=f(:,nsort);
fc0=zeros(1,n);fc0(1)=1e10;fc0(n)=1e10;

for i=2:(n-1)
    for j=1:m
        fc0(i)=fc0(i)+abs((ff(j,i+1)-ff(j,i-1))/delf(j));
    end
end     
fcmax=max([1 (fc0(2:(n-1)))]);
fc0(1)=fcmax+10;fc0(n)=fcmax+10;
fc(nsort)=fc0;
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
function iselect=farchive52(f,narchive)

[m,n]=size(f);

fmax=max(f,[],2);[fmin nmin]=min(f,[],2);
fdel=max(fmax-fmin,1e-5);

for i=1:n
    f(:,i)=(f(:,i)-fmin)./fdel;
end

ii=1:n;
isl0=[];
for i=1:(n-narchive)
    fcrownd=fcrownding52(f(:,ii));
    [fmin,nmin]=min(fcrownd);
    isl0=[isl0 ii(nmin)];
    ii(nmin)=[];
end
iselect=setdiff(1:n,isl0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fc=fcrownding52(f)
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
    fcrownd=fcrownding53(f(:,ii));
    [fmin,nmin]=min(fcrownd);
    isl0=[isl0 ii(nmin)];
    ii(nmin)=[];
end
iselect=setdiff(1:n,isl0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fc=fcrownding53(f)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% End of file %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%