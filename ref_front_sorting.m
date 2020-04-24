function [pareto1,fpareto1,gpareto1,A]=ref_front_sorting(x1,f1,g1,pareto,fpareto,gpareto,A,narchive)
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