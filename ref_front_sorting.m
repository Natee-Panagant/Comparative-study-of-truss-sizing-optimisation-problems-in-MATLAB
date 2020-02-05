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

function iselect=farchive(f,narchive)
    if size(f,1)==2
        iselect=farchive52(f,narchive);
    elseif size(f,1)==3
        iselect=farchive53(f,narchive);
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
        fcrownd=fcrownding(f(:,ii));
        [fmin,nmin]=min(fcrownd);
        isl0=[isl0 ii(nmin)];
        ii(nmin)=[];
    end
    iselect=setdiff(1:n,isl0);
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
        fcrownd=fcrownding(f(:,ii));
        [fmin,nmin]=min(fcrownd);
        isl0=[isl0 ii(nmin)];
        ii(nmin)=[];
    end
    iselect=setdiff(1:n,isl0);
end
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