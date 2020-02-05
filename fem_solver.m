function [f,g,results] = fem_solver(node,ele,A,E,Load,BC,rho,sigma_a,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D/3D Truss Finite Elment Method solver %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% node              = node list (size node = 2 for 2D problem, = 3 for 3D problem
% ele               = element list
% A                 = Cross-section area of member (size A = Nele*1)
% E                 = modulus of elasticity
% Load              = External Load ([node_index Fx Fy . . .])
% BC                = Prescribed displacement (Boundary conditions) [node_index dof_index displacement;...])
% Kaa               = Partitioned Stiffness Matrix
% partitioned_dof   = List of degrees of freedom(DOF) excluded prescribed DOF
% prescribed_dof    = List of prescribed DOF
    
    %%%%%%%%%%%%%%%%%%
    % Initialization %
    %%%%%%%%%%%%%%%%%%
    
    % Get Number of Problem Dimensions (ND)
    if size(node,2)==2
        ND=2; 
    elseif size(node,2)==3
        ND=3;
    else
        error('Invalid node data');
    end
    % Transform A to column vector
    if isrow(A)
        A=A';
    end
    Nnode=size(node,1);
    Ndof=ND*Nnode;
    Nele=size(ele,1);
    L=zeros(Nele,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Stiffness Assembly and Partition operator %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    NK=(2*ND)^2; % number of members of stiffness matrix
    iKsparse=zeros(1,NK*Nele);
    jKsparse=zeros(1,NK*Nele);
    kKsparse=zeros(1,NK*Nele);
    iSsparse=zeros(1,2*ND*Nele);
    jSsparse=zeros(1,2*ND*Nele);
    kSsparse=zeros(1,2*ND*Nele);
    for i=1:Nele
        n1=ele(i,1);
        n2=ele(i,2);
        Ai=A(i,1);
        coor=node([n1,n2],:);
        diff=coor(2,:)-coor(1,:);
        L(i)=norm(coor(2,:)-coor(1,:));
        diff_L=diff/L(i);
        % Compute Transformation matrix
        T=[diff_L       zeros(1,ND);
           zeros(1,ND)	diff_L];
       % Compute B matrix
        B=[-1/L(i) 1/L(i)];
        
        % Compute list of degree of freedom of i-th element
        dof=[ND*n1-(ND-1):ND*n1,ND*n2-(ND-1):ND*n2];
        % Compute Stiffness matrix
        ki=Ai*(T'*B'*E*B*T)*L(i); % K = A*integrate(T'B'EBT,0->L,dx)
        iKsparse(1,NK*i-NK+1:NK*i)=repmat(dof,[1,2*ND]);
        jKsparse(1,NK*i-NK+1:NK*i)=repelem(dof,2*ND);
        kKsparse(1,NK*i-NK+1:NK*i)=reshape(ki,1,[]);
        % Compute S matrix for stress calculation
        iSsparse(1,2*ND*i-2*ND+1:2*ND*i)=i*ones(1,2*ND);
        jSsparse(1,2*ND*i-2*ND+1:2*ND*i)=[ND*n1-ND+1:ND*n1,ND*n2-ND+1:ND*n2];
        kSsparse(1,2*ND*i-2*ND+1:2*ND*i)=E*B*T; % S=EBT, S(Nele x Ndof)*Qele(Ndof x 1) = stress(Nele x 1)
    end
    % Assembly Stiffness (K) and S matrix
    K=sparse(iKsparse,jKsparse,kKsparse,Ndof,Ndof);
    S=sparse(iSsparse,jSsparse,kSsparse,Nele,Ndof);
    % Matrix Partition - partitioned index generation
    partitioned_dof=true(Ndof,1);
    prescribed_dof=ND*(BC(:,1)-1)+BC(:,2);
    partitioned_dof(prescribed_dof)=false;
    prescribed_dof=~partitioned_dof;
    F=cell(numel(Load),1);
    for i=1:numel(Load)
        F{i}=zeros(Ndof,1);
        Find=ND*(Load{i}(:,1)-1)+Load{i}(:,2);
        Fval=Load{i}(:,3);
        F{i}(Find,1)=Fval;
    end
    
        % Solve
    NLC = numel(Load); % Number of Load Cases
    NBC = size(BC,1); % Number of Prescribed displacement (Boundary Conditions)
    [Nnode,ND] = size(node);
    displm_p = zeros(ND*Nnode-NBC,NLC);
    displm = zeros(ND*Nnode,NLC);
    compliance = zeros(1,NLC);
    stress = zeros(size(ele,1),NLC);
    for i=1:NLC    % solve each Load Cases
        displm_p(:,i) = K(partitioned_dof,partitioned_dof)\F{i}(partitioned_dof,1); %Partition dof
        displm(prescribed_dof,i) = BC(:,3);
        displm(partitioned_dof,i) = displm_p(:,i);
        compliance(i) = displm(:,i)'*F{i};
        stress(:,i) = S*displm(:,i);
    end
    max_compliance=max(compliance);
    %
    %
    if isnan(max_compliance)
       0; 
    end
    %
    %
    max_stress=max(abs(stress),[],2);
    % Obtain objective and constrained functions
    mass = rho*sum(A.*L);
    f = [mass;max_compliance];
    g = max(max_stress-sigma_a,0);
    
    % Return FEM results
    results=[];
    if numel(varargin)>0
        if strcmp(varargin,'post') % Return FEM data
            results=struct();
            results.mass=mass;
            results.displm=displm;
            results.stress=stress;
            results.compliance=compliance;
            results.partitioned_dof=partitioned_dof;
            results.prescribed_dof=prescribed_dof;
        end
    end
    %
    %
    %
%     s_a=max(max_stress)/sigma_a
%     a_s=sigma_a/max(max_stress)
    %
    %
    %
end