function [f,g,rst] = f120barSI(x,varargin)
% 120-bar 3D truss optimization problem
% x = design variables (input - 0<=x<=1)
    
    % Decoding design variables
    LB = zeros(7,1); % lower-bound of design variables
    UB = 41*ones(7,1); % upper-bound of design variables
    x = LB+(UB-LB).*x; % truncate x to lower and upper-bound
    x = max(1,ceil(x)); % round x as cross-section size index
    
    % Material properties
    rho=7850; % Density (kg/m^3)
    sigma_a =400e6; %allowable stress(Pa)
    E =200e9; %modulus of elasticity(Pa)
    
    % Geometry (position in m)
    node=zeros(49,3);
    for i=1:49
        %%%%%%%% layer z=0
        if i<=12
            node(i,:)=[15.89*cos((i-1)*pi/6)  15.89*sin((i-1)*pi/6)    0];
        elseif i>12&&i<=36
            node(i,:)=[12.04*cos((i-13)*pi/12)  12.04*sin((i-13)*pi/12)    3];
        elseif i>36&&i<=48
            node(i,:)=[6.94*cos((i-37)*pi/6)  6.94*sin((i-37)*pi/6)    5.85];
        elseif i==49
            node(i,:)=[0  0  7];
        else
        end
    end
    ele=zeros(120,2);
    nele=0;
    for i=1:12 %%%% group 1
        nele=nele+1;
        ele(nele,:) =[49 49-i];
    end
    for i=1:11 %%%% group 2
        nele=nele+1;
        ele(nele,:) =[36+i 36+i+1];
    end
        nele=nele+1;
        ele(nele,:) =[48 37];
    for i=1:12  %%%% group 3
        nele=nele+1;
        ele(nele,:) =[36+i 2*i+11];
    end
    for i=1:11  %%%% group 4
        nele=nele+1;
        ele(nele,:) =[37+i 2*i+12];
        nele=nele+1;
        ele(nele,:) =[37+i 2*i+14];
    end
        nele=nele+1;
        ele(nele,:) =[37 36];
        nele=nele+1;
        ele(nele,:) =[37 14];
    for i=1:23 %%%% group 5
        nele=nele+1;
        ele(nele,:) =[12+i 12+i+1];
    end  
        nele=nele+1;
        ele(nele,:) =[36 13];
    for i=1:12  %%%% group 6
        nele=nele+1;
        ele(nele,:) =[11+2*i i];
    end  
    for i=2:12  %%%% group 7
        nele=nele+1;
        ele(nele,:) =[i 2*i+10];
        nele=nele+1;
        ele(nele,:) =[i 2*i+12];
    end   
    nele=nele+1;
    ele(nele,:) =[1 14];
    nele=nele+1;
    ele(nele,:) =[1 36];
    
    % Cross-section areas (m^2)
    section = (1e-3)*(1:.5:21); % List of discrete cross-section area
    A = zeros(120,1);
    A(1:12,:)=section(x(1))*ones(12,1);
    A(13:24,:)=section(x(2))*ones(12,1);
    A(25:36,:)=section(x(3))*ones(12,1);
    A(37:60,:)=section(x(4))*ones(24,1);
    A(61:84,:)=section(x(5))*ones(24,1);
    A(85:96,:)=section(x(6))*ones(12,1);
    A(97:120,:)=section(x(7))*ones(24,1);
    
    % Loadings
	%          node    dof      load(N)      %dof = 1 for x-axis, 2 for y-axis, 3 for z-axis
    Load{1} = [(13:49)'     3*ones(37,1)    [(-5e5)*ones(24,1);(-15e5)*ones(12,1);-3e6]];
    
    % Prescribed displacement
    %      node    dof  displacement(m)      %dof = 1 for x-axis, 2 for y-axis, 3 for z-axis
    BC=[repelem((1:12)',3,1) repmat((1:3)',12,1) zeros(36,1)];
    
    % Solving with Finit Element Method (FEM)
    rst=[];
    if numel(varargin)>0
        if strcmp(varargin,'post') % Return FEM data
            [f,g,results] = fem_solver(node,ele,A,E,Load,BC,rho,sigma_a,'post');
            rst=struct();
            rst.node=node;
            rst.ele=ele;
            rst.A=A;
            rst.E=E;
            rst.Load=Load;
            rst.BC=BC;
            rst.rho=rho;
            rst.sigma_a=sigma_a;
            
            rst.mass=results.mass;
            rst.displm=results.displm;
            rst.stress=results.stress;
            rst.compliance=results.compliance;
            rst.partitioned_dof=results.partitioned_dof;
            rst.prescribed_dof=results.prescribed_dof;
        end
    else
        [f,g] = fem_solver(node,ele,A,E,Load,BC,rho,sigma_a);
    end
end