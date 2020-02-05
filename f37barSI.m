function [f,g,rst] = f37barSI(x,varargin)
% 37-bar 2D truss optimization problem
% x = design variables (input - 0<=x<=1)
    
    % Decoding design variables
    LB = zeros(15,1); % lower-bound of design variables
    UB = 41*ones(15,1); % upper-bound of design variables
    x = LB+(UB-LB).*x; % truncate x to lower and upper-bound
    x = max(1,ceil(x)); % round x as cross-section size index
    
    % Material properties
    rho=7850; % Density (kg/m^3)
    sigma_a =400e6; %allowable stress(Pa)
    E =200e9; %modulus of elasticity(Pa)
    
    % Geometry (position in m)
    node(1:11,:)=[(0:10)',zeros(11,1)];
    node(12:20,:)=[(1:9)',ones(9,1)];
    
    ele=[1 2           %element1
         2  3           %element2
         3  4           %element3
         4  5           
         5  6           
         6  7           
         7  8           
         8  9           
         9  10           
         10 11           
         12 13           
         13 14           
         14 15           
         15 16           
         16 17           
         17 18           
         18 19           
         19 20           
         1  12           
         2  12           
         3  12           
         3  13           
         4  13           
         4  14           
         5  14           
         5  15           
         6  15           
         6  16           
         6  17           
         7  17           
         7  18           
         8  18           
         8  19          
         9  19           
         9  20          
         10 20           
         11 20           ];%element37
       
    % Cross-section areas (m^2)
    section = (1e-3)*(1:.5:21); % List of discrete cross-section area
    A=zeros(37,1);
    A(1:10,:)=section(x(1))*ones(10,1);
    A([19 37],:)=section(x(2))*ones(2,1);
    A([20 36],:)=section(x(3))*ones(2,1);
    A([21 35],:)=section(x(4))*ones(2,1);
    A([11 18],:)=section(x(5))*ones(2,1);
    A([22 34],:)=section(x(6))*ones(2,1);
    A([23 33],:)=section(x(7))*ones(2,1);
    A([12 17],:)=section(x(8))*ones(2,1);
    A([24 32],:)=section(x(9))*ones(2,1);
    A([25 31],:)=section(x(10))*ones(2,1);
    A([13 16],:)=section(x(11))*ones(2,1);
    A([26 30],:)=section(x(12))*ones(2,1);
    A([27 29],:)=section(x(13))*ones(2,1);
    A([14 15],:)=section(x(14))*ones(2,1);
    A(28,:)=section(x(15))*ones(1,1);
    
    % Loadings
	%          node    dof      load(N)      %dof = 1 for x-axis, 2 for y-axis, 3 for z-axis
    Load{1} = [(2:10)', 2*ones(9,1), -10e4*ones(9,1)];
    
    % Prescribed displacement
    %      node    dof  displacement(m)      %dof = 1 for x-axis, 2 for y-axis, 3 for z-axis
    BC=[	1       1       0;
            1       2       0;
            11      2       0];
        
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