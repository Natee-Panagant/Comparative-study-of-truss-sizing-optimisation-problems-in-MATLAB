function [f,g,rst] = f25barSI(x,varargin)
% 25-bar 3D truss optimization problem
% x = design variables (input - 0<=x<=1)
    
    % Decoding design variables
    LB = zeros(8,1); % lower-bound of design variables
    UB = 41*ones(8,1); % upper-bound of design variables
    x = LB+(UB-LB).*x; % truncate x to lower and upper-bound
    x = max(1,ceil(x)); % round x as cross-section size index
    
    % Material properties
    rho=7850; % Density (kg/m^3)
    sigma_a =400e6; %allowable stress(Pa)
    E =200e9; %modulus of elasticity(Pa)
    
    % Geometry (position in m)
    node=[  -1      0       5;
            1       0       5;
            -1      1       2.5;
            1       1       2.5;
            1       -1      2.5;
            -1      -1      2.5;
            -2.5	2.5     0;
            2.5     2.5     0;
            2.5     -2.5    0;
            -2.5    -2.5    0];
    
    ele=[1 2%element1
   		 1  4%element2
         2  3%element3
         1  5%element4
         2  6%element5
         2  4
         2  5
         1  3
         1  6
         3  6
         4  5
         3  4
         5  6
         3  10
         6  7
         4  9
         5  8
         4  7
         3  8
         5  10
         6  9
         6  10
         3  7
         4  8
         5  9];%element25
   
    % Cross-section areas (m^2)    
    section = (1e-3)*(1:.5:21); % List of discrete cross-section area
    A(1)=section(x(1));
    A(2:5,1)=section(x(2));
    A(6:9,1)=section(x(3));
    A(10:11,1)=section(x(4));
    A(12:13,1)=section(x(5));
    A(14:17,1)=section(x(6));
    A(18:21,1)=section(x(7));
    A(22:25,1)=section(x(8));
    
    % Loadings
	%          node    dof      load(N)      %dof = 1 for x-axis, 2 for y-axis, 3 for z-axis
    Load{1} = [ 1       1       1e5;
                1       2       -1e6
                1       3       -1e6
                2       2       -1e6
                2       3       -1e6
                3       1       5e4
                6       1       6e4];
            
    % Prescribed displacement
    %      node    dof  displacement(m)      %dof = 1 for x-axis, 2 for y-axis, 3 for z-axis
    BC = [	7       1       0;
            7       2       0;
            7       3       0;
            8       1       0;
            8       2       0;
            8       3       0;
            9       1       0;
            9       2       0;
            9       3       0;
            10      1       0;
            10      2       0;
            
            10      3       0];

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