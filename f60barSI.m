function [f,g,rst] = f60barSI(x,varargin)
% 60-bar 2D truss optimization problem
% x = design variables (input - 0<=x<=1)
    
    % Decoding design variables
    LB = zeros(25,1); % lower-bound of design variables
    UB = 41*ones(25,1); % upper-bound of design variables
    x = LB+(UB-LB).*x; % truncate x to lower and upper-bound
    x = max(1,ceil(x)); % round x as cross-section size index
    
    % Material properties
    rho=7850; % Density (kg/m^3)
    sigma_a =400e6; %allowable stress(Pa)
    E =200e9; %modulus of elasticity(Pa)
    
    % Geometry (position in m)
    node=[     2.500000000000000                   0
               2.165063500000000   1.250000000000000
               1.250000000000000   2.165063500000000
                               0   2.500000000000000
              -1.250000000000000   2.165063500000000
              -2.165063500000000   1.250000000000000
              -2.500000000000000                   0
              -2.165063500000000  -1.250000000000000
              -1.250000000000000  -2.165063500000000
                               0  -2.500000000000000
               1.250000000000000  -2.165063500000000
               2.165063500000000  -1.250000000000000
               2.250000000000000                   0
               1.948557250000000   1.125000000000000
               1.125000000000000   1.948557250000000
                               0   2.250000000000000
              -1.125000000000000   1.948557250000000
              -1.948557250000000   1.125000000000000
              -2.250000000000000                   0
              -1.948557250000000  -1.125000000000000
              -1.125000000000000  -1.948557250000000
                               0  -2.250000000000000
               1.125000000000000  -1.948557250000000
               1.948557250000000  -1.125000000000000];
    ele=[1  2;
         2  3;
         3  4;
         4  5;
         5  6;
         6  7;
         7  8;
         8  9;
         9  10;
         10 11;
         11 12;
         12 1;
         13 14;
         14 15;
         15 16;
         16 17;
         17 18;
         18 19;
         19 20;
         20 21;
         21 22;
         22 23;
         23 24;
         24 13;
         1  14;
         2  15;
         3  16;
         4  17;
         5  18;
         6  19;
         7  20;
         8  21;
         9  22;
         10 23;
         11 24;
         12 13;
         13 2;
         14 3;
         15 4;
         16 5;
         17 6;
         18 7;
         19 8;
         20 9;
         21 10;
         22 11;
         23 12;
         24 1;
         13 1;
         14 2;
         15 3;
         16 4;
         17 5;
         18 6;
         19 7;
         20 8;
         21 9;
         22 10;
         23 11;
         24 12]; % Element conectivity 
   
    % Cross-section areas (m^2)
    section = (1e-3)*(1:.5:21); % List of discrete cross-section area
    A = zeros(60,1);
    A(49:60,1) = section(x(1)); 
    A([1 13],1) = section(x(2)); % Assign Group Areas as per symmetry 
    A([2 14],1) = section(x(3));
    A([3 15],1) = section(x(4));
    A([4 16],1) = section(x(5));
    A([5 17],1) = section(x(6)); 
    A([6 18],1) = section(x(7)); 
    A([7 19],1) = section(x(8));
    A([8 20],1) = section(x(9));
    A([9 21],1) = section(x(10));
    A([10 22],1) = section(x(11)); 
    A([11 23],1) = section(x(12)); 
    A([12 24],1) = section(x(13));
    A([25 37],1) = section(x(14));
    A([26 38],1) = section(x(15));
    A([27 39],1) = section(x(16)); 
    A([28 40],1) = section(x(17)); 
    A([29 41],1) = section(x(18));
    A([30 42],1) = section(x(19));
    A([31 43],1) = section(x(20));
    A([32 44],1) = section(x(21)); 
    A([33 45],1) =section( x(22)); 
    A([34 46],1) = section(x(23));
    A([35 47],1) = section(x(24));
    A([36 48],1) = section(x(25));
    
    % Loadings
	%          node    dof      load(N)      %dof = 1 for x-axis, 2 for y-axis, 3 for z-axis
    Load{1} = [ 1       1       -10e5;
                7       1       9e5];
    Load{2} = [ 15      1       -8e5;
                15      2       3e5;
                18      1       -8e5;
                18      2       3e5];
	Load{3} = [ 22      1       -20e5;
                22      2       10e5];
            
    % Prescribed displacement
    %      node    dof  displacement(m)      %dof = 1 for x-axis, 2 for y-axis, 3 for z-axis
    BC=[    10      1       0;
            10      2       0;
            16      1       0];
        
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