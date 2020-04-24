function [f,g,rst] = f72barSI(x,varargin)
% 72-bar 3D truss optimization problem
% x = design variables (input - 0<=x<=1)
    
    % Decoding design variables
    LB = zeros(16,1); % lower-bound of design variables
    UB = 41*ones(16,1); % upper-bound of design variables
    x = LB+(UB-LB).*x; % truncate x to lower and upper-bound
    x = max(1,ceil(x)); % round x as cross-section size index
    
    % Material properties
    rho=7850; % Density (kg/m^3)
    sigma_a =400e6; %allowable stress(Pa)
    E =200e9; %modulus of elasticity(Pa)
    
    % Geometry (position in m)
    node = [0       0       0;
            3       0       0;
            3       3       0;
            0       3       0;
            0       0       1.5;
            3       0       1.5;
            3       3       1.5;
            0       3       1.5;
            0       0       3;
            3       0       3;
            3       3       3;
            0       3       3;
            0       0       4.5;
            3       0       4.5;
            3       3       4.5;
            0       3       4.5;
            0       0       6;
            3       0       6;
            3       3       6;
            0       3       6];
    ele = [ 1   5;
            2   6;
            3   7;
            4   8;
            1   6;
            2   5;
            2   7;
            3   6;
            3   8;
            4   7;
            4   5;
            1   8;
            5   6;
            6   7;
            7   8;
            8   5;
            5   7;
            6   8;
            5   9;
            6   10;
            7   11;
            8   12;
            5   10;
            6   9;
            6   11;
            7   10;
            7   12;
            8   11;
            8   9;
            5   12;
            9   10;
            10  11;
            11  12;
            12  9;
            9   11;
            10  12;
            9   13;
            10  14;
            11  15;
            12  16;
            9   14;
            10  13;
            10  15;
            11  14;
            11  16;
            12  15;
            12  13;
            9   16;
            13  14;
            14  15;
            15  16;
            16  13;
            13  15;
            14  16;
            13  17;
            14  18;
            15  19;
            16  20;
            13  18;
            14  17;
            14  19;
            15  18;
            15  20;
            16  19;
            16  17;
            13  20;
            17  18;
            18  19;
            19  20;
            20  17;
            17  19;
            18  20];
               
    % Cross-section areas (m^2)
    section = (1e-3)*(1:.5:21); % List of discrete cross-section area
    A(1:4) = section(x(1));
    A(5:12) = section(x(2));
    A(13:16) = section(x(3));
    A(17:18) = section(x(4));
    A(19:22) = section(x(5));
    A(23:30) = section(x(6));
    A(31:34) = section(x(7));
    A(35:36) = section(x(8));
    A(37:40) = section(x(9));
    A(41:48) = section(x(10));
    A(49:52) = section(x(11));
    A(53:54) = section(x(12));
    A(55:58) = section(x(13));
    A(59:66) = section(x(14));
    A(67:70) = section(x(15));
    A(71:72) = section(x(16));
    
    % Loadings
	%          node    dof      load(N)      %dof = 1 for x-axis, 2 for y-axis, 3 for z-axis
    Load{1} = [ 17      1       2e6;
                17      2       2e6;
                17      3       -2e6];
	Load{2} = [ 17      3       -2e6;
                18      3       -2e6;
                19      3       -2e6;
                20      3       -2e6];
            
    % Prescribed displacement
    %      node    dof  displacement(m)      %dof = 1 for x-axis, 2 for y-axis, 3 for z-axis
    BC = [	1       1       0;
            1       2       0;
            1       3       0;
            2       1       0;
            2       2       0;
            2       3       0;
            3       1       0;
            3       2       0;
            3       3       0;
            4       1       0;
            4       2       0;
            4       3       0];
        
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