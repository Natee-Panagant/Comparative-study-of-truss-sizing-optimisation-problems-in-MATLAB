function [f,g,rst] = f200barSI(x,varargin)
% 200-bar 2D truss optimization problem
% x = design variables (input - 0<=x<=1)
    
    % Decoding design variables
    LB = zeros(29,1); % lower-bound of design variables
    UB = 41*ones(29,1); % upper-bound of design variables
    x = LB+(UB-LB).*x; % truncate x to lower and upper-bound
    x = max(1,ceil(x)); % round x as cross-section size index
    
    % Material properties
    rho=7850; % Density (kg/m^3)
    sigma_a =400e6; %allowable stress(Pa)
    E =200e9; %modulus of elasticity(Pa)
    
    % Geometry (position in m)
    yh=fliplr([0 linspace(9,9+36,11)]);
    node=[];
    for i=1:11
        if floor(i/2)~=i/2
            node=[node;[linspace(0,24,5)' yh(i)*ones(5,1)]];
        else
            node=[node;[linspace(0,24,9)' yh(i)*ones(9,1)]];
        end
    end
    node=[node;[6 0];[18 0]];
    % element = the metrix of truss elements represented by node number combination
    %  node numbers   Element diameter  
    ele=[];
    for i=1:5
        n1=14*i-13;n2=14*i-12;n3=14*i-11;n4=14*i-10;
        n5=14*i-9;n6=14*i-8;n7=14*i-7;n8=14*i-6;
        n9=14*i-5;n10=14*i-4;n11=14*i-3;n12=14*i-2;
        n13=14*i-1;n14=14*i-0;
        n15=n14+1;n16=n15+1;n17=n16+1;n18=n17+1;n19=n18+1;
        elei1=[n1 n2
            n2 n3
            n3 n4
            n4 n5
            n1 n6
            n1 n7
            n2 n7
            n2 n8
            n2 n9
            n3 n9
            n3 n10
            n3 n11
            n4 n11
            n4 n12
            n4 n13
            n5 n13
            n5 n14];
        elei2=[n6 n7
            n7 n8
            n8 n9
            n9 n10
            n10 n11
            n11 n12
            n12 n13
            n13 n14
            n6 n15
            n7 n15
            n7 n16
            n8 n16
            n9 n16
            n9 n17
            n10 n17
            n11 n17
            n11 n18
            n12 n18
            n13 n18
            n13 n19
            n14 n19];
        ele=[ele;elei1;elei2];
    end
    ele=[ele
        71 72
        72 73
        73 74
        74 75
        71 76
        72 76
        73 76
        73 77
        74 77
        75 77];
    
    % Cross-section areas (m^2)
    section = (1e-3)*(1:.5:21); % List of discrete cross-section area
    A = zeros(200,1);
    dv1=[1, 2, 3, 4];
    dv2=[5, 8, 11, 14, 17];
    dv3=[19, 20, 21, 22, 23, 24];
    dv4=[18, 25, 56, 63, 94, 101, 132, 139, 170, 177];
    dv5=[26, 29, 32, 35, 38];
    dv6=[6, 7, 9, 10, 12, 13, 15, 16, 27, 28, 30, 31, 33, 34, 36, 37];
    dv7=[39, 40, 41, 42];
    dv8=[43, 46, 49, 52, 55];
    dv9=[57, 58, 59, 60, 61, 62];
    dv10=[64, 67, 70, 73, 76];
    dv11=[44, 45, 47, 48, 50, 51, 53, 54, 65, 66, 68, 69, 71, 72, 74,75];
    dv12=[77, 78, 79, 80];
    dv13=[81, 84, 87, 90, 93];
    dv14=[95, 96, 97, 98, 99, 100];
    dv15=[102, 105, 108, 111, 114];
    dv16=[82, 83, 85, 86, 88, 89, 91, 92, 103, 104, 106, 107, 109, 110, 112, 113];
    dv17=[115, 116, 117, 118];
    dv18=[119, 122, 125, 128, 131];
    dv19=[133, 134, 135, 136, 137, 138];
    dv20=[140, 143, 146, 149, 152];
    dv21=[120, 121, 123, 124, 126, 127, 129, 130, 141, 142, 144, 145, 147, 148, 150, 151];
    dv22=[153, 154, 155, 156];
    dv23=[157, 160, 163, 166, 169];
    dv24=[171, 172, 173,174, 175, 176];
    dv25=[178, 181, 184, 187, 190];
    dv26=[158, 159, 161, 162, 164, 165, 167, 168, 179, 180, 182, 183, 185, 186, 188,189];
    dv27=[191, 192, 193, 194];
    dv28=[195, 197, 198, 200];
    dv29=[196, 199];
    
    A(dv1)=section(x(1));
    A(dv2)=section(x(2));
    A(dv3)=section(x(3));
    A(dv4)=section(x(4));
    A(dv5)=section(x(5));
    A(dv6)=section(x(6));
    A(dv7)=section(x(7));
    A(dv8)=section(x(8));
    A(dv9)=section(x(9));
    A(dv10)=section(x(10));
    A(dv11)=section(x(11));
    A(dv12)=section(x(12));
    A(dv13)=section(x(13));
    A(dv14)=section(x(14));
    A(dv15)=section(x(15));
    A(dv16)=section(x(16));
    A(dv17)=section(x(17));
    A(dv18)=section(x(18));
    A(dv19)=section(x(19));
    A(dv20)=section(x(20));
    A(dv21)=section(x(21));
    A(dv22)=section(x(22));
    A(dv23)=section(x(23));
    A(dv24)=section(x(24));
    A(dv25)=section(x(25));
    A(dv26)=section(x(26));
    A(dv27)=section(x(27));
    A(dv28)=section(x(28));
    A(dv29)=section(x(29));
    
    % Loadings
	nf1=[1, 6, 15, 20, 29,34, 43, 48, 57, 62,71];
    nf2=[1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 15, 16, 17, 18,19, 20, 22, 24,...
        26, 28,29 30, 31, 32, 33, 34,36, 38, 40, 42, 43, 44, 45, 46, 47,...
        48, 50, 52, 54, 56, 57, 58, 59, 60, 61, 62, 64, 66, 68, 70, 71,...
        72, 73, 74, 75];
	%          node     dof             load(kips)      %dof = 1 for x-axis, 2 for y-axis, 3 for z-axis
    Load{1} = [ nf1'    ones(11,1)      (1e4)*ones(11,1);
                nf2'     2*ones(55,1)	-(1e5)*ones(55,1)];
    
    % Prescribed displacement
    %      node    dof  displacement(m)      %dof = 1 for x-axis, 2 for y-axis, 3 for z-axis
    BC = [	76      1       0
            76      2       0
            77      1       0
            77      2       0];
    
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