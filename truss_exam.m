function umax = truss_exam(INPUT)
%{
--------------------------------------------------------------------------
Max Ehre
October 2017

adapted from http://elementosfinitosunalmzl.wikispaces.com/
--------------------------------------------------------------------------
2D - plane truss example as given in

Lee and Kwak 2006.
--------------------------------------------------------------------------
input: nsam x 10 matrix containing nsam samples of 10-dimensional input:

E1,E2: Youngs modulus of horizontal and inclined bars resp.  [N/m^2]
A1,A2: cross section of horizontal and inclined bars resp.   [m^2]
P1-P6: vertical forces attacking in nodes 8-13. (1x6 vector) [N]
--------------------------------------------------------------------------
output: 

nsam x 1 vector containing the maximum displacements [m]
--------------------------------------------------------------------------
N: Nodes
R: Rods

       N13__R4___N12__R8___N11__R12__N10__R16__N09__R20__N09  
        /\        /\        /\        /\        /\        /\
       /  \      /  \      /  \      /  \      /  \      /  \
     R01  R3    R5   R7   R9  R11   R13 R15  R17  R19  R21  R23
     /      \  /      \  /      \  /      \  /      \  /      \
    /___R2___\/___R6___\/__R10___\/__R14___\/__R18___\/__R22___\
  N01        N02      N03       N04        N05       N06       N07
--------------------------------------------------------------------------
%}

A1 = INPUT(:,1);
A2 = INPUT(:,2);
E1 = INPUT(:,3);
E2 = INPUT(:,4);
P  = -INPUT(:,5:10);

n_sam = size(INPUT,1);
umax = zeros(n_sam,1);

%% initial input variables
nfe  = 23;          % number of elements (bars)
ned  = 2;          % number of dof per node
nnp  = 13;          % number of nodal points
ndof = ned*nnp;    % number of degrees of freedom (dof)

%% element, nodes and dofs association
% IEN: connectivity matrix, nfe x nodes
IEN = [1 13;   % bar 1 has nodes 1 and 3
       1 2;   % bar 2 has nodes 1 and 4 ...
       13 2;   
       13 12; 
       2 12; 
       2 3;
       12 3;
       12 11;
       3 11;
       3 4;
       11 4;
       11 10;
       4 10;
       4 5;
       10 5;
       10 9;
       5 9;
       5 6;
       9 6;
       9 8;
       6 8;
       6 7;
       8 7]; 

% ID: destination array, nodes vs dof 
ID   = [1 2;  % node 1 has dof 1 and 2
        3 4;  % node 2 has dof 2 and 3 ...
        5 6;  
        7 8;
        9 10;
        11 12;
        13 14;
        15 16;
        17 18;
        19 20;
        21 22;
        23 24;
        25 26]; 

% LM: localization matrix, nfe x dof 
LM = cell(nfe,1);
for e = 1:nfe
   LM{e} = [ID(IEN(e,1),:), ID(IEN(e,2),:)]; %  element 1 has dof 1 2 5 6 ...
end

% deterministic rod properties
ang   = atan2(200,200)*180/pi;                                          % inclination angle of the truss [deg]
theta = repmat([ang     0       -ang    0],1,6); theta(end) = [];       % inclination angle [deg]
leng  = repmat(4*[1/sqrt(2) 1],1,12); leng(end) = [];                   % bar length [m]

%% MC loop

for i=1:n_sam
    
    % stochastic rod properties
    area  = repmat([A2(i) A1(i)],1,12); area(end) = []; % bar cross sectional area [m2]
    E     = repmat([E2(i) E1(i)],1,12); E(end) = [];    % young's modulus [N/m^2]

    % material properties
    k = E.*area./leng;  % stiffness of each bar

    %% stiffness matrix assembly
    K = zeros(ndof);      % stiffness matrix
    T = cell(nfe,1);   % transformation matrix
    for e = 1:nfe 
       % sin and cos of the inclination
       c = cosd(theta(e)); 
       s = sind(theta(e));  

       % coordinate transformation matrix for the bar e
       T{e} = [ c  s  0  0;
               -s  c  0  0;
                0  0  c  s; 
                0  0 -s  c ];

       % local stiffness matrix for the bar e     
       Kloc = k(e)*[ 1  0 -1  0;  
                     0  0  0  0;            
                    -1  0  1  0;  
                     0  0  0  0 ];

       % global stiffness matrix assembly
       K(LM{e},LM{e}) = K(LM{e},LM{e}) + T{e}'*Kloc*T{e}; 
    end

    %% boundary conditions
    % applied loads
    fc = zeros(ndof,1);
    fc(16:2:26) = -P(i,:); % [N]
    
    % supports and free nodes
    c = [1 2 14];             % fixed dofs 
    d = setdiff(1:ndof,c);   % free dofs 

    % remove supports and free nodes from force vector
    fc(c) = [];
    
    % f = equivalent nodal force vector
    % q = equilibrium nodal force vector
    % a = displacements

    %| qd |   | Kcc Kcd || ac |   | fd |
    %|    | = |         ||    | - |    |
    %| qc |   | Kdc Kdd || ad |   | fc |

    %% solve the system of equations
    Kcc = K(c,c); Kcd = K(c,d);
    Kdc = K(d,c); Kdd = K(d,d);

    ac = [0; 0; 0];       % displacements for c are 0
    ad = Kdd\(fc-Kdc*ac); % = linsolve(Kdd, fc-Kdc*ac)
    qd = Kcc*ac + Kcd*ad;

    % 
    a = zeros(ndof,1); q = zeros(ndof,1); 
    a(c) = ac;       q(c) = qd;
    a(d) = ad;     % q(d) = qc = 0

    %% axial loads
    N = zeros(nfe,1);
    for e = 1:nfe
       N(e) = k(e)*[-1 0 1 0]*T{e}*a(LM{e});
    end
    
    umax(i,1) = a(8);
    
end

return

