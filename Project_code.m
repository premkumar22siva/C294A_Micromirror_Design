% This MATLAB script is written by Premkumar Sivakumar (UID: 506367127),
% to calculate the mode shapes and natural frequencies of the flexure system
% synthesized for the resonant translating micromirror

%% 
% Assumed material properties
% for [100] direction single crystal silicon
E = 130 * 10^9; % Modulus of elasticity
G = 50 * 10^9; % Modulus of rigidity

Delta = [zeros(3,3) eye(3,3);
         eye(3,3) zeros(3,3)];

%% Flexure Body 1 (Blade flexure)

% Given dimensions
l1 = 1 * 10^-3; % length
b1 = 1 * 10^-3; % breadth
h1 = 0.01 * 10^-3; % height

% Location & unit vectors
L1 = [0.0015 0.001 -0.000005]';
n11 = [1 0 0]';
n12 = [0 0 -1]';
n13 = [0 1 0]';

% Transformation matrix
N1 = [n11 n12 n13 zeros(3,3); 
      cross(L1,n11) cross(L1,n12) cross(L1,n13) n11 n12 n13];

% To compute polar moment of inertia, J2
sum = 0;
N = 1000;
for n = 1:2:N
    sum = sum + (1/n^5) * tanh(n*b1*pi/(2*h1));
end
J1 = ((h1^3 * b1) / 3) * (1 - (192 / pi^5) * (h1/b1) * sum);

% Moment of inertia, I
Ix1 = b1*h1^3/12; 
Iy1 = h1*b1^3/12;
A1 = b1*h1;

% Temporary variables
u = l1/(E*Ix1);
v = -l1^2/(2*E*Ix1);
w = l1/(E*Iy1);
x = l1^2/(2*E*Iy1); 
y = l1/(G*J1);
z = l1^3/(3*E*Iy1);
r = l1^3/(3*E*Ix1);
s = l1/(E*A1);

% Compliance matrix
C1 = [u 0 0 0 v 0;
      0 w 0 x 0 0;
      0 0 y 0 0 0;
      0 x 0 z 0 0;
      v 0 0 0 r 0;
      0 0 0 0 0 s];

% Stiffness matrix
S1 = C1^-1;

% Twist-Wrench Stiffness matrix
K1 = N1 * Delta * S1 * N1^-1;

%% Flexure Body 2 (Blade flexure)

% Given dimensions
l2 = 1 * 10^-3;
b2 = 1 * 10^-3;
h2 = 0.01 * 10^-3;

% Location & unit vectors
L2 = [0.001 0.0015 -0.000195]';
n21 = [0 -1 0]';
n22 = [0 0 -1]';
n23 = [1 0 0]';

% Transformation matrix
N2 = [n21 n22 n23 zeros(3,3); 
      cross(L2,n21) cross(L2,n22) cross(L2,n23) n21 n22 n23];

% To compute polar moment of inertia, J2
sum = 0;
N = 1000;
for n = 1:2:N
    sum = sum + (1/n^5) * tanh(n*b2*pi/(2*h2));
end
J2 = ((h2^3 * b2) / 3) * (1 - (192 / pi^5) * (h2/b2) * sum);

% Moment of inertia, I
Ix2 = b2*h2^3/12; 
Iy2 = h2*b2^3/12;
A2 = b2*h2;

% Temporary variables
u = l2/(E*Ix2);
v = -l2^2/(2*E*Ix2);
w = l2/(E*Iy2);
x = l2^2/(2*E*Iy2); 
y = l2/(G*J2);
z = l2^3/(3*E*Iy2);
r = l2^3/(3*E*Ix2);
s = l2/(E*A2);

% Compliance matrix
C2 = [u 0 0 0 v 0;
      0 w 0 x 0 0;
      0 0 y 0 0 0;
      0 x 0 z 0 0;
      v 0 0 0 r 0;
      0 0 0 0 0 s];

% Stiffness matrix
S2 = C2^-1;

% Twist-Wrench Stiffness matrix
K2 = N2 * Delta * S2 * N2^-1;

%% Flexure Body 3 (Blade flexure)

% Given dimensions
l3 = 1 * 10^-3;
b3 = 1 * 10^-3;
h3 = 0.01 * 10^-3;

% Location & unit vectors
L3 = [0.0015 0.002 -0.000005]';
n31 = [-1 0 0]';
n32 = [0 0 -1]';
n33 = [0 -1 0]';

% Transformation matrix
N3 = [n31 n32 n33 zeros(3,3); 
      cross(L3,n31) cross(L3,n32) cross(L3,n33) n31 n32 n33];

% To compute polar moment of inertia, J2
sum = 0;
N = 1000;
for n = 1:2:N
    sum = sum + (1/n^5) * tanh(n*b3*pi/(2*h3));
end
J3 = ((h3^3 * b3) / 3) * (1 - (192 / pi^5) * (h3/b3) * sum);

% Moment of inertia, I
Ix3 = b3*h3^3/12; 
Iy3 = h3*b3^3/12;
A3 = b3*h3;

% Temporary variables
u = l3/(E*Ix3);
v = -l3^2/(2*E*Ix3);
w = l3/(E*Iy3);
x = l3^2/(2*E*Iy3); 
y = l3/(G*J3);
z = l3^3/(3*E*Iy3);
r = l3^3/(3*E*Ix3);
s = l3/(E*A3);

% Compliance matrix
C3 = [u 0 0 0 v 0;
      0 w 0 x 0 0;
      0 0 y 0 0 0;
      0 x 0 z 0 0;
      v 0 0 0 r 0;
      0 0 0 0 0 s];

% Stiffness matrix
S3 = C3^-1;

% Twist-Wrench Stiffness matrix
K3 = N3 * Delta * S3 * N3^-1;

%% Flexure Body 4 (Blade flexure)

% Given dimensions
l4 = 1 * 10^-3;
b4 = 1 * 10^-3;
h4 = 0.01 * 10^-3;

% Location & unit vectors
L4 = [0.002 0.0015 -0.000195]';
n41 = [0 1 0]';
n42 = [0 0 -1]';
n43 = [-1 0 0]';

% Transformation matrix
N4 = [n41 n42 n43 zeros(3,3); 
      cross(L4,n41) cross(L4,n42) cross(L4,n43) n41 n42 n43];

% To compute polar moment of inertia, J2
sum = 0;
N = 1000;
for n = 1:2:N
    sum = sum + (1/n^5) * tanh(n*b4*pi/(2*h4));
end
J4 = ((h4^3 * b4) / 3) * (1 - (192 / pi^5) * (h4/b4) * sum);

% Moment of inertia, I
Ix4 = b4*h4^3/12; 
Iy4 = h4*b4^3/12;
A4 = b4*h4;

% Temporary variables
u = l4/(E*Ix4);
v = -l4^2/(2*E*Ix4);
w = l4/(E*Iy4);
x = l4^2/(2*E*Iy4); 
y = l4/(G*J4);
z = l4^3/(3*E*Iy4);
r = l4^3/(3*E*Ix4);
s = l4/(E*A4);

% Compliance matrix
C4 = [u 0 0 0 v 0;
      0 w 0 x 0 0;
      0 0 y 0 0 0;
      0 x 0 z 0 0;
      v 0 0 0 r 0;
      0 0 0 0 0 s];

% Stiffness matrix
S4 = C4^-1;

% Twist-Wrench Stiffness matrix
K4 = N4 * Delta * S4 * N4^-1;

%% Stage 1

% Given material property
rho = 2330; % Silicon

% Given dimensions
l = 0.001;
b = 0.001;
h = 0.0002;
V = l*b*h;

% Location & unit vectors
Lm = [0.0015 0.0015 -0.0001]';
n1m = [1 0 0]';
n2m = [0 1 0]';
n3m = [0 0 1]';

% Transformation matrix
Nm = [n1m n2m n3m zeros(3,3); 
      cross(Lm,n1m) cross(Lm,n2m) cross(Lm,n3m) n1m n2m n3m];

% Mass moments of inertia
Ixm = (rho*V)*(b^2 + h^2)/12; 
Iym = (rho*V)*(l^2 + h^2)/12;
Izm = (rho*V)*(b^2 + l^2)/12;
inertia = [Ixm Iym Izm rho*V rho*V rho*V]';
In = diag(inertia);

%% Computing Mode Shapes and Natural Frequencies

K_TW = K1 + K2 + K3 + K4
M_TW = Nm * Delta * In * Nm^-1

[ModeShapes, EigenValues] = eig(inv(M_TW)*K_TW)

Nat_freq = sqrt(EigenValues)

