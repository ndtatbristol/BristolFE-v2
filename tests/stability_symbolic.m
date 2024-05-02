clear

%test basic procedure - mass spring conservative system

K = sym('K', [1,1]);
M = sym('M', [1,1]);
inv_M = 1 ./ M;
C = sym('C', [1,1]);
C(1,1)= 0;
dt = sym('dt');

A2 = [1 - dt ^ 2 * inv_M * K, dt - dt * inv_M * C
    -dt * inv_M * K, 1 - dt * inv_M * C];

X = [K, 0; 0, M];

s = [1;0];
s.' * (X - A2.' * X * A2) * s


%interface test

clear
K = sym('K', [2,2]);
K(1,2) = 0;
K(2,1) = 0;

M = sym('M', [2,2]);
inv_M = 1 ./ M;
inv_M(1,2) = 0;
inv_M(2,1) = 0;
M(1,2) = 0;
M(2,1) = 0;

C = sym('C', [2,2]);
C(1,2) = -1;
C(2,1) = -1;
C(1,1) = 0;
C(2,2) = 0;

dt = sym('dt');

B1 = 2 * eye(2) - dt * (inv_M * C) - dt ^ 2 * (inv_M * K);
B2 = dt * (inv_M * C) - eye(2);

A = [B1,         B2
             speye(2), zeros(2)];

A2 = [eye(2) - dt ^ 2 * inv_M * K, dt * eye(2) - dt * inv_M * C
    -dt * inv_M * K, eye(2) - dt * inv_M * C];

s1 = sym('s1', [4, 1]);
s2 = A2 * s1;

X = [K, zeros(2); zeros(2), M];

E1 = s1.' * X * s1;
E2 = s2.' * X * s2;

dE = simplify(E2 - E1)

[1,0,0,0] * A2.' * X * A2 * [1,0,0,0].'