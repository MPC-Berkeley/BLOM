

X = eye(3);
Y = eye(4);

c1 = eye(3);
c2 = ones(3,4);
c3=zeros(4,3);
c3(1,3) = 1;
c3(3,3) = 1;
c3(2,3) = 1;

% calculate trace(c1*X*c2*Y'*c3)

[A1 C1] = BLOM_ConstMatConst(c1,X,c2);
[A2 C2] = BLOM_ConstTrMatConst(eye(4),Y,c3);
[Am, Cm] = BLOM_MatMat(A1,C1,[3 4],A2,C2,[4 3],false);
[Af Cf] = BLOM_TraceMat(Am, Cm,[3 3]);

syms x1 x2 x3 x4 x5 x6 x7 x8 x9;
syms y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 y11 y12 y13 y14 y15 y16;
X = [x1 x2 x3 ; x4 x5 x6 ; x7 x8 x9].';
Y = [y1 y2 y3 y4 ;y5 y6 y7 y8 ; y9 y10 y11 y12 ; y13 y14 y15 y16 ].';

Z = expand(trace(c1*X*c2*Y.'*c3))
