function R = rot_mat(a,b)

GG = @(A,B) [ dot(A,B) -norm(cross(A,B)) 0;...
              norm(cross(A,B)) dot(A,B)  0;...
              0              0           1];

FFi = @(A,B) [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];

RR = @(Fi,G) Fi*G*inv(Fi);
R = RR(FFi(a,b), GG(a,b));
 