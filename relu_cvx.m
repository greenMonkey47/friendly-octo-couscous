function y = relu_cvx( x , T ,p)


narginchk(3,3);

y = pow_pos(abs( x ) - T,p);
% pow_pos(x,p) : max{x,0}^p and p >=1; Convex and Non-Decreasing

%y = max((x-T).^2,0);