function [opt_data, vc1_opt, vc2_opt,z] = path_optimization_cvx(array,bool,tbool,lambda0,lambda1,lambda2,lambda3,vc1,vc2,thresh,out_width,per_frameVar,flag)

N = size(array,1);
e = ones(N,1);
D1 = spdiags([-e e], 0:1, N-1, N);
Dt = spdiags([-e e],0:1,N-1,N);
D2 = spdiags([e -2*e e], 0:2, N-2, N);
D3 = spdiags([-e 3*e -3*e e], 0:3, N-3, N);

in = find(bool(:,2)==0);
D1(in,:)=0;
in = find(bool(:,3)==0);
D2(in,:)=0;
in = find(bool(:,4)==0);
D3(in,:)=0;

D1z = spdiags([-e e], 0:1, N-1, N);
D2z = spdiags([e -2*e e], 0:2, N-2, N);
D3z = spdiags([-e 3*e -3*e e], 0:3, N-3, N);

in = find(tbool(:,2)==0);
D1z(in,:)=0;
in = find(tbool(:,3)==0);
D2z(in,:)=0;
in = find(tbool(:,4)==0);
D3z(in,:)=0;



% cvx_begin
% variable g(N,1)
% 
% minimise(lambda0*sum_square(bool(1:N,1).*(g(1:N)-array(1:N)))...
%     + lambda1*norm( (D1*g),1)+ lambda2*norm((D2*g),1) + lambda3*norm((D3*g),1)...
%     )
% subject to
% abs(D1*g) <= vc1;
% abs(D2*g) <= vc2;
% cvx_end


if flag == 0
    z = ones(N,1);

    cvx_begin
    variable g(N,1)    
    minimise(lambda0*sum((bool(1:N,1).*(relu_cvx(g(1:N)-array(1:N),thresh,2))))...
        + lambda1*norm( (D1*g),1)+ lambda2*norm((D2*g),1) + lambda3*norm((D3*g),1)...
        )
    subject to
    abs(D1*g) <= vc1;
    abs(D2*g) <= vc2;
    %abs(Dt*g) >= out_width*0.7;
    % g - (out_width/2) >= 1;
    % g + (out_width/2) <= 1366;
    cvx_end
else
    
    cvx_begin
    variable g(N,1);
    variable temp(N,1);
    variable z(N,1);
    
    minimise(lambda0*sum((bool(1:N,1).*(relu_cvx(g(1:N)-array(1:N),thresh,2))))...
        + lambda1*norm( (D1*g),1)+ lambda2*norm((D2*g),1) + lambda3*norm((D3*g),1)...
        + lambda0*5*sum_square(tbool(1:N,1).*(z(1:N)-per_frameVar(1:N)))+ lambda1*norm((D1*z),1)+ lambda3*norm((D2*z),1) + lambda3*norm((D3*z),1)...
         )
    
    subject to
    abs(D1*g) <= vc1;
    abs(D2*g) <= vc2;
    abs(D1z*z) <= 3;
    cvx_end

end

opt_data=g;
vc1_opt = abs(D1*g);
vc2_opt = abs(D2*g);
