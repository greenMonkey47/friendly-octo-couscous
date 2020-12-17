function  de_noise = cvx_L1opt_with_relaxation(data,l1,points,n)
    

N = size(data,1);
e = ones(N,1);
D1 = spdiags([-e e], 0:1, N-1, N);


flag = ones(size(data));

flag(points) = 1000 ;
flag(points+n) = 1000;


flag1 = 100*ones(size(data));
%flag1(points+n-1) = 1;
%flag1(points+n) = 1;

flag1(points) = 0;
flag1(points-1) = 0;


flag1(points+n) = 0;
flag1(points+n-1) = 0;
%flag1(points+1) = 1;

cvx_begin
variable x(N,1)
size(data)
size(x)
size(flag)

minimize( norm((data-x).*flag(1:N),2)+l1*(norm((D1*x).*flag1(1:N-1),2)) )

cvx_end

de_noise = x;

end