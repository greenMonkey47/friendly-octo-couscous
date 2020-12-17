function cuts_in = cut_detect_cvx(data,cut_dist,n,n3,k)

N = size(data,1);

e = ones(N,1);
temp = [-e repmat(e*0,[1 n-1]) e];
D1 = spdiags([-e repmat(e*0,[1 n-1]) e], 0:n, N-n+1, N);
e1 = ones(N,k);
D2 = spdiags(e1, 0:k-1, N-k+1, N);
D4 = spdiags([repmat(-e/n3,[1 n3]) e*0 repmat(e/n3,[1 n3])], 0:2*n3, N-2*n3, N);


l1 = 1;
l2 = cut_dist;

tic
cvx_begin
variable x(N,1)
variable g(N,1)

variable temp(N,1)

%minimise( l1*norm(((D1*data).*x(1:N-n)),1) + norm(l2.*(-x),1))

%minimise( norm(l1*abs(D1*data).*x(1:N-n) + (1-x(1:N-n))*(l2),1))
minimize(...
    norm(l1*abs(D1*data).*(1-x(1:N-n+1)) + (x(1:N-n+1))*(l2) , 1))% +...
%norm(l1*abs(D4*data).*(1-x(2*n3+1:N)) + (x(2*n3+1:N))*(l2) , 1) ...
%)
%5norm(l1*abs(D4*data).*(1-x(1:N-n3+1)) + (x(1:N-n3+1))*(l2) , 1) + ...

%lambda0*(sum_log(x)) + lambda0*(sum_log((g-array).^2)) + ...
%lambda1*(sum_log(1-x) + sum_log((D1p*g).^2)) + ...
%lambda2*(sum_log(1-x) + sum_log((D2p*g).^2)) + ...
%lambda3*(sum_log(1-x) + sum_log((D3p*g).^2)) ...


subject to
x >= 0 ;
x <= 1 ;
(D2*x) <= 1;
abs(D4*data).*(x(n3+1:N-n3)) + (1-x(n3+1:N-n3))*(l2) >= l2;

cvx_end

cuts_gaze=x;
cuts_in = find(cuts_gaze > 0.5);
toc

% figure,plot(data,'.b')
% hold on;
% scatter(cuts_in,data(cuts_in),20,'*r')
% plot(abs(D1*data))
% plot(1:3600,ones(3600,1)*l2,'-k');
% plot(n3+1:N-n3,abs(D4*data.*(x(n3+1:N-n3))),'-m');
% plot(n3+1:N-n3,abs(D4*data),'-g');
% plot(n3+1:N-n3,abs(D4*data).*(x(n3+1:N-n3)) + (1-x(n3+1:N-n3))*(l2)-l2,'-g')
% axis([0 l 0 1366])
