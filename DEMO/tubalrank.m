function trank = tubalrank(X,tol)

X = fft(X,[],3);
[n1,n2,n3] = size(X);
s = zeros(min(n1,n2),1);

s = s + svd(X(:,:,1),'econ');

halfn3 = round(n3/2);
for i = 2 : halfn3
    s = s + svd(X(:,:,i),'econ')*2;
end

if mod(n3,2) == 0
    i = halfn3+1;
    s = s + svd(X(:,:,i),'econ');
end
s = s/n3;

if nargin==1
   tol = max(n1,n2) * eps(max(s));
end
trank = sum(s > tol);
