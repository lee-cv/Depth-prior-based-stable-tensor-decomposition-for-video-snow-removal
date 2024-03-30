
clear;
clc;
img = './test.jpg';
I = double(imread(img));

I = I/255;
maxP = max(abs(I(:)));
[n1,n2,n3] = size(I);
In = I;
rhos = 0.3;
ind = find(rand(n1*n2*n3,1)<rhos);
In(ind) = rand(length(ind),1);

opts.mu = 1e-4;
opts.tol = 1e-6;
opts.rho = 1.1;
opts.max_iter = 1000;
opts.DEBUG = 1;

[n1,n2,n3] = size(In);
lambda = 1/sqrt(max(n1,n2)*n3);
[hat,E,err,iter] = trpca_tnn(In,lambda,opts);
 
hat = max(hat,0);
hat = min(hat,maxP);
psnr = PSNR(I,hat,maxP);

imshow(hat/max(hat(:)))





