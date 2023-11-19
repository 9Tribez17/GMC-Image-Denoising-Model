function [PSNR0, U] = GMC_ADMM(Img,f,H,lambda1,rho1,rho2,gamma)
%参数设置
%lambda = 0.04; %正则化参数
%lambdaArr = (0.1:0.01:0.2);
T = 100; % 最大迭代次数
%rho1 = 32;
%rho2 = 32;
%gamma = 0.25;
%gammaArr  = (0.4:0.01:0.5);

% Generate the blurring filter
[m1,n1] = size(H);
%Generate the filter for the transform of the blurring operator
[m1,n1] = size(H); 
A = zeros(m1,n1);
H1 = zeros(m1,n1);
for i = 1 : m1
    A(i,:) = H(m1+1-i,:);
end
for j = 1 : n1
    H1(:, n1+1-j) = A(:, j);
end

[m,n]=size(f); 
mn = min(m,n);
H11 = zeros(m,n);
H12 = zeros(m,n);
H11(1:(m1+1)/2, 1:(n1+1)/2)=H((m1+1)/2:m1,(n1+1)/2:n1);
H11(m-(m1-3)/2:m,1:(n1+1)/2)=H(1:(m1-1)/2,(n1+1)/2:n1);
H11(1:(m1+1)/2, n-(n1-3)/2:n)=H((m1+1)/2:m1,1:(n1-1)/2);
H11(m-(m1-3)/2:m,n-(n1-3)/2:n) = H(1:(m1-1)/2,1:(n1-1)/2);

H12(1:(m1+1)/2, 1:(n1+1)/2)=H1((m1+1)/2:m1,(n1+1)/2:n1);
H12(m-(m1-3)/2:m,1:(n1+1)/2)=H1(1:(m1-1)/2,(n1+1)/2:n1);
H12(1:(m1+1)/2, n-(n1-3)/2:n)=H1((m1+1)/2:m1,1:(n1-1)/2);
H12(m-(m1-3)/2:m,n-(n1-3)/2:n) = H1(1:(m1-1)/2,1:(n1-1)/2);

% 初始化变量
U = imfilter(f,H,'circular','conv'); 
V = zeros(m,n); 
D = zeros(m*2,n);
P = zeros(m*2,n);
q1 = zeros(m*2,n);
q2 = zeros(m*2,n);

% Calculate (Det'*Det)_ker
k1 = [0 -1 0; -1 4 -1; 0 -1 0];
[m2,n2] = size(k1);
ker1(1:(m2+1)/2, 1:(n2+1)/2) = k1((m2+1)/2:m2,(n2+1)/2:n2);
ker1(m-(m2-3)/2:m,1:(n2+1)/2) = k1(1:(m2-1)/2,(n2+1)/2:n2);
ker1(1:(m2+1)/2, n-(n2-3)/2:n) = k1((m2+1)/2:m2,1:(n2-1)/2);
ker1(m-(m2-3)/2:m,n-(n2-3)/2:n) = k1(1:(m2-1)/2,1:(n2-1)/2);


% 代入模型
U = admm_gmc(U,V,D,P,lambda1,rho1,rho2,T,f,q1,q2,gamma,H,H1,H11,H12,ker1);
%psnrAr = zeros(1,10);
%for i = 1:10
%   U = admm_gmc(U,V,D,P,lambda,rho1,rho2,T,f,q1,q2,gammaArr(1,i),H,H1,H11,H12,ker1);
%   psnrAr(1,i) = psnr(U,Prev_f)
%end

% 显示结果
PSNR0 = psnr(U,Img);
%PSNR1 = psnr(f,Prev_f);
%figure;
%subplot(1,3,1); imshow(f); title(['Noisy Image PSNR=' num2str(PSNR1)]);
%subplot(1,3,2); imshow(U); title(['Denoised Image PSNR=' num2str(PSNR0)]);
%imshow(U); 
%subplot(1,3,3); imshow(Prev_f); title('Previous Image');

%fprintf('Index      | PSNR     |\n'); 
%fprintf('Rainy      | %.4f  |\n',PSNR0);
%fprintf('Noisy      | %.4f  |\n',PSNR1);
end



