 function [y1,y2] = TGV2L2_ADMM(I,Bnim,H,alpha,miter)
%% Solve TGV2-(shearlet)L1-L2 model:
% $\frac\beta2\|Ku-b\|_2^2+\lambda\|SH(u)\|_1+TGV_\alpha^2(u)$
% In this version, we solve the (u,p)-subproblem directly.
%
% Inputs:
% I - ground truth image
% P - full sampling matrix
% B - full partial k-space data
% para - parameter structure
%     .maxiter
%     .alpha - 2X1 vector
%     .lambda/.beta/
%     .mu - 3X1 vector storing $\mu_1$, $\mu_2$ and $\mu_3$
%     .gamma - scalar, $\mu$ in the paper
%     .skrinkopt - shrink method
%     .period/.contfact - continuation parameters
%     .correctopt - correction option of results
%     .kernels/trans - sparsifying kernels/transform
%     .initial.* - initial intermediate results
%
% Outputs:
% output.y - reconstructed image
%       .err/.enderr/.x/.y/.z/.xx/.yy/.zz
%
% Refer to the paper ``A New Detail-preserving Regularity Scheme" by 
%     Weihong Guo, Jing Qin and Wotao Yin.
%
% Jing, updated on 07-16-2014

%% Read input parameters
 global m n
[m,n]            = size(I);

eta              = alpha;
% eta(2)           = 1;
gamma            = 1./eta;

a                = alpha(1)/eta(1);
b                = alpha(2)/eta(2); 

    
%     u  = imfilter(Bnim,H1,'circular','conv');
%     du(:,:,1) = Dx(u);
%     du(:,:,2) = Dy(u);
%     p = du;
    u  = zeros(m,n);
    p  = zeros(m,n,2);
    y  = zeros(m,n,2);
    du = zeros(m,n,2);
    z  = zeros(m,n,4);
    yy = zeros(m,n,2);
    zz = zeros(m,n,4);
    Epp = zeros(m,n,4);
    
%  B = fft2(u);
err = zeros(miter,1);
errpre = 1;

D1 = abs(psf2otf([-1,1],[m,n])).^2;
D2 = abs(psf2otf([-1;1],[m,n])).^2;
K = abs(psf2otf(H,[m,n])).^2;

d1 = a.*(D1+D2)+K;
d2 = a+b*(D1+.5*D2);
d3 = a+b*(D2+.5*D1);
d4 = psf2otf([1, -1],[m, n]);
d5 = psf2otf([1; -1],[m, n]);
d6 = d5(:,1)*(conj(d4(1,:)));
% d6 = d5(:,1)*d5(:,1)';

d4 = -a*d4;
d5 = -a*d5;
d6 = .5*b*d6;

d4t = conj(d4);
d5t = conj(d5);
d6t = conj(d6);

B = fft2(imfilter(Bnim,H,'circular','conv'));

denom = d1.*d2.*d3+d4t.*d6t.*d5+d4.*d5t.*d6...
    -d2.*d5.*d5t-d1.*d6t.*d6-d3.*d4.*d4t;

psnrM = [];
%% Main loop
% h = waitbar(0,'Please wait ...');
for i = 1:miter
    u2 = u;
    
    % Step 1: shrinkage of shearlets
%     for j = 1:N
%         x(:,:,j) = shrinkmethod(uband(:,:,j)+xx(:,:,j),1/mu(1));
%     end
    
    % Step 2: shrinkage of gradients
%     y = shrink2(du-p-eta(1)*yy,eta(1));
%     y = p_shrinkage_b(du-p-eta(1)*yy,eta(1),1);
      y = shrinkage_a(du-p-eta(1)*yy,eta(1));
    
    % Step 3: shrinkage of second order derivatives
%     Epp = Ep(p);
%     z = shrink2(Epp-eta(2)*zz,eta(2));
      z = shrinkage_a(Epp-eta(2)*zz,eta(2));
%       z = zeros(m,n,4);
    
    % Solve the (u,p1,p2)-subproblem
%     RHS1 = 0;
%     for j = 1:N
%         RHS1 = RHS1 + H(:,:,j).*fft2(x(:,:,j)-xx(:,:,j));
%     end

    FB1 = a*(fft2(Dxt(y(:,:,1)+eta(1)*yy(:,:,1)))...
        +fft2(Dyt(y(:,:,2)+eta(1)*yy(:,:,2))))+B;
    FB2 = -a*fft2(y(:,:,1)+eta(1)*yy(:,:,1))+b*(fft2(Dxt(z(:,:,1)+eta(2)*zz(:,:,1)))...
        +fft2(Dyt(z(:,:,3)+eta(2)*zz(:,:,3))));
    FB3 = -a*fft2(y(:,:,2)+eta(1)*yy(:,:,2))+b*(fft2(Dyt(z(:,:,2)+eta(2)*zz(:,:,2)))...
        +fft2(Dxt(z(:,:,3)+eta(2)*zz(:,:,3))));
    
    % Step 4: Solve for u
    RHS1 = (d2.*d3-d6.*d6t).*FB1-(d3.*d4t-d6.*d5t).*FB2...
        +(d4t.*d6t-d2.*d5t).*FB3;
     u = ifft2(RHS1./denom);
     u = abs(u);
     du(:,:,1) = Dx(u);
     du(:,:,2) = Dy(u);
    
    % Record the error
    err(i) = norm(u-I,'fro')/norm(I,'fro');

    % Step 5: Solve for p1
    RHS2 = (d5.*d6t-d3.*d4).*FB1+(d1.*d3-d5.*d5t).*FB2...
        +(d4.*d5t-d1.*d6t).*FB3;
    p(:,:,1) = ifft2(RHS2./denom);
        
    % Step 6: Solve for p2
    RHS3 = (d4.*d6-d2.*d5).*FB1+(d4t.*d5-d1.*d6).*FB2...
        +(d1.*d2-d4.*d4t).*FB3;
    p(:,:,2) = ifft2(RHS3./denom);

    % Step 7: Bregman update
%     uband = trans(u);
%     xx = xx + gamma*(uband-x);
   
    Epp = Ep(p);
    yy = yy + gamma(1)*(-du+y+p);
    zz = zz + gamma(2)*(z-Epp);
    
    
    % Stopping criterion
    relerr = norm(u(:)-u2(:))/norm(u(:));
%     if relerr < 1e-5 
%         break
%     elseif errpre < err(i)
%         break
%     else 
%         errpre = err(i);
%     end
%    
%     % Display the progress
%     waitbar(i/miter,h,sprintf('iteration %d',i));
    psnr_iter = psnr(u,I);
    psnrM = [psnrM;psnr_iter];
%     disp(['Iteration k = ' num2str(i), ', PSNR = ' num2str(psnr_iter)]);
end
% close(h)

if i<miter
    err(i+1:end) = [];
end

%% Output results
y1 = u;
y2 = psnrM;
% output.u = u;
% output.psnr = psnrM;

% output.u        = abs(u);
% output.err      = err;
% output.enderr   = err(end);
% output.y        = y;
% output.z        = z;
% output.yy       = yy;
% output.zz       = zz;
% output.tv1      = sum(sum(sum(abs(du-p),3)));
% output.tv2      = sum(sum(sum(abs(Ep(p)),3)));

function z = Ep(p)
global m n
z = zeros(m,n,4);
z(:,:,1) = Dx(p(:,:,1));
z(:,:,2) = Dy(p(:,:,2));
z(:,:,3) = (Dy(p(:,:,1))+Dx(p(:,:,2)))./2;
z(:,:,4) = z(:,:,3);

function dxu = Dx(U)
% x-axis forward difference
dxu = U(:,[2:end,1])-U;

function dyu = Dy(U)
% y-axis forward difference
dyu = U([2:end,1],:)-U;           

function dxtu = Dxt(U)
% -(x-axis backward difference)
dxtu = U(:,[end,1:end-1])-U;

function dytu = Dyt(U)
% -(y-axis backward difference)
dytu = U([end,1:end-1],:)-U;