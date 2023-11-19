function y = shrinkage_a(grad_aup,eta)

   [i,j,L] = size(grad_aup);
   y = zeros(i,j,L);
%    n = n/2;
%    grad = [];
%    for l = 1 : L
%        grad = [grad,grad_aup(:,:,l)];
%    end
%    grad = [grad_aup(:,:,1),grad_aup(:,:,2)];
    sn = normgrad_a(grad_aup);
%   sn = normgrad(grad_aup,0);
%     sn = sn.*(sn>eta)+eta*(ones(m,n)-(sn>eta));
%     grad_aun = [u_ax,u_ay];
    msn = max(1-eta./sn,0);
    for l = 1 : L
    y(:,:,l) = grad_aup(:,:,l).*msn;
    end
    