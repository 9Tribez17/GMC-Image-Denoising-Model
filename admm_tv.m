%u = abs(ifft2(fft2(rhs)./(fft2(para*U_ker)+fft2(H12).*fft2(H11))));
function U = admm_tv(U,D,lambda,rho,T,P,f,H1,H11,H12,ker1)
    % define Ux where Ux* U = delta_x(U)
    for i = 1:T
        % 更新u
        rhs = -imDiff(P) + rho * imDiff(D) + imfilter(f,H1,'circular','conv');
        U = abs(ifft2(fft2(rhs)./ (fft2(H12).*fft2(H11) + fft2(rho*ker1))));
        % 更新d
        D = soft_thresh(Diff(U) + P/rho, lambda/rho);
        % 更新潜在变量
        P = P + rho * (Diff(U)-D);
    end
end