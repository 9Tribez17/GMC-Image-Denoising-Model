function U = admm_gmc(U,V,D,P,lambda1,rho1,rho2,T,f,q1,q2,gamma,H,H1,H11,H12,ker1)
    for i = 1:T
        % 更新u
        Av = imfilter(V,H,'circular','conv');
        rhs1 = -imDiff(q1) + rho1 * imDiff(D) - gamma * imfilter(Av,H1,'circular','conv') +imfilter(f,H1,'circular','conv');
        U = abs(ifft2(fft2(rhs1)./(fft2(1-gamma) .* fft2(H12).* fft2(H11)+ fft2(rho1 * ker1))));
        % 更新V
        Au = imfilter(U,H,'circular','conv');
        rhs2 = -imDiff(q2) + rho2 * imDiff(P) + gamma * imfilter(Au,H1,'circular','conv');
        V = abs(ifft2(fft2(rhs2)./(fft2(rho2 * ker1) + fft2(gamma * H12).*fft2(H11))));
        % 更新D
        D = soft(Diff(U)+q1/rho1,lambda1/rho1);
        % 更新P
        P = soft(Diff(V)+q2/rho2,lambda1/rho2);        
        % 更新潜在变量
        q1 = q1 + rho1 * (Diff(U) -D);
        q2 = q2 + rho2 * (Diff(V) -P);
    end
end