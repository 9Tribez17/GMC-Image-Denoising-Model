function y = soft_thresh(x, lambda)
% x为输入信号，lambda为阈值
    y = sign(x) .* max(abs(x) - lambda, 0);
end