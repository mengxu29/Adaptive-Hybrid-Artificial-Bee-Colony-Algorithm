
function y = exp_q(x, q)
if (q == 1)
    y = exp(x);
else
    temp = (1 + (1 - q) * x);
    mask = temp > 0;
    m = 0e-40; % this small m is very small, close to zero, which is to make the log function meaningful
    y = ones(1, length(x)) * m;
    y(mask) = temp(mask).^(1 / (1 - q));
end
end