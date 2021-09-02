function [ y ] = s_shrinkage_b_mu( b,mu )

y=sign(b).*max(zeros(length(b),1),abs(b)-ones(length(b),1)*mu);

end

