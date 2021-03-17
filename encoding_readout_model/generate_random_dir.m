function sign_dir = generate_random_dir(direction0, gamma, n_samples)


N=length(direction0);
sign_dir = zeros(N,n_samples);
% gamma_est = zeros(1,n_samples);

s0 = eye(N);
p = mvnrnd(zeros(N,1),s0,n_samples)';
for i = 1:n_samples
    w = p(:,i)-((p(:,i)'*direction0)/(direction0'*direction0))*direction0;
    w = norm(direction0)*w/norm(w);    
    sign_dir(:,i) = cos(gamma)*direction0 + sin(gamma)*w;
%     gamma_est(i) = md_angle(v,direction0);
end

% figure; plot(gamma_est)