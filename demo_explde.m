% demo_explde.m
% Code examples for explde.m
% Ver 09-Nov-2026

% Solve for F(t): F'(t) = A(t)*F(t), F(0) = F0

% Constant-A case: F(t) = expm(A*t)*F0
% Test case: Compute expm(A) with A = [a,b;b,a].
% expm([a,b;b,a]) = exp(a)*[cosh(b),sinh(b);sinh(b),cosh(b)]
a = 1;
b = 2;
exp_ = exp(a)*[cosh(b),sinh(b);sinh(b),cosh(b)];
err = explde([a,b;b,a])-exp_; err = max(abs(err(:)));
disp(['Constant-A explde error: ' num2str(err)]) % 1.7764e-15

% Variable-A (function handle) test case:
% Define J = [0,-1;1,0], R(t) = expm(t*J) = [cos(t),-sin(t);sin(t),cos(t)]
% F(t) = R(t)*expm(t*B)*F0
% F'(t) = A(t)*F(t) with A(t) = J+R(t)*B*R(-t)

% Plot the exact solution for F(t).
J = [0,-1;1,0];
R = @(t) [cos(t),-sin(t);sin(t),cos(t)];
F0 = [1;0];
B = [1,1;1,1];
t = linspace(0,1,11);
Nt = length(t);
Ft = zeros(2,Nt);
for j = 1:Nt
    Ft(:,j) = R(t(j))*expm(t(j)*B)*F0;
end
figure; plot(t,Ft)
xlabel('t')
legend('F1','F2')
title('Variable-A test case: exact solution')

% Plot the cumulative integration error at t = [.1,.2,...1] with
% integration step .1.
A = @(t) J+R(t)*B*R(-t);
for ord = 2:2:10
    [~,Ft_] = explde(A,t,ord,F0);
    Ft_ = [Ft_{:}];
    figure, plot(t,Ft_-Ft)
    xlabel('t')
    ylabel('error')
    legend('F1','F2')
    title(['Cumulative integration error with step = 0.1, ord = ' ...
        num2str(ord)])
end

% Plot the single-step integration error, normalized to t^(ord+1), with
% integration steps t = [.1,.2,...1].
t(1) = [];
Ft(:,1) = [];
Nt = length(t);
t_ = shiftdim([zeros(1,Nt);t],-1); % size-[1,2,Nt]
% Modify R and A to work with paged t (e.g., size-[1,1,Nt]).
R = @(t) cell2mat({cos(t),-sin(t);sin(t),cos(t)});
A = @(t) J+pagemtimes(pagemtimes(R(t),B),R(-t)); % size-[2,2,Nt]
for ord = 2:2:10
    Ft_ = explde(A,t_,ord,F0); % size-[2,1,Nt]
    Ft_ = reshape(Ft_,2,[]); % size-[2,Nt]
    figure, plot(t,(Ft_-Ft)./t.^(ord+1),'o-')
    xlabel('t')
    ylabel(['error/t\^' num2str(ord+1)])
    legend('F1','F2')
    title(['1-step (integration error)/t\^' num2str(ord+1) ...
        ' with step = t, ord = ' num2str(ord)])
end

