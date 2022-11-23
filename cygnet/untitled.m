syms re ri
syms A I
eq1 = A - pi*(re^2 - ri^2) == 0;
eq2 = I - A*(re^2 + ri^2)/2 == 0;
sols = solve([eq1, eq2], [re ri]);

% [vals, idx] = sort(-(2^(1/2).*((Si(:,1).^2 + 2.*pi.*Si(:,2))./Si(:,1)).^(1/2))./(2*pi^(1/2)))
% plot(Si(idx,1))
% figure()
% plot(Si(idx,2))