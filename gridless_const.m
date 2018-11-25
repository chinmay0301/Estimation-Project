%% input parameters - f1, f2, d1, d2, frac, N
% frac - fraction for setting max and min saturation thresholds
% d1, d2 - elements of the sparse vector to multiplied to the dictionary
% f1, f2 - frequencies for generating the dictionary 
% N - number of time samples
function [y, z, recon, nc_frac] = gridless_unconst(f1,f2,d1,d2,frac, N)
x1 = sine_generator(N,f1);
x2 = sine_generator(N,f2);

y = d1*x1 + d2*x2; max_y = max(y); min_y = min(y);
non_clipped  = find((y <=frac*max_y) .* (y >=frac*min_y));
clipped_max = find(y >= frac*max_y); clipped_min = find(y <= frac*min_y);
y_clip = y;
y_clip(y_clip > frac*max_y) = frac*max_y;
y_clip(y_clip < frac*min_y) = frac*min_y;

cvx_begin SDP quiet
variable z(N,1);
variable x;
variable u(N,1);
minimize x + norm(u,1) + 0.5* sum((y_clip(non_clipped) - z(non_clipped)).^2)
subject tos
[x z';z toeplitz(u)] >=0;
z(clipped_max) >= frac*max_y;
z(clipped_min) <= frac*min_y;
cvx_end;

recon = mse(y - z)/mse(y);
plot(y);
hold on; 
plot(z);
nc_frac = length(non_clipped)/N;
disp(nc_frac);
disp(recon);
end