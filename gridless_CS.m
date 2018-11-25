%% input parameters - f1, f2, d1, d2, frac, N
% frac - fraction for setting max and min saturation thresholds
% d1, d2 - elements of the sparse vector to multiplied to the dictionary
% f1, f2 - frequencies for generating the dictionary 
% N - number of time samples
function [y, z, recon, nc_frac] = gridless_CS(f1,f2,d1,d2,frac, N)
x1 = sine_generator(N,f1);
x2 = sine_generator(N,f2);

y = d1*x1 + d2*x2; max_y = max(y); min_y = min(y); y_org = y;
Phi = dctmtx(N)'; Psi = eye(N); 
non_clipped  = find((y <=frac*max_y) .* (y >=frac*min_y));
clipped_max = find(y >= frac*max_y); clipped_min = find(y <= frac*min_y);
clipped = find( (y >frac*max_y) | (y < frac*min_y));
Psi_clip = Psi; Psi_clip(clipped,:) = [];
A = Psi_clip*Phi;
Psi_max = Psi; Psi_max(clipped_max,:) = []; 
Psi_min = Psi; Psi_min(clipped_min, :) = [];

y_clip = y;
y_clip(y_clip > frac*max_y) = frac*max_y;
y_clip(y_clip < frac*min_y) = frac*min_y;

cvx_begin SDP quiet
variable z(N,1);
minimize norm(z,1)
subject tos
A*z == y(non_clipped)
% for i = 1:size(clipped_max)
% Psi(clipped_max(i),:)*Phi*z >= frac*max_y;
% end 
%  
% for i = 1:size(clipped_min)
% Psi(clipped_min(i),:)*Phi*z <= frac*min_y;
% end

cvx_end;
z = Phi*z;

recon = mse(y_org - z)/mse(y_org);
plot(y_org);
hold on; 
plot(z);
nc_frac = length(non_clipped)/N;
disp(nc_frac);
disp(recon);
end