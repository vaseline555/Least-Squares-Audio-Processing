%% Speech de-clipping
% Estimate speech samples lost due to clipping.
% The lost data is estimated by least squares.

%% Initialize
clc
clear
close all

%% Load Sound Data
[y Fs] = audioread('ah.wav');
sound(y, Fs)

[orig Fs] = audioread('ah.wav');
origl = orig(:,1); origr = orig(:,2);

%%
% Original sound waveform
figure(1)
clf
subplot(2,1,1); hold on; plot(origl, 'black'); ylim([-0.7 0.7]);
title('Original speech waveform - Left');
subplot(2,1,2); hold on; plot(origr, 'black'); ylim([-0.7 0.7]);
title('Original speech waveform - Right');

%%
% Seperate inputs into left and right
yl = y(:,1); yr = y(:,2);
% Find indices of values over or under certain threshold (Set to 0.3)
upidxl = find(yl>0.3); loidxl = find(yl<-0.3);
upidxr = find(yr>0.3); loidxr = find(yr<-0.3);
yl(upidxl) = NaN; yl(loidxl) = NaN;
yr(upidxr) = NaN; yr(loidxr) = NaN;
clipped = horzcat(yl, yr);
% Test playing clipped sound
sound(clipped, Fs);

N = length(y); n = 1:N;

%% Dispaly data
% The NaN's(clipped parts) appeared as discontinuances in the plot
figure(2)
clf
subplot(2,1,1); hold on; plot(yl, 'black'); ylim([-0.33 0.33]);
title('Clipped speech waveform - Left');
subplot(2,1,2); hold on; plot(yr, 'black'); ylim([-0.33 0.33]);
title('Clipped speech waveform - Right');

%% Define matrix D
% D represents the third-order derivitive
% (3rd-order difference).
e = ones(N, 1);
D = spdiags([e -3*e 3*e -e], 0:3, N-3, N);

%% Define matrices S and Sc
% kl,kr : logical vectors (0 if an element is NaN)
kl = isfinite(yl); kr = isfinite(yr); 
% Sl, Sr : sampling matrix
Sl = speye(N); Sr = speye(N);
Sl(~kl, :) = []; Sr(~kr, :) = []; 
% Scl, Scr : complement of Scl, Scr
Scl = speye(N); Scr = speye(N); 
Scl(kl, :) = []; Scr(kr, :) = [];
% Ll, Lr : number of missing values
Ll = sum(~kl); Lr = sum(~kr); 

%% Estimate missing data
% Fast solver is...
vl = -(Scl * (D' * D) * Scl') \ (Scl * D' * D * Sl' * yl(kl));
vr = -(Scr * (D' * D) * Scr') \ (Scr * D' * D * Sr' * yr(kr));

%% Fill in unknown values
% Place the estimated samples into the signal.
xl = zeros(N,1); xl(kl) = yl(kl); xl(~kl) = vl; 
xr = zeros(N,1); xr(kr) = yr(kr); xr(~kr) = vr;
declipped = horzcat(xl, xr);

% Show the retrieved samples   
figure(3)
clf
subplot(211); hold on; plot(n, yl, 'k', n(~kl), xl(~kl) ,'b.');
legend('Known data', 'Estiamted data'); title('Estimated values - Left');
subplot(212); hold on;plot(n, yr, 'k', n(~kr), xr(~kr) ,'b.');
legend('Known data', 'Estiamted data'); title('Estimated values - Right');

figure(4)
subplot(211); hold on; plot(n, xl, 'red', n, yl, 'black', 'linewidth', 2);
title('Estimated Signal - left')
legend('Filled in', 'Clipped data');
subplot(212); hold on; plot(n, xr, 'red', n, yr, 'black', 'linewidth', 2);
title('Estimated Signal - right')
legend('Filled in', 'Clipped data');

sound(declipped, Fs)

%% Calculate RMSE
% Left Signal
RMSE_L = sqrt(mean((xl - origl).^2)) * 100;
% Right Signal
RMSE_R = sqrt(mean((xr - origr).^2)) * 100;
% Print RMSE in percent
sprintf('RMSE of the declipped left signal is %.4f percent', RMSE_L)
sprintf('RMSE of the declipped right signal is %.4f percent', RMSE_R)

%% Smoothing and Plot
% Define 2nd-order Difference matrix
e = ones(N, 1);
D = spdiags([e -2*e e], 0:2, N-2, N);

% Constraint value lambda
lam = 20;
F = speye(N) + lam * D' * D;

% Fast solver is...
z = F \ declipped;

%% Play smoothed sound
figure(5)
clf
plot(n, z); title('Estimated Signal after Smoothing');
print -dpdf declipping_figure_solver
sound(z, Fs)

% Save the estimated signal as file
audiowrite('declipped_and_smoothed_solver.wav', z, Fs);

