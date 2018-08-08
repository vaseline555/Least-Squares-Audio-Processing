%% Missing Sample estiamtion using least square
% Estimate speech samples lost using least squares.

%% Initialize
clc
clear
close all

%% Load Sound Data
% Read sound
[y, Fs] = audioread('chu.wav'); 
% Test playing
sound(y, Fs) 
% Prepare for later processing
N = length(y); n = 1:N; 

%% Make a missing sample
% Seperate sample into left(right) input : sl(sr)
sl = y(:,1); sr = y(:,2); 
lenl = size(sl, 1); lenr = size(sr, 1); 

% Generate random indices to make corresponding elements to NaN (missing sample)
indl = randperm(lenl, int64(lenl * 0.5)); indr = randperm(lenr, int64(lenr * 0.5));
% Randomly eliminate 50% of the sample
sl(indl) = NaN; sr(indr) = NaN; 
% Merge the processed missing sample
missing = horzcat(sl, sr);
% Test playing the missing sample
sound(missing, Fs) 

%% Dispaly data
% The NaN's(missed parts) appear as gaps in the plot
figure(1)
clf
subplot(221); hold on; plot(n, y(:,1), 'black'); title('Original Data - Left');
subplot(222); hold on; plot(n, y(:,2), 'black'); title('Original Data - Right');

subplot(223); hold on; plot(n, missing(:,1), 'black'); title('Missing Data - Left');
subplot(224); hold on; plot(n, missing(:,2), 'black'); title('Missing Data - Right');

%% Define matrix D
% D represents the third-order derivative
% (2nd - order difference).
e = ones(N, 1);
D = spdiags([e, -2 * e, e], 0:2, N-2, N);

%% Define matrices S and Sc
% kl,kr : logical vectors (0 if an element is NaN)
kl = isfinite(sl); kr = isfinite(sr); 
% Sl, Sr : sampling matrix
Sl = speye(N); Sr = speye(N);
Sl(~kl, :) = []; Sr(~kr, :) = []; 
% Scl, Scr : complement of Scl, Scr
Scl = speye(N); Scr = speye(N); 
Scl(kl, :) = []; Scr(kr, :) = [];
% Ll, Lr : number of missing values
Ll = sum(~kl); Lr = sum(~kr); 

%% Estimate missing data
vl = -(Scl * (D' * D) * Scl') \ (Scl * D' * D * Sl' * sl(kl));
vr = -(Scr * (D' * D) * Scr') \ (Scr * D' * D * Sr' * sr(kr));

%% Fill in unknown values
% non-NaN values are from original matrix ; NaNs are from estimated matrix
xl = zeros(N,1); xl(kl) = sl(kl); xl(~kl) = vl; 
xr = zeros(N,1); xr(kr) = sr(kr); xr(~kr) = vr;

% Show the retrieved samples
figure(2)
clf
subplot(211); plot(n, sl, 'k', n(~kl), xl(~kl) ,'k.'); legend('Original data', 'Estiamted data');
subplot(212); plot(n, sr, 'k', n(~kr), xr(~kr) ,'k.'); legend('Origianl data', 'Estiamted data');

%% Show and play finally estimated signal
figure(3)
clf
estimated = horzcat(xl, xr);
subplot(321); hold on; plot(n, y(:,1)); title('Original Data - Left');
subplot(322); hold on; plot(n, y(:,2)); title('Original Data - Right');

subplot(323); hold on; plot(n, xl, 'r'); title('Estimated Data - Left');
subplot(324); hold on; plot(n, xr, 'r'); title('Estimated Data - Right');

subplot(3,2,[5,6]); hold on; plot(n, estimated); title('Estimated Signal');
sound(estimated, Fs)
print -dpdf estimated_figure_solver

% Save the estimated signal as file
audiowrite('estimated_solver.wav',estimated, Fs);

%% Calculate RMSE
% Left Signal
RMSE_L = sqrt(mean((xl - y(:,1)).^2)) * 100;
% Right Signal
RMSE_R = sqrt(mean((xr - y(:,2)).^2)) * 100;
% Print RMSE in percent
sprintf('RMSE of the retrieved left signal is %.4f percent', RMSE_L)
sprintf('RMSE of the retrieved right signal is %.4f percent', RMSE_R)
