%% Speech de-clipping
% Estimate speech samples lost due to clipping.
% The lost data is estimated by least squares.

%% Initialize
clc
clear
close all

%% Load Sound Data
[y, Fs] = audioread('klaxon.wav');
s1 = [1, 0.1 * Fs]; 
clear y Fs
% Read only part with sound
[y, Fs] = audioread('klaxon.wav', s1); 
% Test playing
%sound(y, Fs)

% Prepare for later processing
N = length(y); n = 1:N; 

[orig Fs] = audioread('klaxon.wav', s1);
origl = orig(:,1); origr = orig(:,2);

%%
% Original sound waveform
%{
figure(1)
clf
subplot(2,1,1); hold on; plot(origl, 'blue'); ylim([-0.15 0.15]);
title('Original speech waveform - Left');
subplot(2,1,2); hold on; plot(origr, 'blue'); ylim([-0.15 0.15]);
title('Original speech waveform - Right');
%}
%%
% Seperate inputs into left and right
yl = y(:,1); yr = y(:,2);
% Find indices of values over or under certain threshold (Set to 0.04)
upidxl = find(yl>0.04); loidxl = find(yl<-0.04);
upidxr = find(yr>0.04); loidxr = find(yr<-0.04);
yl(upidxl) = NaN; yl(loidxl) = NaN;
yr(upidxr) = NaN; yr(loidxr) = NaN;
clipped = horzcat(yl, yr);
% Test playing clipped sound
%sound(clipped, Fs);

N = length(y); n = 1:N;

%% Dispaly data
% The NaN's(clipped parts) appeared as discontinuances in the plot
%{
figure(2)
clf
subplot(2,1,1); hold on; plot(yl, 'black'); ylim([-0.04 0.04]);
title('Clipped speech waveform - Left');
subplot(2,1,2); hold on; plot(yr, 'black'); ylim([-0.04 0.04]);
title('Clipped speech waveform - Right');
%}
%% Define matrix D
% D represents the third-order derivitive
% (3rd - order difference).
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
Al = D * Scl.'; Ar = D * Scr.'; 
% QR factorization using Gram-Schmidt algorithm
tic;
[Ql, Rl] = get_inverse_via_GS_QR(Al); [Qr, Rr] = get_inverse_via_GS_QR(Ar); 
% Compose bl and br for convenience in back substitution for getting solution
bl = D * Sl.'* yl(kl); br = D * Sr.'* yr(kr); 
% Get solution (Retrieved samples) via back substitution
vl = -back_substitution(Rl, Ql.'*bl); vr = -back_substitution(Rr, Qr.'*br);
solve = toc;
fprintf('Time passed for solving: %.3f sec \n', solve);

%% Fill in unknown values
% Place the estimated samples into the signal.

xl = zeros(N,1); xl(kl) = yl(kl); xl(~kl) = vl; 
xr = zeros(N,1); xr(kr) = yr(kr); xr(~kr) = vr;
declipped = horzcat(xl, xr);
% Show the retrieved samples
%{
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
%}

%%
%sound(declipped, Fs)

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
lam = 1;
F = speye(N) + lam * D' * D;

% Solving through QR and Back Substitution
tic;
[Q, R] = get_inverse_via_GS_QR(F);                         
zl = back_substitution(R, Q.'*xl); zr = back_substitution(R, Q.'*xr);
solve = toc;

z = horzcat(zl, zr);
fprintf('Time passed for solving: %.3f sec \n', solve);

%% Play smoothed sound and save as a file
%{
figure(5)
clf
plot(n, z); title('Estimated Signal after Smoothing');

%%
%sound(z, Fs)

% Save the estimated signal as file
audiowrite('declipped_and_smoothed_LS.wav', z, Fs);
print -dpdf declipping_figure_LS
%}
