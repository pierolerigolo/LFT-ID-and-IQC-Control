close all;clc;clear
%% open the experimental sinesweep
% Disturbance at the output:  y = G u + d,   u = -C y
%   open-loop sweep  (u = 0)        : y_ol = d
%   closed-loop sweep (same d)      : y_cl = G u_cl + d
% so, with the recorded command u_cl:
%   C_meas = -U_cl / Y_cl                     controller that actually ran
%   G_cl   = (Y_cl - Y_ol) / U_cl             plant seen by the controller
% and without the command, from an identified plant G_id:
%   C_rec  = (Y_ol / Y_cl - 1) / G_id  = (1/S - 1) / G_id
expforlder=fullfile('data','Experiments/','Old Experiments/','controller_sus_7nm_order18_82hz_-34db_nperfhinfnperfoptiIQC_eta9.139390e-01_itr26/');

[S, Uclfrd, Yclfrd, Yolfrd] = utils.sweepS(fullfile(expforlder,'gfrf0_sswep.mat'));

load(fullfile('data','GFRF_iddata_0729_spa.mat'))      % GFRF: one channel per mass (realmass)
Ksave = load(fullfile('data','finalrepportcont', ...
    'controller_sus_7nm_order18_82hz_-34db_nperfhinfnperfoptiIQC_eta9.139390e-01.mat'), 'K_itr');
K   = Ksave.K_itr{26};
iOP = 1;                                               % GFRF(iOP) used for C_rec (1 = 0 kg)

%% Spectra on a common frequency grid
Ts  = S.Ts;
w   = S.Frequency;
sel = w >= max(2*pi*15, min(GFRF.Frequency)) & w <= min(2*pi*450, max(GFRF.Frequency));
w   = w(sel);   f = w/2/pi;
U   = squeeze(Uclfrd.ResponseData);  U  = U(sel);
Yc  = squeeze(Yclfrd.ResponseData);  Yc = Yc(sel);
Yo  = squeeze(Yolfrd.ResponseData);  Yo = Yo(sel);

nb  = 61;                                              % bins averaged (~1 Hz) in each ratio
xr  = @(N, D) movsum(N.*conj(D), nb)./movsum(abs(D).^2, nb);

C_meas = -xr(U, Yc);                                   % controller from the recorded command
G_cl   =  xr(Yc - Yo, U);                              % plant seen in closed loop
Smeas  =  xr(Yc, Yo);                                  % S = y_cl / y_ol
G_id   = squeeze(interp(GFRF(iOP), w).ResponseData);
C_rec  = (1./Smeas - 1)./G_id;                         % controller recovered with GFRF(iOP)
Kw     = squeeze(freqresp(K, w));

band = f >= 20 & f <= 200 & abs(f - 80) > 3;           % away from the controller peak
fprintf('C_meas / K : median %+.2f dB, %+.1f deg\n', median(mag2db(abs(C_meas(band)./Kw(band)))), ...
    median(rad2deg(angle(C_meas(band)./Kw(band)))));
fprintf('C_rec  / K : median %+.2f dB, %+.1f deg (with GFRF(%d), %.2f kg)\n', median(mag2db(abs(C_rec(band)./Kw(band)))), ...
    median(rad2deg(angle(C_rec(band)./Kw(band)))), iOP, realmass(iOP));
pk = f > 25 & f < 75;
[~, k1] = max(abs(G_cl(pk)));   fk = f(pk);
fprintf('plant seen in closed loop: resonance at %.1f Hz\n', fk(k1));
for i = 1:size(GFRF, 2)
    Gi = squeeze(interp(GFRF(i), w).ResponseData);
    [~, k2] = max(abs(Gi(pk)));
    r  = G_cl(band)./Gi(band);
    fprintf('  GFRF(%d) %5.2f kg: resonance %.1f Hz, G_cl/GFRF median %+.1f dB (spread %.1f dB), %+.0f deg (spread %.0f deg)\n', ...
        i, realmass(i), fk(k2), median(mag2db(abs(r))), std(mag2db(abs(r))), median(rad2deg(angle(r))), std(rad2deg(angle(r))));
end

%% Controller: designed, measured from the command, recovered from GFRF(iOP)
figure('Name', 'controller recovery')
subplot(2, 1, 1)
semilogx(f, mag2db(abs(Kw)), 'k', f, mag2db(abs(0.1*Kw)), 'k:', ...
    f, mag2db(abs(C_meas)), 'b', f, mag2db(abs(C_rec)), 'r'); grid on
ylabel('|C| (dB)'); xlim([15 450])
legend('designed K', '0.1 K', 'C_{meas} = -U_{cl}/Y_{cl}', sprintf('C_{rec} with GFRF(%d)', iOP), 'Location', 'best')
subplot(2, 1, 2)
semilogx(f, rad2deg(angle(Kw)), 'k', f, rad2deg(angle(C_meas)), 'b', f, rad2deg(angle(C_rec)), 'r'); grid on
ylabel('phase (deg)'); xlabel('Frequency (Hz)'); xlim([15 450])

%% Plant: seen in closed loop vs identified (GFRF)
figure('Name', 'plant seen in closed loop')
col = lines(size(GFRF, 2));
subplot(2, 1, 1); hold on
semilogx(f, mag2db(abs(G_cl)), 'k', 'LineWidth', 1.5, 'DisplayName', 'G_{cl} = (Y_{cl} - Y_{ol})/U_{cl}')
for i = 1:size(GFRF, 2)
    Gi = squeeze(interp(GFRF(i), w).ResponseData);
    semilogx(f, mag2db(abs(Gi)), 'Color', col(i, :), 'DisplayName', sprintf('GFRF(%d) %.2f kg', i, realmass(i)))
end
set(gca, 'XScale', 'log'); grid on; ylabel('|G| (dB)'); xlim([15 450]); legend('Location', 'best')
subplot(2, 1, 2); hold on
semilogx(f, rad2deg(unwrap(angle(G_cl))), 'k', 'LineWidth', 1.5)
for i = 1:size(GFRF, 2)
    Gi = squeeze(interp(GFRF(i), w).ResponseData);
    semilogx(f, rad2deg(unwrap(angle(Gi))), 'Color', col(i, :))
end
set(gca, 'XScale', 'log'); grid on; ylabel('phase (deg)'); xlabel('Frequency (Hz)'); xlim([15 450])

%% Sensitivity: measured vs theory with K, 0.1 K and the plant seen in closed loop
figure('Name', 'sensitivity')
semilogx(f, mag2db(abs(Smeas)), 'k', ...
    f, mag2db(abs(1./(1 + Kw.*G_id))), 'b--', ...
    f, mag2db(abs(1./(1 + 0.1*Kw.*G_id))), 'r--', ...
    f, mag2db(abs(1./(1 + 0.1*Kw.*G_cl))), 'g--'); grid on
legend('measured y_{cl}/y_{ol}', sprintf('1/(1+K GFRF(%d))', iOP), sprintf('1/(1+0.1K GFRF(%d))', iOP), ...
    '1/(1+0.1K G_{cl})', 'Location', 'best')
ylabel('|S| (dB)'); xlabel('Frequency (Hz)'); xlim([15 450])

%% Theoretical command: K applied to the recorded acceleration (lsim) vs recorded command
[y, ~, u] = utils.loadexp(fullfile(expforlder, 'gfrf0_sswep.mat'));
k0   = find(abs(u) > 1e-3*max(abs(u)), 1);                 % controller switched on
idx  = (k0:numel(y)).';
t    = (idx - 1)*Ts;
u_th = lsim(K, -y(idx), t - t(1));                          % u = -K y from the recorded y
ur   = u(idx);

% best gain for lags of -3..3 samples (a delay would move the best fit off lag 0)
lags = -3:3;  gl = zeros(size(lags));  res = gl;
for j = 1:numel(lags)
    a = ur(max(1, 1 + lags(j)):end + min(0, lags(j)));
    b = u_th(max(1, 1 - lags(j)):end + min(0, -lags(j)));
    gl(j)  = (b.'*a)/(b.'*b);
    res(j) = norm(a - gl(j)*b)/norm(a);
end
disp(table(lags.', gl.', mag2db(abs(gl.')), 100*res.', ...
    'VariableNames', {'lag_samples', 'gain', 'gain_dB', 'residual_pct'}))
[~, jb] = min(res);
fprintf('recorded command = %.4f x lsim(K) with lag %+d samples (residual %.2f %%)\n', gl(jb), lags(jb), 100*res(jb));

figure('Name', 'command: recorded vs lsim(K)')
[~, kp] = max(abs(u_th));                                   % zoom where the command is largest
zm = t(kp) + [-0.05 0.05];
subplot(2, 1, 1)
plot(t, ur, 'b', t, u_th, 'k:', t, gl(jb)*u_th, 'r--'); grid on
legend('recorded command', 'lsim(K, -y)', sprintf('%.3f x lsim(K, -y)', gl(jb)), 'Location', 'best')
ylabel('controller output'); title('closed-loop part of the run')
subplot(2, 1, 2)
plot(t, ur, 'b', t, gl(jb)*u_th, 'r--'); grid on; xlim(zm)
legend('recorded command', sprintf('%.3f x lsim(K, -y)', gl(jb)), 'Location', 'best')
xlabel('time (s)'); ylabel('controller output'); title('zoom')
