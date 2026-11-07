function esc_simulation
% ESC (Electronic Stability Control) - Control Systems + DSP simulation
% Alan A Samkutty, S5 ECA - Topic 12
%
% Chain: yaw-rate sensor (noise + vibration) -> 200 Hz sampling ->
%        8 Hz Butterworth low-pass (DSP) -> PI yaw-rate controller ->
%        corrective yaw moment Mz -> single-track vehicle (plant)
%
% Cases compared: no ESC, ESC with RAW sensor feedback, ESC with FILTERED feedback.
% Run:  esc_simulation      (MATLAB R2016b or later; needs Signal Processing Toolbox for butter)

clc; close all;

%% ---------------- Vehicle and simulation parameters ----------------
p.m  = 1500;        % mass (kg)
p.Iz = 2500;        % yaw inertia (kg m^2)
p.a  = 1.2;         % CG to front axle (m)
p.b  = 1.5;         % CG to rear axle (m)
p.g  = 9.81;
p.L  = p.a + p.b;
p.mu = 0.45;        % road friction (slippery road)
p.Fzf = p.m*p.g*p.b/p.L;      % front axle load (N)
p.Fzr = p.m*p.g*p.a/p.L;      % rear axle load (N)
p.Cf = 60000;       % front cornering stiffness (N/rad)
p.Cr = 60000;       % rear cornering stiffness (N/rad)
p.vx = 25;          % constant forward speed (m/s)
p.Kus = p.m/p.L*(p.b/p.Cf - p.a/p.Cr);   % understeer gradient

p.dt   = 1e-3;      % plant integration step (s)
p.T    = 8;         % simulation time (s)
p.fs   = 200;       % sensor / controller sampling rate (Hz)
p.step = round(1/(p.fs*p.dt));
p.fc   = 8;         % low-pass cut-off (Hz)
p.Kp   = 9000;      % PI gains
p.Ki   = 14000;
p.MzMax = 6000;     % yaw-moment limit (N m)

%% ---------------- Run the three cases ----------------
R.off  = run_case('off',  p);
R.raw  = run_case('raw',  p);
R.filt = run_case('filt', p);

%% ---------------- Performance numbers ----------------
rmsf = @(x) sqrt(mean(x.^2));
names = {'off','raw','filt'};
labels = {'No ESC','ESC raw feedback','ESC filtered feedback'};
fprintf('\n%-24s %-18s %-16s %-18s\n','Case','RMS yaw err (rad/s)','Peak beta (deg)','Mz chatter (RMS step)');
for k = 1:3
    c = R.(names{k});
    fprintf('%-24s %-18.4f %-16.2f %-18.1f\n', labels{k}, rmsf(c.rdes-c.r), max(abs(c.beta)), rmsf(diff(c.Mz)));
end

%% ---------------- Figure 1: DSP stage (time + FFT) ----------------
f = R.filt; t = f.t;
figure('Name','DSP stage','Position',[100 100 1100 380]);
subplot(1,2,1);
idx = t>1 & t<3;
plot(t(idx), f.rmeas(idx), 'Color',[0.7 0.7 0.7]); hold on;
plot(t(idx), f.rfb(idx),   'Color',[0.75 0.22 0.17], 'LineWidth',1.6);
xlabel('Time (s)'); ylabel('Yaw rate (rad/s)');
legend('Raw yaw-rate sensor','After 8 Hz Butterworth LPF','Location','best');
title('Sensor signal before / after filtering'); grid on;

subplot(1,2,2);
[fr, X1] = spec(f.rmeas, p.fs);
[~,  X2] = spec(f.rfb,   p.fs);
semilogy(fr, X1, 'Color',[0.6 0.6 0.6]); hold on;
semilogy(fr, X2, 'Color',[0.75 0.22 0.17]);
xline(p.fc,'--k',sprintf('f_c = %d Hz',p.fc));
xlim([0 100]); xlabel('Frequency (Hz)'); ylabel('Magnitude');
legend('Raw','Filtered'); title('FFT spectrum (f_s = 200 Hz)'); grid on;

%% ---------------- Figure 2: filter frequency response ----------------
[bf, af] = butter(2, p.fc/(p.fs/2));
figure('Name','Filter response','Position',[100 100 700 350]);
freqz(bf, af, 1024, p.fs); title('2nd-order Butterworth LPF (8 Hz) - magnitude and phase (lag)');

%% ---------------- Figure 3: closed-loop results ----------------
figure('Name','Closed loop','Position',[100 100 1300 380]);
subplot(1,3,1);
plot(t, R.off.rdes,'k--'); hold on;
plot(t, R.off.r,'Color',[0.75 0.22 0.17]);
plot(t, R.filt.r,'Color',[0.10 0.50 0.22]);
xlabel('Time (s)'); title('Yaw rate (rad/s)');
legend('Desired','No ESC','ESC (filtered feedback)','Location','best'); grid on;

subplot(1,3,2);
plot(t, R.off.beta,'Color',[0.75 0.22 0.17]); hold on;
plot(t, R.filt.beta,'Color',[0.10 0.50 0.22]);
xlabel('Time (s)'); title('Body sideslip angle (deg)');
legend('No ESC','ESC','Location','best'); grid on;

subplot(1,3,3);
plot(t, R.raw.Mz,'Color',[0.90 0.50 0.13]); hold on;
plot(t, R.filt.Mz,'Color',[0.10 0.37 0.71],'LineWidth',1.2);
xlabel('Time (s)'); title('Corrective yaw moment M_z (N m)');
legend('Raw feedback','Filtered feedback','Location','best'); grid on;
end

%% =====================================================================
function out = run_case(mode, p)
% Simulates one case. mode: 'off' | 'raw' | 'filt'
rng(3);                                   % same noise sequence for every case
N  = round(p.T/p.dt);
nrec = ceil(N/5);
rec.t=zeros(1,nrec); rec.r=rec.t; rec.rdes=rec.t; rec.rmeas=rec.t;
rec.rfb=rec.t; rec.Mz=rec.t; rec.beta=rec.t; rec.d=rec.t;

[bf, af] = butter(2, p.fc/(p.fs/2));      % DSP: low-pass filter
zf = zeros(max(length(af),length(bf))-1, 1);

tire = @(alpha, C, Fz) p.mu*Fz*tanh(C*alpha/(p.mu*Fz));   % saturating tyre force

vy = 0; r = 0; Mz = 0; I = 0; rfb = 0; rm = 0;
Ts = 1/p.fs; k = 0;
rlim = 0.85*p.mu*p.g/p.vx;

for i = 0:N-1
    t = i*p.dt;

    % driver steering input: 0.5 Hz sine, 2.6 deg amplitude, 0.5 s to 6 s
    if t > 0.5 && t < 6
        d = 0.045*sin(2*pi*0.5*t)*min(1, t/0.3);
    else
        d = 0;
    end

    % tyre slip angles and lateral forces (single-track model)
    af_ = d - (vy + p.a*r)/p.vx;
    ar_ =   - (vy - p.b*r)/p.vx;
    Ff = tire(af_, p.Cf, p.Fzf);
    Fr = tire(ar_, p.Cr, p.Fzr);

    % ---------- sampled sensor + DSP + controller (every Ts) ----------
    if mod(i, p.step) == 0
        % sensor: true yaw rate + white noise + 45 Hz vibration
        rm = r + 0.012*randn + 0.015*sin(2*pi*45*t);

        % DSP stage
        if strcmp(mode,'filt')
            [rfb, zf] = filter(bf, af, rm, zf);   % one-sample IIR update
        else
            rfb = rm;
        end

        % reference yaw rate from driver input, limited by friction
        rdes = p.vx*d/(p.L + p.Kus*p.vx^2);
        rdes = min(max(rdes, -rlim), rlim);

        % PI controller with anti-windup (integrator clamp) and saturation
        e = rdes - rfb;
        if strcmp(mode,'off')
            Mz = 0;
        else
            I  = min(max(I + e*Ts, -0.3), 0.3);
            Mz = min(max(p.Kp*e + p.Ki*I, -p.MzMax), p.MzMax);
        end
    end

    % ---------- plant: yaw dynamics ----------
    vy_dot = (Ff*cos(d) + Fr)/p.m - p.vx*r;
    r_dot  = (p.a*Ff*cos(d) - p.b*Fr + Mz)/p.Iz;
    vy = vy + vy_dot*p.dt;
    r  = r  + r_dot*p.dt;

    % ---------- logging (every 5 ms) ----------
    if mod(i,5) == 0
        k = k + 1;
        rd = p.vx*d/(p.L + p.Kus*p.vx^2);
        rec.t(k)=t; rec.r(k)=r; rec.rdes(k)=min(max(rd,-rlim),rlim);
        rec.rmeas(k)=rm; rec.rfb(k)=rfb; rec.Mz(k)=Mz;
        rec.beta(k)=atan2d(vy, p.vx); rec.d(k)=d;
    end
end
out = structfun(@(x) x(1:k), rec, 'UniformOutput', false);
end

%% =====================================================================
function [fr, X] = spec(x, fs)
% Single-sided amplitude spectrum with Hann window
x = x(:) - mean(x);
w = hann(length(x));
X = abs(fft(x.*w))/sum(w)*2;
n = length(x);
X = X(1:floor(n/2)+1);
fr = (0:floor(n/2))'*fs/n;
end
