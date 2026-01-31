% Parameters
c = 3e8; % Speed of light
fc = 76.5e9; % Carrier frequency = 76.5 GHz
Tch = 2.1e-6;
PRI = 8.4e-6;
PRF = 1/PRI;
fmin = -0.5e9;
fmax = 0.5e9;
BW = fmax - fmin;
s = (fmax - fmin)/Tch;
N_pulses = 1880; 
Fs = 2e9; 
dt = 1/Fs; 
t = 0:dt:PRI-dt;


%% Targets ------> [ positive = Going away -----------  negative = Approaching ]
R_true = [50, 150];    
V_true = [-50, 80];     
A_true = [0.8, 0.5];  % Amplitude of targets    
L = length(R_true); 

%% Transmitted Pulse
P = exp(1j * 2*pi*( fmin*t + 0.5*s*t.^2 )) .* (t<=Tch);
Tx = repmat(P, N_pulses, 1); 


%% Received Signal 
Rx = zeros(N_pulses, length(t));
SNR_dB = 15; 
noise_power = 10^(-SNR_dB/10);

fprintf('Generating Signals for %d Fixed Targets...\n', L);

for n = 1:N_pulses
    pulse_echo = zeros(size(t));
    for k = 1:L
        RTT = 2*R_true(k)/c; 
        Ns = round(RTT/dt);
        if Ns+1 <= length(t)
            fd = 2*V_true(k)*fc/c;
            doppler_phase = exp(1j*2*pi*fd*n*PRI);
            shifted = [zeros(1,Ns), P(1:end-Ns)] * A_true(k);
            pulse_echo = pulse_echo + shifted*doppler_phase;
        end
    end
    Rx(n,:) = pulse_echo + sqrt(noise_power/2)*(randn(size(t)) + 1j*randn(size(t))); 
end

% Signal Processing & Detection
timer_start = tic;
% Mixing 
Mix = Tx .* conj(Rx); % The Most Important 
% Range FFT 
N_samples_active = round(Tch * Fs); 
Mix_active = Mix(:, 1:N_samples_active);
n_fft_range = 2^nextpow2(N_samples_active) ; 
range_spectrum_matrix = fft(Mix_active, n_fft_range, 2);
 
% Doppler FFT
RDM = fftshift(fft(range_spectrum_matrix, N_pulses, 1), 1);
RDM_dB = 10*log10(abs(RDM));
processing_time = toc(timer_start);
% Axis Generation
freq_range = (0:n_fft_range-1) * (Fs / n_fft_range);
rng_axis = (freq_range * c) / (2 * s);
lambda = c / fc;
freq_doppler = -(-N_pulses/2 : N_pulses/2 - 1) * (PRF / N_pulses);
vel_axis = (freq_doppler * lambda) / 2;
%% PEAK DETECTION 
% Detect Range Peaks
first_pulse_fft = abs(range_spectrum_matrix(1, :));
threshold = max(first_pulse_fft) * 0.1; 
range_resolution = rng_axis(2) - rng_axis(1);
min_dist_meters = 10; 
min_dist_indices = round(min_dist_meters / range_resolution);
[pks, locs] = findpeaks(first_pulse_fft, ...
                        'MinPeakHeight', threshold, ...
                        'MinPeakDistance', min_dist_indices);
detected_ranges = rng_axis(locs);
% Detect Velocity
detected_velocities = zeros(1, length(detected_ranges));
for i = 1:length(detected_ranges)
    dop_cut = abs(RDM(:, locs(i))); 
    [~, v_idx] = max(dop_cut);
    detected_velocities(i) = vel_axis(v_idx);
end


%%  Plotting Section
N_show = 3; 
Tx_show = reshape(Tx(1:N_show, :).', 1, []);
Rx_show = reshape(Rx(1:N_show, :).', 1, []);
t_show  = linspace(0, N_show*PRI, length(Tx_show));
% Frequency Sawtooth
t_one_pulse = linspace(0, PRI, 1000);
freq_one_pulse = zeros(size(t_one_pulse));
active_idx = t_one_pulse <= Tch;
freq_one_pulse(active_idx) = linspace(fmin, fmax, sum(active_idx));
t_sawtooth = linspace(0, N_show*PRI, 1000*N_show);
freq_sawtooth = repmat(freq_one_pulse, 1, N_show);


% Transmitted Signal
figure;
plot(t_show*1e6, real(Tx_show)); 
title(sprintf('2. Transmitted Signal (Real Part) - Showing %d Pulses', N_show));
ylabel('Amplitude');
xlabel('Time (\mus)');
xlim([0, N_show*PRI*1e6]);
grid on;
figure('Name', 'FMCW Sawtooth', 'Color', 'w');
plot(t_sawtooth*1e6, freq_sawtooth/1e9, 'b', 'LineWidth', 2);

title('FMCW Sawtooth Waveform (Frequency vs Time)');
ylabel('Frequency (GHz)');
xlabel('Time (\mus)');
grid on;
xlim([0, N_show*PRI*1e6]);

% Add Annotations
text(Tch*1e6/2, 0, 'Active Chirp', 'HorizontalAlignment', 'center', ...
    'BackgroundColor', 'w', 'EdgeColor', 'k');
text(PRI*1e6 - 2, 0, 'Idle Time', 'HorizontalAlignment', 'center', ...
    'BackgroundColor', 'w', 'EdgeColor', 'k');
% Received Signal
figure;
plot(t_show*1e6, real(Rx_show)); 
title(sprintf('3. Received Signal (Real Part) - Showing %d Pulses', N_show));
ylabel('Amplitude');
xlabel('Time (\mus)');
xlim([0, N_show*PRI*1e6]);
grid on;

% Figure 2: Range FFT Detection
figure('Name', 'Range Detection', 'Color', 'w');
plot(rng_axis, 20*log10(abs(range_spectrum_matrix(1,:)))); hold on;
plot(detected_ranges, 20*log10(pks), 'rv', 'MarkerFaceColor', 'r');
yline(20*log10(threshold), '--k', 'Threshold'); 
title('Range FFT');
xlabel('Range (m)'); ylabel('dB'); grid on; xlim([0, 300]);

% Figure 3: Doppler Spectrum
figure('Name', 'Doppler Spectrum', 'Color', 'w');
num_plots = length(detected_ranges);
for i = 1:num_plots
    subplot(num_plots, 1, i);
    dop_cut_dB = RDM_dB(:, locs(i));
    plot(vel_axis, dop_cut_dB, 'LineWidth', 1.5); hold on;
    [min_val, true_idx] = min(abs(R_true - detected_ranges(i)));
    if min_val < 10 
        xline(V_true(true_idx), '--g', 'LineWidth', 2);
    end
    
    title(sprintf('Target detected at %.1f m', detected_ranges(i)));
    xlabel('Velocity (m/s)'); ylabel('dB');
    grid on; xlim([-100, 100]);
end

%  FIGURE 4: Range-Doppler Map
figure('Name', 'Range-Doppler Map (Direction)', 'Color', 'w');
imagesc(rng_axis, vel_axis, RDM_dB);
axis xy; colormap('jet');
clim([max(RDM_dB(:))-45, max(RDM_dB(:))]); 

hold on;
plot(R_true, V_true, 'wo', 'MarkerSize', 14, 'LineWidth', 1.5);

for m = 1:length(detected_ranges)
    r_val = detected_ranges(m);
    v_val = detected_velocities(m);
    
    if v_val < -0.5 
        % NEGATIVE = APPROACHING
        plot(r_val, v_val, '^', 'Color', 'g', 'MarkerSize', 12, ...
            'LineWidth', 2, 'MarkerFaceColor', 'g');
        text(r_val, v_val - 15, 'Approaching', 'Color', 'g', ...
            'FontWeight', 'bold', 'FontSize', 9, 'HorizontalAlignment', 'center');
            
    elseif v_val > 0.5
        % POSITIVE = RECEDING
        plot(r_val, v_val, 'v', 'Color', 'r', 'MarkerSize', 12, ...
            'LineWidth', 2, 'MarkerFaceColor', 'r');
        text(r_val, v_val + 15, 'Receding', 'Color', 'r', ...
            'FontWeight', 'bold', 'FontSize', 9, 'HorizontalAlignment', 'center');
    end
end

title('Range-Doppler Map (Neg=Approaching, Pos=Receding)');
xlabel('Range (m)'); 
ylabel('Velocity (m/s)');
xlim([0, 300]); 
ylim([-100, 100]);
h1 = plot(NaN,NaN,'^g', 'MarkerFaceColor','g');
h2 = plot(NaN,NaN,'vr', 'MarkerFaceColor','r');
h3 = plot(NaN,NaN,'wo', 'LineWidth', 1.5);
legend([h1, h2, h3], 'Approaching (-)', 'Receding (+)', 'Ground Truth', 'Location', 'best');

%  FIGURE 5: 3D Range-Doppler Map 
figure('Name', '3D Range-Doppler', 'Color', 'w');
[R_grid, V_grid] = meshgrid(rng_axis, vel_axis);
% Plot Surface
surf(R_grid, V_grid, RDM_dB, 'EdgeColor', 'none');
colormap('jet'); colorbar;
view(45, 45); 
xlim([0, 300]); ylim([-100, 100]);
zlim([max(RDM_dB(:))-60, max(RDM_dB(:))]); 
title('3D Range-Doppler Map');
xlabel('Range (m)'); ylabel('Velocity (m/s)'); zlabel('Magnitude (dB)');


%% Printing in Command Window
fprintf('RANGE DETECTION RESULTS\n');
fprintf('Table 1 : Processing Time: %.4f seconds\n', processing_time);
fprintf('%-12s | %-12s | %-10s\n', 'True R (m)', 'Det R (m)', 'Error (m)');
[R_true_sorted, sort_idx] = sort(R_true);
V_true_sorted = V_true(sort_idx);

for k = 1:L
   [val, idx] = min(abs(detected_ranges - R_true_sorted(k)));
   if val < 10 
       det_r = detected_ranges(idx);
       err_r = det_r - R_true_sorted(k);
       fprintf('%-12.2f | %-12.2f | %-10.4f\n', R_true_sorted(k), det_r, err_r);
   else
       fprintf('%-12.2f | %-12s | %-10s\n', R_true_sorted(k), 'MISSED', '-');
   end
end

fprintf(' Table 2 : VELOCITY DETECTION RESULTS \n');
fprintf('%-12s | %-12s | %-10s | %-15s\n', 'True V', 'Det V', 'Error', 'Direction');
for k = 1:L
   [val, idx] = min(abs(detected_ranges - R_true_sorted(k)));
   if val < 10
       det_v = detected_velocities(idx);
       err_v = det_v - V_true_sorted(k);
       if det_v < -0.1
           dir_str = 'Approaching (-)'; 
       elseif det_v > 0.1
           dir_str = 'Receding (+)';
       else
           dir_str = 'Static';
       end
       
       fprintf('%-12.2f | %-12.2f | %-10.4f | %-15s\n', V_true_sorted(k), det_v, err_v, dir_str);
   else
       fprintf('%-12.2f | %-12s | %-10s | %-15s\n', V_true_sorted(k), 'MISSED', '-', '-');
   end
end