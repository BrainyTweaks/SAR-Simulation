clear;
clc;
close all;

%% SAR Simulation
% LFM chirp -> raw SAR data -> range compression
% -> backprojection -> focused image

%% Radar parameters

c = 3e8;                    % speed of light
fc = 5.3e9;                 % carrier frequency
lambda = c/fc;              % wavelength

B = 40e6;                   % bandwidth
Tp = 5e-6;                  % pulse duration
K = B/Tp;                   % chirp rate

v = 100;                    % radar velocity
R0 = 5000;                  % reference range

fprintf("Wavelength = %.4f m\n", lambda);
fprintf("Range resolution = %.2f m\n", c/(2*B));
fprintf("Chirp rate = %.2e Hz/s\n", K);


%% Sampling

fs = 80e6;
dt = 1/fs;


%% Transmitted chirp

t_tx = -Tp/2 : dt : Tp/2-dt;
Nt = length(t_tx);

tx = exp(1j*pi*K*t_tx.^2);

% Normalize the chirp
tx = tx / sqrt(sum(abs(tx).^2));

figure;
plot(t_tx*1e6, real(tx), 'LineWidth', 1);

xlabel('Time (\mus)');
ylabel('Amplitude');
title('Transmitted LFM Chirp');
grid on;


%% Radar movement

aperture = 200;             % aperture length
Na = 401;                   % number of radar positions

x_radar = linspace(-aperture/2, aperture/2, Na);

figure;
plot(x_radar, zeros(size(x_radar)), 'o-');

xlabel('Along-track position (m)');
ylabel('Cross-track position (m)');
title('Synthetic Aperture');
grid on;
axis equal;


%% Targets
% [x position, range, reflectivity]

targets = [
    -40, 5000, 1.0;
      0, 5000, 0.7;
     40, 5050, 0.4
];

fprintf("\nTargets:\n");

for n = 1:size(targets,1)
    fprintf("Target %d: X = %.1f m, Range = %.1f m, Reflectivity = %.2f\n", ...
        n, targets(n,1), targets(n,2), targets(n,3));
end


%% Receive signal

t_rx = -12e-6 : dt : 12e-6;
Nr = length(t_rx);


%% Generate raw SAR data

raw = zeros(Nr, Na);

fprintf("\nGenerating raw SAR data...\n");

for k = 1:Na

    echo_total = zeros(1, Nr);

    % Current radar position
    radar_x = x_radar(k);

    for n = 1:size(targets,1)

        target_x = targets(n,1);
        target_y = targets(n,2);
        target_amp = targets(n,3);

        % Distance between radar and target
        R = sqrt((radar_x-target_x)^2 + target_y^2);

        % Delay relative to the reference range
        tau = 2*(R-R0)/c;

        t_shift = t_rx - tau;

        echo = zeros(1, Nr);

        % Only keep the part where the chirp exists
        valid = abs(t_shift) <= Tp/2;

        echo(valid) = exp(1j*pi*K*t_shift(valid).^2);

        % Phase change due to propagation
        propagation_phase = ...
            exp(-1j*4*pi*(R-R0)/lambda);

        echo_total = echo_total + ...
            target_amp .* echo .* propagation_phase;
    end

    raw(:,k) = echo_total(:);
end

fprintf("Raw SAR data generated.\n");


%% Add noise

noise_level = 0.02;

noise = noise_level/sqrt(2) * ...
    (randn(size(raw)) + 1j*randn(size(raw)));

raw_noisy = raw + noise;


%% Raw data plot

range_raw = R0 + c*t_rx/2;

figure;
imagesc(x_radar, range_raw, abs(raw_noisy));

xlabel('Along-track position (m)');
ylabel('Range (m)');
title('Raw SAR Data');

colorbar;
axis xy;


%% Range compression

fprintf("\nPerforming range compression...\n");

% Matched filter for the transmitted chirp
matched_filter = conj(fliplr(tx));

% FFT size for convolution
Nfft = 2^nextpow2(Nr + Nt - 1);

RAW = fft(raw_noisy, Nfft, 1);
H = fft(matched_filter, Nfft).';

% Apply matched filter
range_compressed_full = ifft(RAW .* H, Nfft, 1);


% Time axis after convolution
t_conv = t_rx(1) + t_tx(1) + ...
    (0:Nfft-1)*dt;

% Keep the useful part
valid_conv = t_conv >= -12e-6 & t_conv <= 12e-6;

compressed = range_compressed_full(valid_conv,:);

t_compressed = t_conv(valid_conv);

range_compressed = R0 + c*t_compressed/2;

fprintf("Range compression complete.\n");


%% Range-compressed data

figure;

imagesc(x_radar, range_compressed, abs(compressed));

xlabel('Along-track position (m)');
ylabel('Range (m)');
title('Range-Compressed SAR Data');

colorbar;
axis xy;


%% Image grid

x_image = linspace(-100, 100, 161);
range_image = linspace(4950, 5100, 121);

Nx = length(x_image);
Ny = length(range_image);

fprintf("\nImage size = %d x %d\n", Ny, Nx);


%% Backprojection

fprintf("\nPerforming SAR focusing...\n");

sar_image = zeros(Ny, Nx);

% Range spacing after compression
dr = range_compressed(2) - range_compressed(1);

Nrc = length(range_compressed);

% Grids for the radar positions and image range
[range_grid, radar_grid] = ndgrid(range_image, x_radar);

% Used for indexing all aperture positions at once
column_offset = (0:Na-1)*Nrc;


for ix = 1:Nx

    image_x = x_image(ix);

    % Distance from each radar position to each image pixel
    R_pixel = sqrt((radar_grid-image_x).^2 + range_grid.^2);

    % Corresponding round-trip delay
    tau_pixel = 2*(R_pixel-R0)/c;

    % Convert delay to range
    sample_range = R0 + c*tau_pixel/2;

    % Fractional sample position
    sample_index = ...
        (sample_range-range_compressed(1))/dr + 1;

    index_low = floor(sample_index);
    index_high = index_low + 1;

    fraction = sample_index-index_low;

    % Check that the interpolation indices are valid
    valid = index_low >= 1 & index_high <= Nrc;

    % Avoid invalid MATLAB indexing
    index_low(~valid) = 1;
    index_high(~valid) = 1;
    fraction(~valid) = 0;

    % Linear interpolation
    data_low = compressed(index_low + column_offset);
    data_high = compressed(index_high + column_offset);

    signal = data_low + ...
        fraction .* (data_high-data_low);

    signal(~valid) = 0;

    % Compensate the propagation phase
    phase_correction = ...
        exp(1j*4*pi*(R_pixel-R0)/lambda);

    signal = signal .* phase_correction;

    % Sum the signals from the whole aperture
    sar_image(:,ix) = sum(signal,2);

    if mod(ix,20) == 0
        fprintf("Processing image column %d / %d\n", ix, Nx);
    end
end

fprintf("SAR focusing complete.\n");


%% Normalize image

sar_abs = abs(sar_image);
sar_abs = sar_abs / max(sar_abs(:));


%% Focused SAR image

figure;

imagesc(x_image, range_image, sar_abs);

xlabel('Along-track position (m)');
ylabel('Range (m)');
title('Focused SAR Image');

colorbar;
axis xy;


%% Image in dB

sar_image_db = 20*log10(sar_abs + eps);

figure;

imagesc(x_image, range_image, sar_image_db);

xlabel('Along-track position (m)');
ylabel('Range (m)');
title('Focused SAR Image (dB)');

colorbar;
caxis([-40 0]);
axis xy;


%% Mark the actual target positions

hold on;

plot(targets(:,1), targets(:,2), ...
    'wx', 'MarkerSize', 12, 'LineWidth', 2);

for n = 1:size(targets,1)

    text(targets(n,1)+3, targets(n,2), ...
        sprintf('T%d',n), ...
        'Color','w', 'FontWeight','bold');

end

hold off;


%% Range profile

% Take the range profile through the strongest azimuth location
[~, strongest_x_index] = max(max(sar_abs,[],1));

range_profile = sar_abs(:,strongest_x_index);

range_profile_db = 20*log10(range_profile + eps);

figure;

plot(range_image, range_profile_db, 'LineWidth', 1.5);

xlabel('Range (m)');
ylabel('Magnitude (dB)');
title('SAR Range Profile');

grid on;
ylim([-40 5]);


%% Azimuth profile

% Take the azimuth profile through the strongest range location
[~, strongest_range_index] = max(max(sar_abs,[],2));

azimuth_profile = sar_abs(strongest_range_index,:);

azimuth_profile_db = 20*log10(azimuth_profile + eps);

figure;

plot(x_image, azimuth_profile_db, 'LineWidth', 1.5);

xlabel('Along-track position (m)');
ylabel('Magnitude (dB)');
title('SAR Azimuth Profile');

grid on;
ylim([-40 5]);


%% Find the three strongest peaks

fprintf("\nDetected image peaks:\n");

% Sort all pixels from strongest to weakest
[~, sorted_indices] = sort(sar_abs(:), 'descend');

num_peaks = 3;
reported = zeros(num_peaks,2);

count = 0;

for k = 1:length(sorted_indices)

    [row,col] = ind2sub(size(sar_abs), sorted_indices(k));

    candidate_x = x_image(col);
    candidate_range = range_image(row);

    % Don't report another peak too close to an existing one
    if count == 0

        count = 1;
        reported(count,:) = [candidate_x candidate_range];

    else

        distance = sqrt( ...
            (reported(1:count,1)-candidate_x).^2 + ...
            (reported(1:count,2)-candidate_range).^2);

        if all(distance > 10)

            count = count + 1;
            reported(count,:) = ...
                [candidate_x candidate_range];

        end
    end

    if count == num_peaks
        break;
    end
end


for n = 1:count

    fprintf("Peak %d: X = %.1f m, Range = %.1f m\n", ...
        n, reported(n,1), reported(n,2));

end


%% Final results

fprintf("\n========================================\n");
fprintf("SAR PROCESSING COMPLETE\n");
fprintf("========================================\n");

fprintf("Carrier frequency  = %.2f GHz\n", fc/1e9);
fprintf("Bandwidth          = %.2f MHz\n", B/1e6);
fprintf("Wavelength         = %.4f m\n", lambda);
fprintf("Range resolution   = %.2f m\n", c/(2*B));
fprintf("Synthetic aperture = %.1f m\n", aperture);
fprintf("Radar positions    = %d\n", Na);
fprintf("Image size         = %d x %d\n", Ny, Nx);
fprintf("FFT size           = %d\n", Nfft);

fprintf("\nExpected target positions:\n");

for n = 1:size(targets,1)

    fprintf("Target %d: X = %6.1f m, Range = %6.1f m, Reflectivity = %.2f\n", ...
        n, targets(n,1), targets(n,2), targets(n,3));

end

fprintf("\nSimulation finished successfully.\n");