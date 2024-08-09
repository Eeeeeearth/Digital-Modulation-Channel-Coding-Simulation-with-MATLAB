tic; % timer start

img = imread('C:\Users\earth\OneDrive\UNSW Course\2024-T2\TELE4652 Mobile & Satellite Comm System\4652 LAB myLOG\Lab4A&B_MATLAB_Code\low_pixel_test.png');

% converts the image to byte stream
img_bytes = typecast(img(:), 'uint8');

% pre-allocate bit_stream
bit_stream = zeros(1, length(img_bytes) * 8);

% converts byte stream to bit stream
for i = 1:length(img_bytes)
    byte = img_bytes(i);
    bits = bitget(byte, 8:-1:1); % convert bytes to binary bits (high first)
    bit_stream((i-1)*8+1:i*8) = bits;
end

% Ensure bit_stream is a column vector
if isrow(bit_stream)
    bit_stream = bit_stream.';  % Transpose to column vector
end

% Define BCH code
m = 3;          % Degree of the generator polynomial
n = 2^m - 1;    % Codeword length
k = 4;          % Data length
t = 1;          % Number of correctable errors

% Create BCH encoder and decoder objects
bchEncoder = comm.BCHEncoder(n, k);
bchDecoder = comm.BCHDecoder(n, k);

% Define Eb/No range
EbNo_dB_range = 0:5:20;
BER = zeros(size(EbNo_dB_range)); % Pre-allocate BER array
sps = 4; % Samples per symbol
M = 2; % Number of modulation symbols

for idx = 1:length(EbNo_dB_range)
    EbNo_dB = EbNo_dB_range(idx);
    SNR = EbNo_dB + 10*log10(log2(M)) - 10*log10(sps);

    % BCH encoding
    code_bits = step(bchEncoder, bit_stream);

    % BPSK Modulation
    mod_bits = 2*code_bits - 1; % Map 0 to -1 and 1 to 1

    % Transmit through AWGN channel
    rxSignal = awgn(mod_bits, SNR, 'measured');

    % demodulation
    received_bits = rxSignal > 0;

    % decode
    decoded_bits = step(bchDecoder, received_bits);

    % Calculate number of bit errors
    numErrors = sum(bit_stream ~= decoded_bits);
    BER(idx) = numErrors / length(bit_stream);
end

% Interpolation for smooth curve
EbNo_dB_fine = linspace(min(EbNo_dB_range), max(EbNo_dB_range), 100); % finer points for smoother curve
BER_fine = interp1(EbNo_dB_range, BER, EbNo_dB_fine, 'pchip'); % pchip interpolation for smooth curve

% Plot BER vs Eb/No with smooth curve
figure;
plot(EbNo_dB_fine, BER_fine, '-','LineWidth', 2);
hold on;
plot(EbNo_dB_range, BER, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'r'); % plot original points
grid on;
xlabel('Eb/No (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs Eb/No with BCH Coding');
xlim([0 20]);
ylim([0 1]); % Set y-axis to range from 0 to 1

toc; % timer end
