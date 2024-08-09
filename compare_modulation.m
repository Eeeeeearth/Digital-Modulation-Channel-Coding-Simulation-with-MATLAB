tic; % timer start

img = imread('C:\Users\earth\OneDrive\UNSW Course\2024-T2\TELE4652 Mobile & Satellite Comm System\4652 LAB myLOG\Lab4A&B_MATLAB_Code\low_pixel_test.png');

% Converts the image to byte stream
img_bytes = typecast(img(:), 'uint8');

% Pre-allocate bit_stream
bit_stream = zeros(1, length(img_bytes) * 8);

% Converts byte stream to bit stream
for i = 1:length(img_bytes)
    byte = img_bytes(i);
    bits = bitget(byte, 8:-1:1); % Convert bytes to binary bits (high first)
    bit_stream((i-1)*8+1:i*8) = bits;
end

% Ensure bit_stream is a column vector
if isrow(bit_stream)
    bit_stream = bit_stream.';  % Transpose to column vector
end

% Define BCH code parameters
m = 3;          % Degree of the generator polynomial
n = 2^m - 1;    % Codeword length
k = 4;          % Data length
t = 1;          % Number of correctable errors

% Create BCH encoder and decoder objects
bchEncoder = comm.BCHEncoder(n, k);
bchDecoder = comm.BCHDecoder(n, k);

% Common parameters
EbNo_dB_range = 0:5:20;

% Allocate BER arrays
BER_bpsk = zeros(size(EbNo_dB_range));
BER_16psk = zeros(size(EbNo_dB_range));
BER_16qam = zeros(size(EbNo_dB_range));

% 1. BPSK Modulation
for idx = 1:length(EbNo_dB_range)
    EbNo_dB = EbNo_dB_range(idx);
    SNR = EbNo_dB + 10*log10(log2(2)) - 10*log10(4);

    code_bits = step(bchEncoder, bit_stream);

    % BPSK Modulation
    mod_bits = 2*code_bits - 1; % Map 0 to -1 and 1 to 1

    % AWGN channel
    rxSignal = awgn(mod_bits, SNR, 'measured');

    % Demodulation
    received_bits = rxSignal > 0;

    % Release and reset decoder before usage
    release(bchDecoder);
    % Decode
    decoded_bits = step(bchDecoder, received_bits);

    % Calculate BER
    numErrors = sum(bit_stream ~= decoded_bits);
    BER_bpsk(idx) = numErrors / length(bit_stream);
end

% 2. 16PSK Modulation
M_psk = 16; % 16PSK modulation
for idx = 1:length(EbNo_dB_range)
    EbNo_dB = EbNo_dB_range(idx);
    SNR = EbNo_dB + 10*log10(log2(M_psk)) - 10*log10(1);

    code_bits = step(bchEncoder, bit_stream);

    if mod(length(code_bits), log2(M_psk)) ~= 0
        error('The number of bits is not compatible with the symbol size.');
    end

    % 16PSK Modulation
    bit_grouping = reshape(code_bits, log2(M_psk), []).'; % Reorganizes the array to match 4 bits per symbol
    symbols = bi2de(bit_grouping, 'left-msb');  % Converts a binary bit group to a decimal integer
    mod_bits = pskmod(symbols, M_psk, pi/16, 'gray');  % Modulation

    % AWGN channel
    rxSignal = awgn(mod_bits, SNR, 'measured');

    % Demodulation
    received_grouping = pskdemod(rxSignal, M_psk, pi/16, 'gray');
    bit_matrix = de2bi(received_grouping, log2(M_psk), 'left-msb'); % Converts decimal integers to a binary bit matrix
    received_bits = reshape(bit_matrix.', [], 1);  % The two-dimensional bit matrix is reconstructed into a one-dimensional bit stream

    % Release and reset decoder before usage
    release(bchDecoder);
    % Decode
    decoded_bits = step(bchDecoder, received_bits);

    % Calculate BER
    numErrors = sum(bit_stream ~= decoded_bits);
    BER_16psk(idx) = numErrors / length(bit_stream);
end

% 3. 16QAM Modulation
M_qam = 16; % 16QAM modulation
for idx = 1:length(EbNo_dB_range)
    EbNo_dB = EbNo_dB_range(idx);
    SNR = EbNo_dB + 10*log10(log2(M_qam)) - 10*log10(1);

    code_bits = step(bchEncoder, bit_stream);

    if mod(length(code_bits), log2(M_qam)) ~= 0
        error('The number of bits is not compatible with the symbol size.');
    end

    % 16QAM Modulation
    bit_grouping = reshape(code_bits, log2(M_qam), []).'; % Reorganizes the array to match 4 bits per symbol
    symbols = bi2de(bit_grouping, 'left-msb');  % Converts a binary bit group to a decimal integer
    mod_bits = qammod(symbols, M_qam, 'UnitAveragePower', true);  % Modulation

    % AWGN channel
    rxSignal = awgn(mod_bits, SNR, 'measured');

    % Demodulation
    received_grouping = qamdemod(rxSignal, M_qam, 'UnitAveragePower', true);
    bit_matrix = de2bi(received_grouping, log2(M_qam), 'left-msb'); % Converts decimal integers to a binary bit matrix
    received_bits = reshape(bit_matrix.', [], 1);  % The two-dimensional bit matrix is reconstructed into a one-dimensional bit stream

    % Release and reset decoder before usage
    release(bchDecoder);
    % Decode
    decoded_bits = step(bchDecoder, received_bits);

    % Calculate BER
    numErrors = sum(bit_stream ~= decoded_bits);
    BER_16qam(idx) = numErrors / length(bit_stream);
end

% Interpolation for smooth curves
EbNo_dB_fine = linspace(min(EbNo_dB_range), max(EbNo_dB_range), 100);
BER_bpsk_fine = interp1(EbNo_dB_range, BER_bpsk, EbNo_dB_fine, 'pchip');
BER_16psk_fine = interp1(EbNo_dB_range, BER_16psk, EbNo_dB_fine, 'pchip');
BER_16qam_fine = interp1(EbNo_dB_range, BER_16qam, EbNo_dB_fine, 'pchip');

% Plot BER vs Eb/No for all three modulation schemes
figure;
plot(EbNo_dB_fine, BER_bpsk_fine, '-r','LineWidth', 2);
hold on;
plot(EbNo_dB_fine, BER_16psk_fine, '-g','LineWidth', 2);
plot(EbNo_dB_fine, BER_16qam_fine, '-b','LineWidth', 2);
plot(EbNo_dB_range, BER_bpsk, 'or', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
plot(EbNo_dB_range, BER_16psk, 'og', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
plot(EbNo_dB_range, BER_16qam, 'ob', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
grid on;
xlabel('Eb/No (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs Eb/No for Different Modulation Schemes');
legend('BPSK', '16PSK', '16QAM');
xlim([0 20]);
ylim([0 1]);

toc; % timer end
