tic; % timer start

img = imread('C:\Users\earth\OneDrive\UNSW Course\2024-T2\TELE4652 Mobile & Satellite Comm System\4652 LAB myLOG\Lab4A&B_MATLAB_Code\low_pixel_test.png');

img_bytes = typecast(img(:), 'uint8');

bit_stream = zeros(1, length(img_bytes) * 8);

for i = 1:length(img_bytes)
    byte = img_bytes(i);
    bits = bitget(byte, 8:-1:1);
    bit_stream((i-1)*8+1:i*8) = bits;
end

m = 3;
n = 2^m - 1;
k = 4; 
t = 1;

if isrow(bit_stream)
    bit_stream = bit_stream.'; 
end

bchEncoder = comm.BCHEncoder(n, k);
bchDecoder = comm.BCHDecoder(n, k);

% Define Eb/No range
EbNo_dB_range = 0:5:20;
BER = zeros(size(EbNo_dB_range)); % Pre-allocate BER array
sps = 1; % Samples per symbol
M = 16; % 16PSK modulation

for idx = 1:length(EbNo_dB_range)
    EbNo_dB = EbNo_dB_range(idx);
    SNR = EbNo_dB + 10*log10(log2(M)) - 10*log10(sps);
    
    code_bits = step(bchEncoder, bit_stream);
    
    if mod(length(code_bits), log2(M)) ~= 0
        error('The number of bits is not compatible with the symbol size.');
    end

    % 16PSK Modulation
    bit_grouping = reshape(code_bits, log2(M), []).'; % Reorganizes the array to match 4 bits per symbol
    symbols = bi2de(bit_grouping, 'left-msb');  % Converts a binary bit group to a decimal integer
    mod_bits = pskmod(symbols, M, pi/16, 'gray');  % Modulation

    % AWGN channel
    rxSignal = awgn(mod_bits, SNR, 'measured');

    % demodulation
    received_grouping = pskdemod(rxSignal, M, pi/16, 'gray');
    bit_matrix = de2bi(received_grouping, log2(M), 'left-msb'); % Converts decimal integers to a binary bit matrix
    received_bits = reshape(bit_matrix.', [], 1);  % The two-dimensional bit matrix is reconstructed into a one-dimensional bit stream

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
title('BER vs Eb/No with 16PSK and BCH Coding');
xlim([0 20]);
ylim([0 1]); % Set y-axis to range from 0 to 1

toc; % timer end
