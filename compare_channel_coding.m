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

% Common parameters
EbNo_dB_range = 0:5:20;
sps = 4; % samples per symbol
M = 2; % M-ary

% 1. (7,4) Hamming code
n = 7;
k = 4;
num_blocks = ceil(length(bit_stream) / k);
code_bits_hamming = zeros(1, num_blocks * n);

for i = 1:k:length(bit_stream)
    if i+k-1 <= length(bit_stream)
        data_bits = bit_stream(i:i+k-1);
    else
        data_bits = [bit_stream(i:end), zeros(1, i+k-1-length(bit_stream))]; % insert zero
    end
    
    encoded_bits = encode(data_bits, n, k, 'hamming/binary'); % Hamming coding
    start_index = ((i-1)/k) * n + 1;
    end_index = start_index + n - 1;
    code_bits_hamming(start_index:end_index) = encoded_bits'; % Store code_bits
end

% 2. Convolutional coding
trellis = poly2trellis(7, [171 133]);
code_bits_conv = convenc(bit_stream, trellis);

% 3. BCH coding
m = 3;
n_bch = 2^m - 1;
k_bch = 4;
bchEncoder = comm.BCHEncoder(n_bch, k_bch);
bchDecoder = comm.BCHDecoder(n_bch, k_bch);
code_bits_bch = step(bchEncoder, bit_stream(:));

% Allocate BER arrays
BER_hamming = zeros(size(EbNo_dB_range));
BER_conv = zeros(size(EbNo_dB_range));
BER_bch = zeros(size(EbNo_dB_range));

for idx = 1:length(EbNo_dB_range)
    EbNo_dB = EbNo_dB_range(idx);
    SNR = EbNo_dB + 10*log10(log2(M)) - 10*log10(sps);

    % Hamming code simulation
    mod_bits_hamming = 2*code_bits_hamming - 1; % BPSK modulation
    rxSignal_hamming = awgn(mod_bits_hamming, SNR, 'measured');
    received_bits_hamming = rxSignal_hamming > 0;
    decoded_bits_hamming = zeros(1, length(bit_stream));
    for i = 1:n:length(received_bits_hamming)
        if i+n-1 <= length(received_bits_hamming)
            encoded_bits = received_bits_hamming(i:i+n-1);
        else
            encoded_bits = [received_bits_hamming(i:end), zeros(1, i+n-1-length(received_bits_hamming))];
        end
        decoded_bits_hamming(((i-1)/n)*k+1:((i-1)/n+1)*k) = decode(encoded_bits, n, k, 'hamming/binary');
    end
    decoded_bits_hamming = decoded_bits_hamming(1:length(bit_stream));
    numErrors_hamming = sum(bit_stream ~= decoded_bits_hamming);
    BER_hamming(idx) = numErrors_hamming / length(bit_stream);

    % Convolutional code simulation
    mod_bits_conv = 2*code_bits_conv - 1; % BPSK modulation
    rxSignal_conv = awgn(mod_bits_conv, SNR, 'measured');
    received_bits_conv = rxSignal_conv > 0;
    decoded_bits_conv = vitdec(received_bits_conv, trellis, 34, 'trunc', 'hard');
    numErrors_conv = sum(bit_stream ~= decoded_bits_conv);
    BER_conv(idx) = numErrors_conv / length(bit_stream);

    % BCH code simulation
    mod_bits_bch = 2*code_bits_bch - 1; % BPSK modulation
    rxSignal_bch = awgn(mod_bits_bch, SNR, 'measured');
    received_bits_bch = rxSignal_bch > 0;
    decoded_bits_bch = step(bchDecoder, received_bits_bch);
    numErrors_bch = sum(bit_stream ~= decoded_bits_bch(:)');
    BER_bch(idx) = numErrors_bch / length(bit_stream);
end

% Interpolation for smooth curves
EbNo_dB_fine = linspace(min(EbNo_dB_range), max(EbNo_dB_range), 100);
BER_hamming_fine = interp1(EbNo_dB_range, BER_hamming, EbNo_dB_fine, 'pchip');
BER_conv_fine = interp1(EbNo_dB_range, BER_conv, EbNo_dB_fine, 'pchip');
BER_bch_fine = interp1(EbNo_dB_range, BER_bch, EbNo_dB_fine, 'pchip');

% Plot BER vs Eb/No for all three coding schemes
figure;
plot(EbNo_dB_fine, BER_hamming_fine, '-r','LineWidth', 2);
hold on;
plot(EbNo_dB_fine, BER_conv_fine, '-g','LineWidth', 2);
plot(EbNo_dB_fine, BER_bch_fine, '-b','LineWidth', 2);
plot(EbNo_dB_range, BER_hamming, 'or', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
plot(EbNo_dB_range, BER_conv, 'og', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
plot(EbNo_dB_range, BER_bch, 'ob', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
grid on;
xlabel('Eb/No (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs Eb/No for Different Coding Schemes');
legend('(7,4) Hamming Code', 'Convolutional Code', 'BCH Code');
xlim([0 20]);
ylim([0 1]);

toc; % timer end
