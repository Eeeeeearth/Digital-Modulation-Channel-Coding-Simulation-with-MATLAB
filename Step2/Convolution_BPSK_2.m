tic; % timer start

img = imread('C:\Users\earth\OneDrive\UNSW Course\2024-T2\TELE4652 Mobile & Satellite Comm System\4652 LAB LOG\low_pixel_test.png');

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

% Define convolutional code
trellis = poly2trellis(7, [171 133]);  % The encoder structure is defined using the poly2trellis function

% convolutional coding
code_bits = convenc(bit_stream, trellis);

% BPSK modulation
mod_bits = 2*code_bits - 1; % Map 0 to -1 and 1 to 1

% Define Eb/No values
EbNo_dB = 0:5:20; % Eb/No values in dB
sps = 4; % samples per symbol
M = 2; % M-ary
ber = zeros(1, length(EbNo_dB)); % Pre-allocate BER array

for idx = 1:length(EbNo_dB)
    % Calculate SNR
    SNR = EbNo_dB(idx) + 10*log10(log2(M)) - 10*log10(sps);

    % AWGN channel
    rxSignal = awgn(mod_bits, SNR, 'measured');

    % demodulation
    received_bits = rxSignal > 0;

    % Viterbi decode
    decoded_bits = vitdec(received_bits, trellis, 34, 'trunc', 'hard');  % 34 here is the tracking depth of the decoder

    % cut exceed bits
    decoded_bits = decoded_bits(1:length(bit_stream));

    % calculate BER
    numErrors = sum(bit_stream ~= decoded_bits);
    ber(idx) = numErrors / length(bit_stream);
end

% Plot BER vs Eb/No
figure;
semilogy(EbNo_dB, ber, 'b-o');
grid on;
xlabel('Eb/No (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs Eb/No');

% reconstruct img
num_bytes = length(decoded_bits) / 8;
reconstructed_bytes = uint8(zeros(1, num_bytes));
for i = 1:num_bytes
    bits = int32(decoded_bits((i-1)*8 + 1:i*8));
    reconstructed_bytes(i) = bitshift(bits(1), 7) + ...
                             bitshift(bits(2), 6) + ...
                             bitshift(bits(3), 5) + ...
                             bitshift(bits(4), 4) + ...
                             bitshift(bits(5), 3) + ...
                             bitshift(bits(6), 2) + ...
                             bitshift(bits(7), 1) + ...
                             bits(8);
end

img_reconstructed = typecast(reconstructed_bytes, class(img));
img_reconstructed = reshape(img_reconstructed, size(img));

% display original img
figure;
subplot(1, 2, 1);
imshow(img);
title('Original Image');

% display reconstructed img
subplot(1, 2, 2);
imshow(img_reconstructed);
title('Reconstructed Image');

% Display Signal
figure;
subplot(3, 1, 1);
stem(mod_bits(1:50), 'filled');
title('Modulated Signal');
xlabel('Bit Number');
ylabel('Amplitude');

subplot(3, 1, 2);
stem(rxSignal(1:50*sps), 'filled');
title('Received Signal after AWGN');
xlabel('Symbol Number');
ylabel('Amplitude');

subplot(3, 1, 3);
stem(received_bits(1:50), 'filled');
title('Decoded Bits');
xlabel('Bit Number');
ylabel('Amplitude');

toc; % timer end
