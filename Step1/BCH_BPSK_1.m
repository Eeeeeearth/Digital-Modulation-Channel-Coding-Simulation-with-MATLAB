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

% Define BCH code
m = 3;          % Degree of the generator polynomial
n = 2^m - 1;    % Codeword length
k = 4;          % Data length
t = 1;          % Number of correctable errors

% Ensure bit_stream is a column vector
if isrow(bit_stream)
    bit_stream = bit_stream.';  % Transpose to column vector
end

% Create BCH encoder and decoder objects
bchEncoder = comm.BCHEncoder(n, k);
bchDecoder = comm.BCHDecoder(n, k);

% BCH encoding
code_bits = step(bchEncoder, bit_stream);

% BPSK Modulation
mod_bits = 2*code_bits - 1; % Map 0 to -1 and 1 to 1

% Define channel
EbNo_dB = 10; % Eb/No values in dB
sps = 4; % Samples per symbol
M = 2; % Number of modulation symbols
SNR = EbNo_dB + 10*log10(log2(M)) - 10*log10(sps);

% Transmit through AWGN channel
rxSignal = awgn(mod_bits, SNR, 'measured');

% demodulation
received_bits = rxSignal > 0;

% decode
decoded_bits = step(bchDecoder, received_bits);

% Calculate number of bit errors
numErrors = sum(bit_stream ~= decoded_bits);
ber = numErrors / length(bit_stream);
disp(['Bit Error Rate (BER): ', num2str(ber)]);

% reconstruct image
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

% display origin img
figure;
subplot(1, 2, 1);
imshow(img);
title('Origin Img');

% display reconstructed img
subplot(1, 2, 2);
imshow(img_reconstructed);
title('Reconstructed Img');

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