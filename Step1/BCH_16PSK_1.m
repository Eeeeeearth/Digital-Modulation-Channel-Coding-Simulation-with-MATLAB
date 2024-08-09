tic;

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

code_bits = step(bchEncoder, bit_stream);

if mod(length(code_bits), log2(M)) ~= 0
    error('The number of bits is not compatible with the symbol size.');
end

% 16PSK
% The bitstream is divided into 4-bit groups for 16PSK modulation
M = 16;
bit_grouping = reshape(code_bits, log2(M), []).'; % Reorganizes the array to match 4 bits per symbol

% Converts a binary bit group to a decimal integer
symbols = bi2de(bit_grouping, 'left-msb');  % 'left-msb' The leftmost bit is the highest bit

% Modulation
mod_bits = pskmod(symbols, M, pi/16, 'gray');

EbNo_dB = 10;
sps = 1;
M = 16;
SNR = EbNo_dB + 10*log10(log2(M)) - 10*log10(sps);

rxSignal = awgn(mod_bits, SNR, 'measured');

% demodulation
received_grouping = pskdemod(rxSignal, M, pi/16, 'gray');

% Converts decimal integers to a binary bit matrix
bit_matrix = de2bi(received_grouping, log2(M), 'left-msb');

% The two-dimensional bit matrix is reconstructed into a one-dimensional bit stream
received_bits = reshape(bit_matrix.', [], 1);  % Use transposition to ensure the correct bit order

% decode
decoded_bits = step(bchDecoder, received_bits);

numErrors = sum(bit_stream ~= decoded_bits);

ber = numErrors / length(bit_stream);
disp(['Bit Error Rate (BER): ', num2str(ber)]);

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

% display constellation
figure;
subplot(1,2,1);
plot(real(mod_bits), imag(mod_bits), 'bo');
title('Modulated Signal Constellation');
xlabel('real');
ylabel('image');
axis square;
grid on;

subplot(1,2,2);
plot(real(rxSignal), imag(rxSignal), 'ro');
title('Received Signal Constellation');
xlabel('real');
ylabel('image');
axis square;
grid on;
toc; % timer end