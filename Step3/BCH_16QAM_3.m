tic; % Start the timer

% Read the image
img = imread('C:\Users\earth\OneDrive\UNSW Course\2024-T2\TELE4652 Mobile & Satellite Comm System\4652 LAB myLOG\Lab4A&B_MATLAB_Code\low_pixel_test.png');

% Convert the image to a 1D byte stream
img_bytes = typecast(img(:), 'uint8');

% Convert the byte stream to a bit stream
% Pre-allocate the bit stream array
bit_stream = zeros(1, length(img_bytes) * 8);

% Convert the byte stream to a bit stream
for i = 1:length(img_bytes)
    byte = img_bytes(i);
    bits = bitget(byte, 8:-1:1); % Convert the byte to binary bits (most significant bit first)
    bit_stream((i-1)*8+1:i*8) = bits;
end

% Define BCH code parameters
m = 3;          % Degree of the generator polynomial
n = 2^m - 1;    % Codeword length
k = 4;          % Data length
t = 1;          % Number of correctable errors

% Ensure bit_stream is a column vector
if isrow(bit_stream)
    bit_stream = bit_stream.';  % Transpose to a column vector
end

% Create BCH encoder and decoder objects
bchEncoder = comm.BCHEncoder(n, k);
bchDecoder = comm.BCHDecoder(n, k);

% Encode the bit stream
code_bits = step(bchEncoder, bit_stream);

% Ensure the length of the bit stream is a multiple of num_bits_per_symbol
if mod(length(code_bits), log2(M)) ~= 0
    error('The number of bits is not compatible with the symbol size.');
end

% 16PSK modulation
% Group the bit stream into groups of 4 bits for 16PSK modulation
M = 16; % M value for 16PSK
bit_grouping = reshape(code_bits, log2(M), []).'; % Reshape the array to match 4 bits per symbol

% Convert binary bit groups to decimal integers
symbols = bi2de(bit_grouping, 'left-msb');  % 'left-msb' indicates the leftmost bit is the most significant bit

% Perform 16PSK modulation
mod_bits = qammod(symbols, M, 'UnitAveragePower', true);

% Define Signal-to-Noise Ratio (SNR)
EbNo_dB = 15; % Eb/No value in dB
sps = 1; % Samples per symbol
M = 16; % Number of modulation symbols
SNR = EbNo_dB + 10*log10(log2(M)) - 10*log10(sps);

% Transmit through AWGN channel
rxSignal = awgn(mod_bits, SNR, 'measured');

% Demodulate and decode
received_grouping = qamdemod(rxSignal, M, 'UnitAveragePower', true);

% Convert decimal integers to binary bit matrix
bit_matrix = de2bi(received_grouping, log2(M), 'left-msb');

% Reshape the 2D bit matrix back into a 1D bit stream
received_bits = reshape(bit_matrix.', [], 1);  % Transpose to ensure the bit order is correct

% Decode the bit stream
decoded_bits = step(bchDecoder, received_bits);

% Calculate the number of bit errors
numErrors = sum(bit_stream ~= decoded_bits);

% Calculate the Bit Error Rate (BER)
ber = numErrors / length(bit_stream);
disp(['Bit Error Rate (BER): ', num2str(ber)]);

% Convert the bit stream back to a byte stream and reconstruct the image
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

% Display the original image
figure;
subplot(1, 2, 1);
imshow(img);
title('Original Image');

% Display the reconstructed image
subplot(1, 2, 2);
imshow(img_reconstructed);
title('Reconstructed Image');

% Display constellation diagram
figure;
subplot(1,2,1);
plot(real(mod_bits), imag(mod_bits), 'bo');
title('Constellation Diagram of Modulated Signal');
xlabel('In-Phase');
ylabel('Quadrature');
axis square;
grid on;

subplot(1,2,2);
plot(real(rxSignal), imag(rxSignal), 'ro');
title('Constellation Diagram of Received Signal');
xlabel('In-Phase');
ylabel('Quadrature');
axis square;
grid on;

toc; % End the timer and display the elapsed time
