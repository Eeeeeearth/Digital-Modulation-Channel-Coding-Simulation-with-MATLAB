tic; % Start timer

img = imread('C:\Users\earth\OneDrive\UNSW Course\2024-T2\TELE4652 Mobile & Satellite Comm System\4652 LAB LOG\low_pixel_test.png'); % Read image

% Convert image to a one-dimensional byte stream
img_bytes = typecast(img(:), 'uint8');

% Preallocate bit_stream
bit_stream = zeros(1, length(img_bytes) * 8);

% Convert byte stream to bit stream
for i = 1:length(img_bytes)
    byte = img_bytes(i);
    bits = bitget(byte, 8:-1:1); % Convert byte to binary bits (MSB first)
    bit_stream((i-1)*8+1:i*8) = bits;
end

% Define BCH code parameters
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

% Encoding
code_bits = step(bchEncoder, bit_stream);

% Ensure the length of the bit stream is a multiple of the number of bits per symbol
M = 16; % 16PSK modulation order
if mod(length(code_bits), log2(M)) ~= 0
    error('The number of bits is not compatible with the symbol size.');
end

% Reshape the bit stream into groups of 4 bits for 16PSK modulation
bit_grouping = reshape(code_bits, log2(M), []).'; % Reshape array to match 4 bits per symbol

% Convert binary bit groups to decimal integers
symbols = bi2de(bit_grouping, 'left-msb');  % 'left-msb' indicates the leftmost bit is the most significant

% Perform 16PSK modulation
mod_bits = pskmod(symbols, M, pi/16, 'gray');

% Define Root Raised Cosine Filter parameters
rolloff = 0.25;
span = 10;
sps = 4; % Samples per symbol
rrcFilter = rcosdesign(rolloff, span, sps);

% Filter the signal using Root Raised Cosine Filter (Tx filter)
txSignal = upfirdn(mod_bits, rrcFilter, sps, 1);

% Define Eb/No values
EbNo_dB = 0:5:20; % Eb/No values in dB
ber = zeros(1, length(EbNo_dB)); % Preallocate BER array

for idx = 1:length(EbNo_dB)
    % Calculate SNR
    SNR = EbNo_dB(idx) + 10*log10(log2(M)) - 10*log10(sps);

    % Transmit through AWGN channel
    rxSignal = awgn(txSignal, SNR, 'measured');

    % Filter the received signal using Root Raised Cosine Filter (Rx filter)
    rxFilteredSignal = upfirdn(rxSignal, rrcFilter, 1, sps);
    rxFilteredSignal = rxFilteredSignal(span+1:end-span); % Remove filter transients

    % Demodulate
    received_grouping = pskdemod(rxFilteredSignal, M, pi/16, 'gray');

    % Convert decimal integers to binary bit matrix
    bit_matrix = de2bi(received_grouping, log2(M), 'left-msb');

    % Reshape the 2D bit matrix back into a 1D bit stream
    received_bits = reshape(bit_matrix.', [], 1);  % Use transpose to ensure correct bit order

    % Decode
    decoded_bits = step(bchDecoder, received_bits);

    % Calculate number of bit errors
    numErrors = sum(bit_stream ~= decoded_bits);

    % Calculate BER
    ber[idx] = numErrors / length(bit_stream);
end

% Plot BER vs Eb/No
figure;
semilogy(EbNo_dB, ber, 'b-o');
grid on;
xlabel('Eb/No (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs Eb/No');

% Reconstruct image
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

% Display original image
figure;
subplot(1, 2, 1);
imshow(img);
title('Original Image');

% Display reconstructed image
subplot(1, 2, 2);
imshow(img_reconstructed);
title('Reconstructed Image');

% Display constellations
figure;
subplot(1, 2, 1);
plot(real(mod_bits), imag(mod_bits), 'bo');
title('Constellation of Modulated Signal');
xlabel('In-Phase');
ylabel('Quadrature');
axis square;
grid on;

subplot(1, 2, 2);
plot(real(rxSignal), imag(rxSignal), 'ro');
title('Constellation of Received Signal');
xlabel('In-Phase');
ylabel('Quadrature');
axis square;
grid on;

toc; % End timer
