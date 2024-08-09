tic; % Start the timer

% Read the image and convert it to a byte stream
img = imread('C:\Users\earth\OneDrive\UNSW Course\2024-T2\TELE4652 Mobile & Satellite Comm System\4652 LAB myLOG\Lab4A&B_MATLAB_Code\low_pixel_test.png');
img_bytes = typecast(img(:), 'uint8');

% Pre-allocate the bit stream
bit_stream = zeros(1, length(img_bytes) * 8);

% Convert byte stream to bit stream
for i = 1:length(img_bytes)
    byte = img_bytes(i);
    bits = bitget(byte, 8:-1:1); % Convert bytes to binary bits (most significant bit first)
    bit_stream((i-1)*8+1:i*8) = bits;
end

% (7,4) Hamming code
n = 7;
k = 4;
bit_stream_length = length(bit_stream);
num_blocks = ceil(bit_stream_length / k);
code_bits = zeros(1, num_blocks * n); % Pre-allocate code bits

for i = 1:k:length(bit_stream)
    if i+k-1 <= length(bit_stream)
        data_bits = bit_stream(i:i+k-1);
    else
        data_bits = [bit_stream(i:end), zeros(1, i+k-1-length(bit_stream))]; % Pad with zeros if necessary
    end
    
    encoded_bits = encode(data_bits, n, k, 'hamming/binary'); % Hamming encoding
    
    start_index = ((i-1)/k) * n + 1;
    end_index = start_index + n - 1;
    code_bits(start_index:end_index) = encoded_bits'; % Store encoded bits
end

% BPSK modulation
mod_bits = 2*code_bits - 1; % Map 0 to -1 and 1 to 1

% Raised Cosine filter parameters
rolloff = 0.35; % Roll-off factor
span = 10; % Filter span in symbols
sps = 4; % Samples per symbol

% Create a Root Raised Cosine filter for the transmitter
tx_filter = comm.RaisedCosineTransmitFilter('RolloffFactor', rolloff, 'FilterSpanInSymbols', span, 'OutputSamplesPerSymbol', sps);

% Create a Root Raised Cosine filter for the receiver
rx_filter = comm.RaisedCosineReceiveFilter('RolloffFactor', rolloff, 'FilterSpanInSymbols', span, 'InputSamplesPerSymbol', sps, 'DecimationFactor', 1);

% Pass the modulated bits through the transmit filter
txSignal = tx_filter(mod_bits');

% Define the Eb/No range
EbNo_dB_range = 0:5:20;
BER = zeros(size(EbNo_dB_range)); % Pre-allocate BER array
M = 2; % M-ary (for BPSK, M = 2)

for idx = 1:length(EbNo_dB_range)
    EbNo_dB = EbNo_dB_range(idx);
    SNR = EbNo_dB + 10*log10(log2(M)) - 10*log10(sps);
    
    % Add AWGN to the transmitted signal
    rxSignal = awgn(txSignal, SNR, 'measured');
    
    % Pass the received signal through the receive filter
    rx_filtered_signal = rx_filter(rxSignal);
    
    % Demodulate the received signal
    received_bits = rx_filtered_signal > 0;
    
    % (7,4) Hamming decode
    decoded_bits = zeros(1, length(bit_stream));
    for i = 1:n:length(received_bits)
        if i+n-1 <= length(received_bits)
            encoded_bits = received_bits(i:i+n-1);
        else
            encoded_bits = [received_bits(i:end), zeros(1, i+n-1-length(received_bits))];
        end
        decoded_bits(((i-1)/n)*k+1:((i-1)/n+1)*k) = decode(encoded_bits, n, k, 'hamming/binary');
    end
    
    % Truncate excess bits
    decoded_bits = decoded_bits(1:length(bit_stream));
    
    % Calculate BER
    numErrors = sum(bit_stream ~= decoded_bits);
    BER(idx) = numErrors / length(bit_stream);
end

% Interpolate for a smooth curve
EbNo_dB_fine = linspace(min(EbNo_dB_range), max(EbNo_dB_range), 100); % Finer points for smoother curve
BER_fine = interp1(EbNo_dB_range, BER, EbNo_dB_fine, 'pchip'); % pchip interpolation for smooth curve

% Plot BER vs Eb/No curve
figure;
plot(EbNo_dB_fine, BER_fine, '-','LineWidth', 2);
hold on;
plot(EbNo_dB_range, BER, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'r'); % Plot original points
grid on;
xlabel('Eb/No (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs Eb/No for (7,4) Hamming Code with Raised Cosine Filter');
xlim([0 20]);
ylim([0 1]); % Set y-axis to range from 0 to 1

toc; % End the timer

% Plot the transmitted signal
figure;
subplot(2,1,1); % Create two subplots
plot(txSignal);
title('Transmitted Signal After Raised Cosine Filtering');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

% Plot the received signal (after noise and matched filtering)
subplot(2,1,2);
plot(rx_filtered_signal);
title('Received Signal After AWGN Channel and Raised Cosine Filtering');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;
