clear all;
close all;

%% 1. Setting up the modulation parameters

% M-ary PSK modulation orders. The values must be the power of 2.
M_PSK = [2];

% The range of Eb/No values to observe
EbNo_dB = [0:5:20]; % Eb/No values in dB

numBits = 1e6; % Number of bits to simulate. Not used in the simulation.

numSymbolsForEye = 1000; % Number of symbols for eye diagram

% Setting up the filter parameters
rolloff = 0.35;
span = 10; % Filter span in symbols
sps = 4; % Samples per symbol

% Function to simulate and calculate BER
simulateBER = @(M, modType, useFilter, codeType) calculateBER(M, modType, EbNo_dB, numBits, rolloff, span, sps, useFilter, codeType, numSymbolsForEye);

% Channel code types. Do not change this in Lab 4A.
codeTypes = {'None', 'Hamming', 'Convolutional', 'CRC'};

%% 2. Simulate BER for PSK with and without the filters.
BER_PSK = cell(length(M_PSK), length(codeTypes), 2);
for i = 1:length(M_PSK)
    for j = 1:length(codeTypes)
        BER_PSK{i, j, 1} = simulateBER(M_PSK(i), 'psk', false, codeTypes{j});
    end
end

%% 3. Plot results
figure;
for i = 1:length(M_PSK)
    for j = 1:length(codeTypes)
        semilogy(EbNo_dB, BER_PSK{i, j, 1}, 'DisplayName', ['PSK ' num2str(M_PSK(i)) ' ' codeTypes{j} ' No Filter']);
        hold on;
    end
end

grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs. E_b/N_0 for Different Channel Codes');
legend;
hold off;

%% 4. Function to calculate BER
function BER = calculateBER(M, modType, EbNo_dB, numBits, rolloff, span, sps, useFilter, codeType, numSymbolsForEye)
% Pre-allocate a vector to store the BER to be calculated later.
BER = zeros(size(EbNo_dB));

% Run the simulation for each Eb/No value in the vector "EbNo_dB"
for i = 1:length(EbNo_dB)
    % Task 1: Prepare the message to send
    img = imread('C:\Users\earth\OneDrive\UNSW Course\2024-T2\TELE4652 Mobile & Satellite Comm System\4652 LAB myLOG\Lab4A&B_MATLAB_Code\low_pixel_test.png'); % Read the image
    [rows, cols, ~] = size(img);
    
    % Convert each RGB channel to a bit stream
    bitsR = reshape(de2bi(img(:,:,1), 8)', [], 1);
    bitsG = reshape(de2bi(img(:,:,2), 8)', [], 1);
    bitsB = reshape(de2bi(img(:,:,3), 8)', [], 1);
    
    % Initialize encoded bit streams for each channel
    encodedR = bitsR;
    encodedG = bitsG;
    encodedB = bitsB;
    
    % Task 2: Encoding for each channel
    switch codeType
        case 'Hamming'
            fprintf("Encoding the message with Hamming Code...\n");
            encodedR = encodeHamming(bitsR);
            encodedG = encodeHamming(bitsG);
            encodedB = encodeHamming(bitsB);
        case 'Convolutional'
            fprintf("Encoding the message with Convolutional Code...\n");
            encodedR = encodeConvolutional(bitsR);
            encodedG = encodeConvolutional(bitsG);
            encodedB = encodeConvolutional(bitsB);
        case 'CRC'
            fprintf("Encoding the message with CRC Code...\n");
            encodedR = encodeCRC(bitsR);
            encodedG = encodeCRC(bitsG);
            encodedB = encodeCRC(bitsB);
    end

    % Task 3: Modulation for each channel
    txSymbolsR = pskmod(double(encodedR), M, pi/M); % BPSK modulation for R
    txSymbolsG = pskmod(double(encodedG), M, pi/M); % BPSK modulation for G
    txSymbolsB = pskmod(double(encodedB), M, pi/M); % BPSK modulation for B

    % Task 4: Transmit over AWGN channel for each channel
    SNR = EbNo_dB(i) + 10*log10(log2(M)) - 10*log10(sps);
    rxSymbolsR = awgn(txSymbolsR, SNR, 'measured');
    rxSymbolsG = awgn(txSymbolsG, SNR, 'measured');
    rxSymbolsB = awgn(txSymbolsB, SNR, 'measured');
    
    % Task 4: Demodulation for each channel
    rxBitsR = pskdemod(rxSymbolsR, M, pi/M); % BPSK demodulation for R
    rxBitsG = pskdemod(rxSymbolsG, M, pi/M); % BPSK demodulation for G
    rxBitsB = pskdemod(rxSymbolsB, M, pi/M); % BPSK demodulation for B

    % Task 5: Decoding for each channel
    switch codeType
        case 'Hamming'
            fprintf("Decoding the message with Hamming Code...\n");
            decodedR = decodeHamming(rxBitsR);
            decodedG = decodeHamming(rxBitsG);
            decodedB = decodeHamming(rxBitsB);
        case 'Convolutional'
            fprintf("Decoding the message with Convolutional Code...\n");
            decodedR = decodeConvolutional(rxBitsR);
            decodedG = decodeConvolutional(rxBitsG);
            decodedB = decodeConvolutional(rxBitsB);
        case 'CRC'
            fprintf("Decoding the message with CRC Code...\n");
            decodedR = decodeCRC(rxBitsR);
            decodedG = decodeCRC(rxBitsG);
            decodedB = decodeCRC(rxBitsB);
        otherwise
            decodedR = rxBitsR;
            decodedG = rxBitsG;
            decodedB = rxBitsB;
    end

    % Task 6: Calculate the Bit Error Rate (BER) for each channel
    numErrorsR = sum(bitsR ~= decodedR(1:length(bitsR)));
    numErrorsG = sum(bitsG ~= decodedG(1:length(bitsG)));
    numErrorsB = sum(bitsB ~= decodedB(1:length(bitsB)));
    BER(i) = (numErrorsR + numErrorsG + numErrorsB) / (length(bitsR) + length(bitsG) + length(bitsB));

    % Task 7: Showing your results
    receivedR = reshape(decodedR(1:rows*cols*8), 8, [])';
    receivedR = uint8(bi2de(receivedR));
    receivedR = reshape(receivedR, rows, cols);
    
    receivedG = reshape(decodedG(1:rows*cols*8), 8, [])';
    receivedG = uint8(bi2de(receivedG));
    receivedG = reshape(receivedG, rows, cols);
    
    receivedB = reshape(decodedB(1:rows*cols*8), 8, [])';
    receivedB = uint8(bi2de(receivedB));
    receivedB = reshape(receivedB, rows, cols);
    
    % Combine the RGB channels into one image
    receivedImage = cat(3, receivedR, receivedG, receivedB);
    
    % Display the received image
    figure;
    imshow(receivedImage);
    title(['Received Image for E_b/N_0 = ' num2str(EbNo_dB(i)) ' dB']);
end
end

%% Supporting functions

% Manual Hamming (7,4) encoding function
function encodedBits = encodeHamming(bits)
    % Group bits into blocks of 4
    bits = reshape(bits, 4, [])';
    
    % Convert bits to double for matrix multiplication
    bits = double(bits);
    
    % Generator matrix for Hamming (7,4)
    G = [1 0 0 0 1 1 0;
         0 1 0 0 1 0 1;
         0 0 1 0 0 1 1;
         0 0 0 1 1 1 1];
     
    % Encode the bits using matrix multiplication
    encodedBits = mod(bits * G, 2);
    
    % Reshape to a column vector and convert back to logical
    encodedBits = logical(encodedBits');
    encodedBits = encodedBits(:);
end

% Manual Hamming (7,4) decoding function
function decodedBits = decodeHamming(bits)
    % Reshape bits into blocks of 7
    bits = reshape(bits, 7, [])';
    
    % Convert bits to double for matrix multiplication
    bits = double(bits);
    
    % Parity-check matrix for Hamming (7,4)
    H = [1 1 1 0 1 0 0;
         1 0 1 1 0 1 0;
         1 1 0 1 0 0 1];
     
    % Syndrome calculation
    syndrome = mod(bits * H', 2);
    
    % Error correction
    errorTable = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1];
    errorPositions = [7, 6, 5, 4, 3, 2, 1];  % Correct error positions (index in bits)
    for i = 1:size(bits, 1)
        % Find the position of the error based on the syndrome
        [~, idx] = ismember(syndrome(i, :), errorTable, 'rows');
        if idx > 1
            errorPos = errorPositions(idx - 1);
            bits(i, errorPos) = mod(bits(i, errorPos) + 1, 2);
        end
    end
    
    % Extract the original 4 bits and convert back to logical
    decodedBits = logical(bits(:, 1:4)');
    decodedBits = decodedBits(:);
end

% Convolutional encoding function
function encodedBits = encodeConvolutional(bits)
    trellis = poly2trellis(7, [171 133]); % Define a rate 1/2 trellis
    encodedBits = convenc(bits, trellis); % Perform convolutional encoding
end

% Convolutional decoding function
function decodedBits = decodeConvolutional(bits)
    trellis = poly2trellis(7, [171 133]); % Define a rate 1/2 trellis
    decodedBits = vitdec(bits, trellis, 34, 'trunc', 'hard'); % Perform Viterbi decoding
end

% CRC encoding function
function encodedBits = encodeCRC(bits)
    crcGen = comm.CRCGenerator('Polynomial', [1 0 0 1]); % Simple CRC polynomial x^3 + 1
    encodedBits = crcGen(bits);
end

% CRC decoding function
function decodedBits = decodeCRC(bits)
    crcDet = comm.CRCDetector('Polynomial', [1 0 0 1]); % Simple CRC polynomial x^3 + 1
    [decodedBits, ~] = crcDet(bits);
end
