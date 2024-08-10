clear;
clc;
% Parameter setup
% M-ary PSK modulation orders. The numbers here must be the power of 2.
M_PSK = [16]; 
% M-ary QAM modulation orders. The numbers here must be the power of 2.
M_QAM = [16]; 
EbNo_dB = [0:5:20]; % Eb/No values in dB
numSymbolsForEye = 1000; % Number of symbols for eye diagram

% Filter parameters
rolloff = 0.35;
span = 10; % Filter span in symbols
sps = 4; % Samples per symbol

% Function to simulate and calculate BER
simulateBER = @(M, modType, useFilter, codeType) calculateBER(M, modType, EbNo_dB, rolloff, span, sps, useFilter, codeType, numSymbolsForEye);

% Channel code types
codeTypes = {'None'};

% Simulate BER for PSK with and without SRRC filters and different codes
BER_PSK = cell(length(M_PSK), length(codeTypes), 2);
for i = 1:length(M_PSK)
    for j = 1:length(codeTypes)
        BER_PSK{i, j, 1} = simulateBER(M_PSK(i), 'psk', false, codeTypes{j});
        BER_PSK{i, j, 2} = simulateBER(M_PSK(i), 'psk', true, codeTypes{j});
    end
end

% Simulate BER for QAM with and without SRRC filters and different codes
BER_QAM = cell(length(M_QAM), length(codeTypes), 2);
for i = 1:length(M_QAM)
    for j = 1:length(codeTypes)
        BER_QAM{i, j, 1} = simulateBER(M_QAM(i), 'qam', false, codeTypes{j});
        BER_QAM{i, j, 2} = simulateBER(M_QAM(i), 'qam', true, codeTypes{j});
    end
end

% Plot results
figure;
for i = 1:length(M_PSK)
    for j = 1:length(codeTypes)
        semilogy(EbNo_dB, BER_PSK{i, j, 1}, 'DisplayName', ['PSK ' num2str(M_PSK(i)) ' ' codeTypes{j} ' No Filter']);
        hold on;
        semilogy(EbNo_dB, BER_PSK{i, j, 2}, '--', 'DisplayName', ['PSK ' num2str(M_PSK(i)) ' ' codeTypes{j} ' With Filter']);
    end
end

for i = 1:length(M_QAM)
    for j = 1:length(codeTypes)
        semilogy(EbNo_dB, BER_QAM{i, j, 1}, 'DisplayName', ['QAM ' num2str(M_QAM(i)) ' ' codeTypes{j} ' No Filter']);
        hold on;
        semilogy(EbNo_dB, BER_QAM{i, j, 2}, '--', 'DisplayName', ['QAM ' num2str(M_QAM(i)) ' ' codeTypes{j} ' With Filter']);
    end
end

grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs. E_b/N_0 for PSK and QAM Modulations with and without Filters');
legend;
hold off;

% Function to calculate BER
function BER = calculateBER(M, modType, EbNo_dB, rolloff, span, sps, useFilter, codeType, numSymbolsForEye)
BER = zeros(size(EbNo_dB));
tx_constDiagram = comm.ConstellationDiagram('SamplesPerSymbol',sps, ...
    'SymbolsToDisplaySource','Property','SymbolsToDisplay',100);
tx_constDiagram.ShowReferenceConstellation = true;
rx_constDiagram = comm.ConstellationDiagram('SamplesPerSymbol',sps, ...
    'SymbolsToDisplaySource','Property','SymbolsToDisplay',100);
rx_constDiagram.ShowReferenceConstellation = true;
for i = 1:length(EbNo_dB)
    % Lab 4B, Task 1: Prepare the message to send.
    % In this task, you are required to prepare the "message" to transmit
    % during the simulation. This message is the image named
    % "image4Exp4.jpg". Since you are required to transmit an image, you
    % may consider the following issues:
    % - how to read this image into MATLAB?
    % - how the image is stored in MATLAB?
    % - Do I need to convert the read image to other form(s) or data types
    % to align with the next step of process?

    % Hint: 1. you may use the built-in function "imread" to read an image
    % file to MATLAB.
    % 2. You may re-use your codes if you have completed Lab 4A.

    %%% Your codes for Task 1 starts from below: %%%
    
    img = imread('C:\Users\earth\OneDrive\UNSW Course\2024-T2\TELE4652 Mobile & Satellite Comm System\4652 LAB myLOG\Lab4A&B_MATLAB_Code\low_pixel_test.png');   % 读取图像
    grayImg = rgb2gray(img);          % 转为灰度图像
    [rows, cols] = size(grayImg);
    bits = de2bi(grayImg(:), 8, 'left-msb');   % 转换为二进制比特流
    bits = bits(:);   % 展平成一维数组

    %%% Your codes for Task 1 ends here. %%%

    % Apply channel coding
    % Do not apply any channel coding for now.
   switch codeType
        case 'Hamming'
            fprintf("Encoding the message with Hamming Code...\n");
        case 'codeTypeA'
            fprintf("Encoding the message with %s...\n", 'codeTypeA');
        case 'codeTypeB'
            fprintf("Encoding the message with %s...\n", 'codeTypeB');
       otherwise
            % no channel coding used.
            codedBits = bits;
    end

    % % Lab 4B, Task 2: Modulation.
    % In this task, you are required to implement the following types of 
    % digtial modulation that are widely used in modern communication 
    % systems:
    % - BPSK
    % - M-ary PSK
    % - M-ary QAM

    % For each of the modulation types, your implementation should allow
    % you to configure the following properties:
    % - A fixed Phase Offset
    % - natrual binary labelling or Gray labelling (see Hint 2 for more details)
    %

    % The result of your modulation should be an array (or a matrix) of 
    % complex valued symbols that are ready for the next step of
    % processing.

    % Hints: 
    % 1. When implementing the BPSK modulator, it is recommended to create your own 
    % function in a separate .m file and call this function here when you
    % need to modulate the bits with BPSK. The reason is that having a
    % function defined in a separate file allows you to test this function
    % easily. It is also recommended to do the same when implementing the 
    % BPSK demodulator in the later task.
    %
    % 2. What modulation does is essentially mapping a bit (or a group of
    % bits) into a number on a complex plane. However, there are many ways to create
    % such a mapping and you need to know which bit (or a group of bits) will be mapped
    % to which complex number - the bit labelling. You are recommended to
    % do a bit research about the Gray labelling (or sometimes referred to
    % as "Gray Coding").
    % 
    % 3. You are allowed to implement the M-ary PSK and QAM (but not BPSK) with the
    % built-in PSK and QAM modulators/demodulators from the Communication Toolbox. However, you should
    % read through their manuals/documentation and their theories behind.

    if strcmp(modType, 'psk')
        % --- Your response to Task 2 for BPSK/M-ary PSK modulation starts from below: ---
        
        if M == 2
            symbols = 2*bits-1; % BPSK modulation
        else
            symbols = pskmod(bits, M, pi/M, 'gray'); % M-ary PSK modulation with Gray coding
        end

        % --- Your response to Task 2 for PSK modulation ends here. ---
    elseif strcmp(modType, 'qam')
        % --- Your response to Task 2 for M-ary QAM starts from below: ---
        
        symbols = qammod(bits, M, 'gray'); % M-ary QAM modulation with Gray coding
        
        % --- Your response to Task 2 for M-ary QAM ends here. ---    
    end

    % Lab 4B, Task 3: Pulse-shaping filter at the transmitter
    % A pulse shaping filter is a device or algorithm that smooths the 
    % transition of signals between different states, thereby reducing the 
    % bandwidth of the transmitted signal and mitigating inter-symbol 
    % interference (ISI).
    % 
    % In this task, you are required to investigate different types of 
    % pulse-shaping filters (choosing from Gaussian or Root Raised Cosine) 
    % and how these filters affect the overall performance of the 
    % communication system.
    % 

    % Hint:
    % 1. You can call the built-in pulse-shaping functions in the Communicaiton (or 
    % Signal Processing) Toolbox. Feel free to use the MATLAB commands 
    % "doc " or "help " to see how to use the built-in functions correctly. 
    % 2. Plot your pulse-shaping filter in both time and frequency domain
    % to check if your designed filter is exactly what you expect.

    if useFilter
        % --- Your response to Task 3 start from below: ---
        
        txFilter = comm.RaisedCosineTransmitFilter( ...
            'RolloffFactor', rolloff, ...
            'FilterSpanInSymbols', span, ...
            'OutputSamplesPerSymbol', sps);
        txSignal = txFilter(symbols);
        
        % --- Your response ends here. ---
    else
        % No pulse-shaping filter applied
        txSignal = symbols;
    end

    % Add AWGN noise
    SNR = EbNo_dB(i) + 10*log10(log2(M)) - 10*log10(sps);
    rxSignal = awgn(txSignal, SNR, 'measured');

    % Lab 4B, Task 4: Pulse-shaping filter at the receiver
    % In this task, you are required to apply the pulse-shaping filter of
    % your choice from Task 3 at the receiver.
    % 

    % Hint:
    % 1. You can call the built-in functions in the Communicaiton (or 
    % Signal Processing) Toolbox. Feel free to use the MATLAB commands 
    % "doc " or "help " to see how to use the built-in functions correctly. 
    % 2. Plot your pulse-shaping filter in both time and frequency domain
    % to check if your designed filter is exactly what you expect.

    if useFilter
        % --- Your response to Task 4 start from below: ---
        
        rxFilter = comm.RaisedCosineReceiveFilter( ...
            'RolloffFactor', rolloff, ...
            'FilterSpanInSymbols', span, ...
            'InputSamplesPerSymbol', sps, ...
            'DecimationFactor', sps);
        rxSignal = rxFilter(rxSignal);
        
        % --- Your response ends here. ---
    end

    % Lab 4B, Task 5: Demodulation
    %
    % In this task, you are required to demodulate the received symbols
    % (after the pulse-shaping) and convert them back to bits. You may
    % follow the hints in Task 2 to guide your implementation.

    if strcmp(modType, 'psk')
        % --- Your response to Task 5 for BPSK/M-ary PSK modulation starts from below: ---
        
        if M == 2
            demoded_bits = real(rxSignal) > 0; % BPSK demodulation
        else
            demoded_bits = pskdemod(rxSignal, M, pi/M, 'gray'); % M-ary PSK demodulation with Gray coding
        end
        
        % --- Your response to Task 5 for PSK modulation ends here. ---
    elseif strcmp(modType, 'qam')
        % --- Your response to Task 5 for M-ary QAM starts from below: ---
        
        demoded_bits = qamdemod(rxSignal, M, 'gray'); % M-ary QAM demodulation with Gray coding
        
        % --- Your response to Task 5 for M-ary QAM ends here. ---
    end

    % Calculate bit errors
    [~, BER(i)] = biterr(bits, demoded_bits);

    %%% Plot eye diagram for the given Eb/No value
    if useFilter
        eyediagram(rxSignal(1:numSymbolsForEye*sps), 2*sps);
        title(['Eye Diagram: ' modType ' ' num2str(M) '-ary with ' codeType ' coding']);
    end
    
end
end
