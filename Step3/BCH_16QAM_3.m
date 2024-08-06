tic; % 开始计时

img = imread('C:\Users\earth\OneDrive\UNSW Course\2024-T2\TELE4652 Mobile & Satellite Comm System\4652 LAB\Experiment 4A & 4B\Lab 4 Package for Students\1.png'); % 读取图像

% 将图像转换为一维字节流
img_bytes = typecast(img(:), 'uint8');

% 将字节流转换为位流
% 预分配数组
bit_stream = zeros(1, length(img_bytes) * 8);

% 将字节流转换为位流
for i = 1:length(img_bytes)
    byte = img_bytes(i);
    bits = bitget(byte, 8:-1:1); % 将字节转换为二进制位（高位在前）
    bit_stream((i-1)*8+1:i*8) = bits;
end

% 定义BCH码参数
m = 3;          % 生成多项式的次数
n = 2^m - 1;    % 码长
k = 4;          % 数据长度
t = 1;          % 可纠正的错误数

% 确保 bit_stream 是列向量
if isrow(bit_stream)
    bit_stream = bit_stream.';  % 转置为列向量
end

% 创建BCH码对象
bchEncoder = comm.BCHEncoder(n, k);
bchDecoder = comm.BCHDecoder(n, k);

% 编码
code_bits = step(bchEncoder, bit_stream);

% 确保比特流的长度是 num_bits_per_symbol 的整数倍
if mod(length(code_bits), log2(M)) ~= 0
    error('The number of bits is not compatible with the symbol size.');
end

% 使用16PSK
% 将比特流分为4比特一组进行16PSK调制
M = 16; % 16PSK的M值
bit_grouping = reshape(code_bits, log2(M), []).'; % 重组数组以匹配每个符号4比特

% 将二进制比特组转换为十进制整数
symbols = bi2de(bit_grouping, 'left-msb');  % 'left-msb' 表示最左边的比特是最高位

% 使用16PSK调制
mod_bits = qammod(symbols, M, 'UnitAveragePower', true);

% 定义信噪比 (SNR)
EbNo_dB = 15; % Eb/No values in dB
sps = 1; % 每符号的采样数
M = 16; % 调制的符号数
SNR = EbNo_dB + 10*log10(log2(M)) - 10*log10(sps);

% 通过AWGN信道传输
rxSignal = awgn(mod_bits, SNR, 'measured');

% 解调和解码
received_grouping = qamdemod(rxSignal, M, 'UnitAveragePower', true);

% 将十进制整数转换为二进制比特矩阵
bit_matrix = de2bi(received_grouping, log2(M), 'left-msb');

% 将二维比特矩阵重构为一维比特流
received_bits = reshape(bit_matrix.', [], 1);  % 使用转置以确保比特顺序正确

% 解码
decoded_bits = step(bchDecoder, received_bits);

% 计算误码数
numErrors = sum(bit_stream ~= decoded_bits);

% 计算误码率
ber = numErrors / length(bit_stream);
disp(['误码率 (BER): ', num2str(ber)]);

% 将位流转换回字节流并重组图像
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

% 显示原图像
figure;
subplot(1, 2, 1);
imshow(img);
title('原图像');

% 显示重组的图像
subplot(1, 2, 2);
imshow(img_reconstructed);
title('重组后的图像');

% 显示星座图
figure;
subplot(1,2,1);
plot(real(mod_bits), imag(mod_bits), 'bo');
title('调制信号的星座图');
xlabel('实部');
ylabel('虚部');
axis square;
grid on;

subplot(1,2,2);
plot(real(rxSignal), imag(rxSignal), 'ro');
title('接收信号的星座图');
xlabel('实部');
ylabel('虚部');
axis square;
grid on;
toc; % 结束计时并显示耗时