
SNR = -4.55; %dB operating point
N = 1e8; % number of generated quantum states

sigma_square = 10^(-SNR/10);

random_numbers = randn(2,N);

Alice = random_numbers(1,:);

Bob = Alice + sqrt(sigma_square)*random_numbers(2,:);

%save the values to the files

filename = sprintf('Alice_%.2f_dB_sigma_square_%.3f_N_%i.txt',SNR,sigma_square,N);

fileID = fopen(filename, 'w');

% Write the vector to the file
fprintf(fileID, '%d ', Alice);

% Close the file
fclose(fileID);

filename = sprintf('Bob_%.2f_dB_sigma_square_%.3f_N_%i.txt',SNR,sigma_square,N);

fileID = fopen(filename, 'w');

% Write the vector to the file
fprintf(fileID, '%d ', Bob);

% Close the file
fclose(fileID);
