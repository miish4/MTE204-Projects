% circuit #1, matrix equation solver code 
clear; 
clc; 
A = [ -1 1 1 0 0;
    0 -1 0 1 1;
    5 0 15 0 0;
    0 10 -15 25 0;
    0 0 0 -25 20];

b = [0; 0; 80; 0; -50];

i = A\b;   % solves the system numerically

for j = 1:length(i) % display solutions as i1, i2 .. = 
    disp("i" + j + " = " + num2str(i(j))); 
end 

% circuit #4, matrix equation solver code 
A = [ 0 -1 1 1 0;
    1 0 1 0 -1;
    0 5 0 3 0;
    0 0 9 -3 6;
    7 -5 -9 0 0];

b = [0; 0; 24; 0; 0];

i = A\b;   % solves the system numerically

for v = 1:length(i)
    disp("i" + v + " = " + num2str(i(v)));
end 
