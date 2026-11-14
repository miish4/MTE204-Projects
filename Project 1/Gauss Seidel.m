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


% gauss-seidel w/ relaxation for circuit #1 & #4 ===

% solves A*x = b iteratively; w is the relaxation factor 
% w < 1 under-relaxation, w = 1 standard gauss-seidel, w > 1 over-relaxation
clear; 
clc;

es = 0.05;          % stopping criteria: max approx. relative error (%)
maxIter = 200;      % iteration limit

% circuit #1 - - -  
A = [ -1 1 1 0 0;
    0 -1 0 1 1;
    5 0 15 0 0;
    0 10 -15 25 0;
    0 0 0 -25 20];

b = [0; 0; 80; 0; -50];

[A, b] = pivot_rows(A, b, 'zero'); % swap only if a diagonal entry is zero
x0 = 2 * ones(length(b), 1);       % initial guess
ws = [0.5 0.8 1.0 1.02 1.1];       % relaxation factors to compare 

run_case('Circuit #1', A, b, x0, ws, es, maxIter);

% circuit #4 - - - 
A = [ 0 -1 1 1 0;
    1 0 1 0 -1;
    0 5 0 3 0;
    0 0 9 -3 6;
    7 -5 -9 0 0];

b = [0; 0; 24; 0; 0];

[A, b] = pivot_rows(A, b, 'max');   % zeroes on diagonal, put largest entry in each column on diagonal
x0 = 2 * ones(length(b), 1);        % initial guess
ws = [0.5 0.7 0.8 0.9 1.0 1.1];     %% relaxation factors to compare

run_case('Circuit #4', A, b, x0, ws, es, maxIter);

% local functions - - - 

% runs gauss-seidel for each w and prints a summary table 
function run_case(name, A, b, x0, ws, es, maxIter)
    fprintf('\n==== %s ====\n', name);
    fprintf(' w | iterations | converged | final error (%%)\n');
    fprintf('--------------------\n');

    for j = 1:length(ws)
        w = ws(j);
        [x, iters, conv, errHist] = gauss_seidel_relax(A, b, x0, w, es, maxIter);
        fprintf('%5.2f | %4d | %d | %.3g\n', w, iters, conv, errHist(end));
        
        % only print solution if this w converged 
        if conv 
           fprintf('   x = [%s]\n', num2str(x', '%.4f   '));
        end 
    end
    
    % compare with matlab's direct solver 
    fprintf('\n %s MATLAB direct solution (A\\b):\n', name);
    disp((A\b)');

end 

% gauss-seidel sweep w/ relaxation weighting w 
function [x, iters, converged, errHist] = gauss_seidel_relax(A, b, x0, w, es, maxIter)
    n = length(b);
    x = x0;
    converged = false; 
    errHist = zeros(maxIter, 1);        % approx. relative error at each iteration 

    for iters = 1:maxIter 
        x_old = x;                      % save previous iteration

        % one sweep: solve eq i for x(i), using newest values available 
        for i = 1:n
            % off-diagonal terms (left: new values, right: old values)
            sigma = A(i,1:i-1) * x(1:i-1) + A(i, i+1:n) * x(i+1:n);
            x_gs = (b(i) - sigma) / A(i, i);            % standard gauss-seidel value 
            x(i) = w * x_gs + (1-w) * x_old(i);         % relaxed (weighted) update
        end

        % approx. relative error (%): largest change over all unknowns 
        ea = max(abs((x-x_old) ./ (x + eps))) * 100; 
        errHist(iters) = ea;

        if ea < es                      % converged
            converged = true;
            break;
        end
        if ~isfinite(ea) || ea > 1e10   % diverging, stop early 
            break;
        end 

        errHist = errHist(1:iters);     % drop unused entries (covers every exit path) 
    end
end 

% row swapping (partial pivoting) to get rid of zeroes on the diagonal 
% 'zero' : swap only when A(k,k) is zero
% 'max' : always move largest |entry| of column k onto the diagonal 
function [A, b] = pivot_rows(A, b, mode)
    n = length(b);
    for k = 1:n
        % in 'zero' mode, skip columns whose diagonal is already nonzero 
        if strcmp(mode, 'zero') && abs(A(k,k)) > eps
            continue;
        end
        
        [~, p] = max(abs(A(k:n, k)));       % row w/ largest entry in column k 
        p = p + k - 1;                      % convert to full-matrix row index
        
        if p ~= k
            A([k, p], :) = A([p k], :);     % swap rows of A 
            b([k, p]) = b([p, k]);          % swap same rows of b 
            fprintf('Pivoting: swapped rows %d and %d\n', k, p);
        end
    end
    
    % safety check: stop if pivoting could not fix the diagonal 
    if any(abs(diag(A)) < eps)
        error('zero remains on the diagonal after pivoting.');
    end
end
