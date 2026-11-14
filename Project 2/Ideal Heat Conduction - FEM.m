% gauss-seidel w/ relaxation for project 2, question 2 ===
% solves A*x = b iteratively; w is the relaxation factor 
% w < 1 under-relaxation, w = 1 standard gauss-seidel, w > 1 over-relaxation
clear; 
clc;

es = 0.05;          % stopping criteria: max approx. relative error (%)
maxIter = 200;      % iteration limit

% matrix setup   

% A - global stiffness matrix 
A = [10 -10 0 0 0;  
    -10 20 -10 0 0;
    0 -10 20 -10 0;
    0 0 -10 20 -10;
    0 0 0 -10 20];

% b - source term vector
b = [-140; 120; 120; 120; 220];

[A, b] = pivot_rows(A, b, 'zero'); % swap only if a diagonal entry is zero
T0 = 2 * ones(length(b), 1);       % initial guess
ws = [0.5 0.8 1.0 1.02 1.1];       % relaxation factors to compare 

run_case('Question 2', A, b, T0, ws, es, maxIter);

% local functions - - -

% runs gauss-seidel for each w and prints a summary table 
function run_case(name, A, b, T0, ws, es, maxIter)
    fprintf('\n==== %s ====\n', name);
    fprintf(' w | iterations | converged | final error (%%)\n');
    fprintf('--------------------\n');

    for j = 1:length(ws)
        w = ws(j);
        [T, iters, conv, errHist] = gauss_seidel_relax(A, b, T0, w, es, maxIter);
        fprintf('%5.2f | %4d | %d | %.3g\n', w, iters, conv, errHist(end));
        
        % only print solution if this w converged 
        if conv 
           fprintf('   T = [%s]\n', num2str(T', '%.4f   '));

           T5 = T(5);   % extract T5, the 5th element of T 
           T6 = 10;     % given/known boundary value @ node 6 
            
           % pass T5 into gradient calculating function (refer to bottom) 
           dT_dx_end = compute_right_gradient(T5, T6);
           fprintf(' Calculated gradient at right boundary (dT/dx  at x=0.5): %.4f C/m\n', dT_dx_end);
        end 
    end
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
        error('zero remains on the diagonal after pivoting.')
    end
end

function dT_dx_6 = compute_right_gradient(T5, T6)
    % from element 5 / row 6 analysis:
    % 10 * (-T5 + T6) = (dT/dx)_end + 60
    % rearranging to solve for (dT/dx)_end: 
    
    dx = 0.1; 
    f = 1200; 
    source_term_half = (f * dx) / 2; 
    dT_dx_6 = 10 * (T6 - T5) - source_term_half;
end
