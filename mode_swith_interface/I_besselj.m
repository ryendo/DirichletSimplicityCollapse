function output = I_besselj(n,x)
   global INTERVAL_MODE;
   if INTERVAL_MODE
      output = intval_besselj(n,x)';
   else
      % Use MATLAB's built-in Bessel function
      output = besselj(n,x);
   end
end

function val = intval_besselj(n,x)

    N = 30;       % Number of terms in the series expansion

    x = x(:);     % Convert x to a column vector

    % Vector of indices k
    k = 0:N; % Row vector of size 1 x (N+1)

    % Compute coefficients
    numerator = (-1).^k; % 1 x (N+1)
    denominator = factorial(k) .* gamma(n + k + 1); % 1 x (N+1)
    coefficients = numerator ./ denominator; % 1 x (N+1)

    % Compute exponents
    exponents = 2 * k + n; % 1 x (N+1)

    % Compute powers of (x/2)
    % powers = (x/2) .^ exponents; % M x (N+1)
    a = x/2; b = exponents;
    A = repmat(a, 1, length(b));
    B = repmat(b, length(a), 1);
    powers = A .^ B;

    % % Compute a_k

    S_N = powers* coefficients';

    % Compute the error term
    error_term = (1./(factorial(N+1).*gamma(n+N+2))).*(x/2).^(2*(N+1)+n); % M x 1

    % Create an interval including the error
    val = S_N + error_term; % Interval of size M x 1

    % Return the output as a column vector
    val = val(:);
end
