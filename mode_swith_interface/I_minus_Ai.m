function val = I_minus_Ai(x)
    % class(x)
    x = x(:)';

    % Initialize the output vector to maintain vectorization
    [m,n]=size(x);

    gi = isa(x, 'gradient');

    val = GI_zeros(m,n,gi);

    if isa(x, 'gradient')
        val = [x,gradientinit(2)];
        val = val(m,n)';
        val = val(:)';
    end

    % Get indices based on the condition
    idx_near_zero = (inf(abs(x)) < 0.01); % Elements with absolute values less than 0.1
    idx_outside_zero = (~idx_near_zero); % Elements with absolute values greater than or equal to 0.1

     % Compute values for elements satisfying each condition
    if any(idx_near_zero)
        val(idx_near_zero) = GI_minus_Ai_near_zero(x(idx_near_zero),gi);
    end

    % if isa(x, 'gradient')
    %     disp('a')
    % end

    if any(idx_outside_zero)
        val(idx_outside_zero) = GI_minus_Ai_outside_zero(x(idx_outside_zero),gi);
    end

    val = val(:)';
end

function val = GI_minus_Ai_outside_zero(x,gi)
    gi = isa(x, 'gradient');
    % Compute Ai(-x)
    x = x(:); % Ensure x is a column vector

    % Compute common terms
    arg = GI_intval(2,gi)/GI_intval(3,gi) * x.^(GI_intval(3,gi)/GI_intval(2,gi)); % Argument for the Bessel function

    % Compute the Bessel function term (element-wise)
    bessel_term = I_besselj(GI_intval(1,gi)/3, arg) + I_besselj(-GI_intval(1,gi)/3, arg);

    % Compute the final value (element-wise)

    bessel_term = bessel_term(:);
    sqx_9 = sqrt(x / 9); sqx_9 = sqx_9(:);

    [m,n]=size(x);
    val = bessel_term.*sqx_9; % Removed transpose to match element-wise computation
    val = val(:)';

end


function val = GI_minus_Ai_near_zero(x,gi)
    gi = isa(x, 'gradient');
    % Compute Ai(-x)
    % x = x(:); % Ensure x is a column vector

    % Compute common terms
    Ai_0 = GI_intval(1,gi) / (3^(GI_intval(2,gi)/3) * gamma(GI_intval(2,gi)/3)); % Initial value at zero

    % Compute the error term (element-wise)
    error_term = abs(x).^4 ./ (3^(GI_intval(1,gi)/3) * GI_intval(12,gi) * gamma(GI_intval(1,gi)/3));
    % Compute the final value (element-wise)
    val = Ai_0 + error_term;
end



% function val = I_minus_Ai(x)
%     gi = isa(x, 'gradient');
%     % Compute Ai(-x)
%     x = x(:); % Ensure x is a column vector
% 
%     % Compute common terms
%     arg = GI_intval(2,gi)/GI_intval(3,gi) * x.^(GI_intval(3,gi)/GI_intval(2,gi)); % Argument for the Bessel function
% 
%     % Compute the Bessel function term (element-wise)
%     bessel_term = I_besselj(GI_intval(1,gi)/3, arg) + I_besselj(-GI_intval(1,gi)/3, arg);
% 
%     % Compute the final value (element-wise)
% 
%     bessel_term = bessel_term(:);
%     sqx_9 = sqrt(x / 9); sqx_9 = sqx_9(:);
% 
%     [m,n]=size(x);
%     val = bessel_term.*sqx_9; % Removed transpose to match element-wise computation
%     val = val(:)';
% end
% 
% 
% 
