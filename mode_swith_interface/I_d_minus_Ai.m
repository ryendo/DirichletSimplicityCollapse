function val = I_d_minus_Ai(x)

    % Initialize the output vector to maintain vectorization
    x = x(:)';
    [m,n]=size(x);
    
    gi = isa(x, 'gradient');

    val = GI_zeros(m,n,gi);

    if isa(x, 'gradient')
        val = [x,gradientinit(2)];
        val = val(m,n)';
        val = val(:)';
    end

    % Get indices based on the condition
    idx_near_zero = inf(abs(x)) < 0.01; % Elements with absolute values less than 0.1
    idx_outside_zero = ~idx_near_zero; % Elements with absolute values greater than or equal to 0.1

    if any(idx_near_zero)
        val(idx_near_zero) = GI_d_minus_Ai_near_zero(x(idx_near_zero),gi);
    end

    if any(idx_outside_zero)
        val(idx_outside_zero) = GI_d_minus_Ai_outside_zero(x(idx_outside_zero),gi);
    end

    val=val(:)';
end

function val = GI_d_minus_Ai_outside_zero(x,gi)
    gi = isa(x, 'gradient');

    % Compute Ai(-x)
    x = x(:);

    % Compute common terms
    % Element-wise square root
    x_power = x.^(GI_intval(3,gi)/GI_intval(2,gi));                  % Element-wise power
    arg = (GI_intval(2,gi) / GI_intval(3,gi)) * x_power;            % Argument for the Bessel function

    % Compute Bessel function terms
    term1_bessel = I_besselj(GI_intval(2,gi)/GI_intval(3,gi), arg);
    term2_bessel = I_besselj(GI_intval(4,gi)/GI_intval(3,gi), arg);
    term3_bessel = I_besselj(GI_intval(1,gi)/GI_intval(3,gi), arg);

    % Compute individual terms
    term1 = (GI_intval(1,gi)/GI_intval(3,gi) * x)' .* term1_bessel;
    term2 = (GI_intval(1,gi)/GI_intval(3,gi) * x)' .* term2_bessel;
    term3 = (GI_intval(1,gi)/GI_intval(3,gi) ./ sqrt(x))' .* term3_bessel;

    % Compute the final value
    val = -term1 + -term2 + term3;

end

function val = GI_d_minus_Ai_near_zero(x,gi)
    gi = isa(x, 'gradient');

    % Compute Ai(-x)
    x = x(:);
    % Compute common terms

    Ai_0 = 1/(3^(GI_intval(2,gi)/GI_intval(3,gi))*gamma(GI_intval(2,gi)/GI_intval(3,gi)));
    Aid_0 = 1/(3^(GI_intval(1,gi)/GI_intval(3,gi))*gamma(GI_intval(1,gi)/GI_intval(3,gi)));

    ep = max(abs(x));

    % Compute the Bessel function term
    error_term = ep^GI_intval(2,gi)/GI_intval(2,gi)*Ai_0 + ep^GI_intval(4,gi)/(GI_intval(3,gi)^(GI_intval(1,gi)/GI_intval(3,gi))*12*gamma(GI_intval(1,gi)/GI_intval(3,gi)));

    % Compute the final value
    val = Aid_0 + error_term;
end
