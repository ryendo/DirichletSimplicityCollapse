function val = func_left_hand_side_singularity(sigma,x)
    
    x = x(:)';
    gi = isa(x, 'gradient');
    c1 = GI_intval(1,gi);
    c2 = GI_intval(2,gi);
    c3 = GI_intval(3,gi);

    sigma = GI_intval(sigma,gi);
    
    val = (c2-sigma.^c3).^(c1/c3) * I_minus_Ai(((c2-sigma^c3).^(c2/c3)) .* x) .* I_d_minus_Ai(sigma.^c2 .* x) +sigma* I_minus_Ai(sigma.^c2 .* x) .* I_d_minus_Ai(((c2-sigma.^c3).^(c2/c3)) .* x);
    val = val(:)';
   
end