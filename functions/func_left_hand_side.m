function val = func_left_hand_side(s,x)
    
    x = x(:)';
    gi = isa(x, 'gradient');
    c1 = GI_intval(1,gi);
    c2 = GI_intval(2,gi);
    c3 = GI_intval(3,gi);
    
    cs = GI_intval(s,gi);

    val = (c1 + cs).^(c1/c3) * I_minus_Ai(((c1 + cs).^(c2/c3)) .* x) .* I_d_minus_Ai(((c1 - cs).^(c2/c3)) .* x) +(c1 - cs)^(c1/c3) * I_minus_Ai(((c1 - cs).^(c2/c3)) .* x) .* I_d_minus_Ai(((c1 + cs).^(c2/c3)) .* x);
    
    val = val(:)';
end


