function A=GI_ones(m,n,grad_or_intval)
   global INTERVAL_MODE;
   if INTERVAL_MODE
      if grad_or_intval
        A = gradientinit(I_ones(m,n));
      else
        A = I_intval(ones(m,n));
      end
   else
      A = ones(m,n);
   end
end



