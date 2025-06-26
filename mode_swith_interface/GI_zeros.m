function A=GI_zeros(m,n,grad_or_intval)
   global INTERVAL_MODE;
   if INTERVAL_MODE
      if grad_or_intval
        A = gradientinit(I_zeros(m,n));
      else
        A = I_intval(zeros(m,n));
      end
   else
      A = zeros(m,n);
   end
end



