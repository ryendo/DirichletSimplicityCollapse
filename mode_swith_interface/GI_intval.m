function output = GI_intval(var,grad_or_intval)
   % switchable wrapper for gradient and intval
   global INTERVAL_MODE;
   if INTERVAL_MODE
      if grad_or_intval && ~isa(var, 'gradient')
        output = gradientinit(var);
      else
        output = intval(var);
      end
   else
      if ischar(var)
        output = str2double(var);
      else
        output = var;
      end
   end
end






