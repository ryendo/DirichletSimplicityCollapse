function output=ID_hull(I1,I2,di)
   global INTERVAL_MODE;
   if INTERVAL_MODE
      if ~di
        output = hull(I1,I2);
      else
        output = (I1+I2)/2;
      end
   else
      output = (I1+I2)/2;
   end
end



