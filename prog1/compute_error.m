function [err] = compute_error(analytical, appx, gp)
  err = 0;
  for i=1:length(gp)-1
    piece = @(x) ((appx(i+1)-appx(i)) ./ (gp(i+1)-gp(i))) .* (x-gp(i)) + appx(i);
    err += integral(@(x) abs(analytical(x) - piece(x)) .^ 2, gp(i),gp(i+1));
  end 
  err = sqrt(err);
end
