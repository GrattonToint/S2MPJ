function varargout = your_obj_function( x, nargs )

D  = diag(1:length(x));
Dx = D*x;
varargout{1} = 0.5*x'*Dx;
if ( nargs > 1 )
   varargout{2} = Dx;
   if ( nargs > 2 )
      varargout{3} = D;
   end
end

return

end
