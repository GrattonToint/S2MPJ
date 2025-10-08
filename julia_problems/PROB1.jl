function PROB1( action, args... )

if action == "setup"
   name = "PROB1"
   call = eval( Meta.parse( name ) )
   return call
elseif action == "act"
   return 42
else
   return args[1]("act")
end

end
