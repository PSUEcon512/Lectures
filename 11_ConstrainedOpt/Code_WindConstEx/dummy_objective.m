function [f, g] = dummy_objective(x0)
   %This is a dummy objective function which is used to let ktrlink serve
   %as an equation solver rather than an optimization routine.  We simply
   %supply the contraints and use this function as the "objective"...
   f = 0;
   g = zeros(size(x0));
end