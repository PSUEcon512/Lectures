function nodes = haltonNormShuffle(nPoints, nDims, seed)
% haltonNormShuffle(nPoints, nDims) returns a nDimensional set of nodes
% derived from shuffled halton sequence draws that are now normally
% distributed. This code calls haltonseq and makes two changes to the
% draws:
%
% 1 The draws are shuffled dimension by dimension to avoid colinearity in
% draws when using a large number of prime numbers
%
% 2 Standard halton draws are in quasi-MC for Uniform(0,1) we achieve
% an approximation for the standard normal distribution by passing these
% draws through the inverse normal cdf. 
%
  
   assert(nargin>=2 & nargin <= 4);
   if nargin==2
       seed = 101;
   end

   %Note that we transpose here so that each column is a draw and the
   %number of rows is the number of dimensions.
   UnifStraight = haltonseq(nPoints, nDims)';
   
   %Now we shuffle each dimension we will be careful to replace the seed.
   sprev = rng();
   myrng = rng(seed,'twister');
   UnifShuffle = zeros(nDims,nPoints);
   for row = 1:nDims
       shuffle = randperm(nPoints);
       UnifShuffle(row,:) = UnifStraight(row,shuffle);
   end
   
   %Replace the seed
   rng(sprev)
   
   %Now convert to Standard normal
   nodes = norminv(UnifShuffle);
end