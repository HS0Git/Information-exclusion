%T:the number of probability distributions
%d:length of distribution
%c:the sum of indexes of coincidence of T distributions
function [entropyr]=SentropyIC(c,T)
      L=ceil(T./c);                     % maximal number of nonzero probabilities in a dsitribution
      k=floor(c.*L.*(L-1)-T.*(L-1));    % number of uniform distributions over L-1 outcomes 
      temp=c-k./(L-1+eps)-(T-k-1)./L;   % index of coincidence of the nonuniform distribution
      p1=(1+sqrt((L.*temp-1)./(L-1+eps)))./L+eps;p2=(1-sqrt((L.*temp-1).*(L-1+eps)))./L+eps;
      entropyr=k.*log2(L-1+eps)+(T-k-1).*log2(L)-(L-1).*p1.*log2(p1)-p2.*log2(p2);
end