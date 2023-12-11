function Ans = TDMA(W,P,E,Q)
%W = west boundary: 1D array, len = n, W(1) = 0 
%P = Middle node: 1D array, len = n, all index are defined
%E = East boundary: 1D array, len = n, E(n) = 0
%Q = Answer: 1D array, len = n, all index are defined
%X = Solution of the tri-diognal matrix.

n = length(Q);
X = zeros(1,n);
for i = 2:n                                %for loop is extra. It might be removed. But it is easier to understand in this way. 
    P(i) = P(i) - E(i-1) * W(i) / P(i-1);
    Q(i) = Q(i) - Q(i-1) * W(i) / P(i-1);
end
X(end) = Q(end)/P(end);
for i = n-1:-1:1
    X(i) = (Q(i) - E(i) * X(i+1))/P(i);
end
Ans = X;
end
