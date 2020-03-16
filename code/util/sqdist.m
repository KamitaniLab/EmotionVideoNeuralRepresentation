function m = sqdist(p, q, A)
% SQDIST      Squared Euclidean or Mahalanobis distance.
% SQDIST(p,q)   returns m(i,j) = (p(:,i) - q(:,j))'*(p(:,i) - q(:,j)).
% SQDIST(p,q,A) returns m(i,j) = (p(:,i) - q(:,j))'*A*(p(:,i) - q(:,j)).

% Written by Tom Minka

[d, pn] = size(p);
[d, qn] = size(q);

if pn == 0 || qn == 0
  m = zeros(pn,qn);
  return
end

if nargin == 2
  pmag = sum(p .* p,1);
  qmag = sum(q .* q,1);
  m = repmat(qmag, pn, 1) + repmat(pmag', 1, qn) - 2*p'*q;
  %m = ones(pn,1)*qmag + pmag'*ones(1,qn) - 2*p'*q;
  
else
    if all(size(A)==size(p))
        Ap = A.*p;
        Aq = A.*q;
    else
        Ap = A*p;
        Aq = A*q;
    end
        pmag = sum(p .* Ap,1);
        qmag = sum(q .* Aq,1);
        m = repmat(qmag, pn, 1) + repmat(pmag', 1, qn) - 2*p'*Aq;
end
