function val = BLOM_EvalPolyBlock(A,C,X)

  A = A(any(C,1),:); % remove unused terms
  C = C(:,any(C,1));
  [Arows Acols Avals] = find(A);
  vX = X(Acols).^Avals; % powers of input variables
  expbool = (Avals == inf);
  vX(expbool) = exp(X(Acols(expbool))); % exponentials
  logbool =  (Avals == -inf);
  vX(logbool) = log(X(Acols(logbool))); % logarithms
  prods = ones(size(A,1),1);
  for v = 1:length(Avals)
      prods(Arows(v)) = prods(Arows(v)) * vX(v); % compute products for each term
  end
  val = C*prods; 
