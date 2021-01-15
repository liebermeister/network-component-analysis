[A,B,Amin,Amax,Bmin,Bmax] = nca_bootstrap(X,Wblock);

nit = 100;

[A,B] = network_component_analysis(X,Wblock);

X_rec = A*B;

errors = X(:)-X_rec(:);

for it = 1: nit,
  X_boot = X_rec + reshape( errors(ceil(length(errors)*rand(length(errors),1))), size(X,1),size(X,2));
  [A{it},B{it}] = network_component_analysis(X_boot,Wblock);
end
