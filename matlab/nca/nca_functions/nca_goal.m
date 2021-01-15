function [goal,residual,prior] = nca_goal(X,A,B,general,regularise_flag,lambda_dep,lambda_A,lambda_B)

residual = nanmean(nanmean((X-A*B).^2)');
goal     = residual;

if regularise_flag,
  if length(lambda_dep) == 1,
    Ared = A; Ared(:,general)=0;
    Bred = B; Bred(general,:)=0;
    goal = residual + lambda_dep^2 * nanmean(nanmean((Ared*Bred).^2)');
  else,
    goal = residual + nanmean(nanmean( (A * diag(lambda_dep) * B).^2)' );
  end
end

prior = 0;

if lambda_A ~= 0, prior = prior + lambda_A^2 * nanmean(nanmean(A.^2)); end
if lambda_B ~= 0, prior = prior + lambda_B^2 * nanmean(nanmean(B.^2)); end

goal = goal + prior;