function MM = shift_matrix(M,s)

% MM = shift_matrix(M,s)
% shift values in matrix by s columns to the right

MM = nan*ones(size(M));

if s>=0,
    MM(:,s+1:end) = M(:,1:end-s);
else
MM(:,1:end+s) =      M(:,-s+1:end);
end
