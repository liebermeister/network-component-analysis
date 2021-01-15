function [row_training,row_test,column_training,column_test] = make_crossval_groups(data,strategy)

% [row_training,row_test,column_training,column_test] = make_crossval_groups(data,strategy)
%
% definiere trainings-und testset fuer samples
% strategy 'random', 'half'
[n_rows,n_columns] = size(data);

switch strategy,

    case 'random',

row_training = randperm(n_rows);
row_test     = row_training(ceil(n_rows/2):end);
row_training = row_training(1:ceil(n_rows/2)-1);

column_training = randperm(n_columns);
column_test     = column_training(ceil(n_columns/2):end);
column_training = column_training(1:ceil(n_columns/2)-1);

case 'half',

% first and second half
row_training = 1:ceil(n_rows/2)-1;
row_test     = ceil(n_rows/2):n_rows;

column_training = 1:ceil(n_columns/2)-1;
column_test     = ceil(n_columns/2):n_columns;

end
