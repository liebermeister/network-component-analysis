function cd_data(filename,type);

% cd_data(filename,type);

global data_directory

switch type,
  case 'network',
    eval(['cd ' data_directory '/genetic_networks/']);
  case 'expression',
    eval(['cd ' data_directory '/expression_data/' filename]);
  case 'analysis',
    eval(['cd ' data_directory '/analysis/' filename]);
  case 'graphics',
    eval(['cd ' data_directory '/graphics/' filename]);
  case 'graphics',
    eval(['cd ' data_directory '/graphics/' filename]);
end