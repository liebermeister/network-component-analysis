function nca_dependencies()

% nca_dependencies()
%
% Check whether the packages used by the NCA toolbox are installed

if ~exist('matlab_utils_version','file')
  error('Please install the Matlab utility functions (https://github.com/liebermeister/matlab-utils) - They are required by the NCA toolbox.');
end

package_version_check('The Matlab utility functions', 'matlab_utils_version', '1.0.0');
