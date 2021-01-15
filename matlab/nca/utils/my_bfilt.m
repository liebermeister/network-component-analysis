function r = my_bfilt(data)

% r = my_bfilt(data)

if length(find(isfinite(data))),
  data =  fill_nans_by_polynomial_interpolation(data,(1:size(data,2))');
end

if size(data,2)>=5,
  
  for it =1:size(data,1),
    
    this_data = data(it,:);
    dum2      = [];  
    
    for it2 =3:length(this_data)-2,
      dum = this_data(it2-2:it2+2);
      [d1,order] = sort(dum);
      if d1(1)==d1(end),
        dum2 = [dum2; d1(1:3)];
      else,
        dum([order(1),order(end)])  = nan;
        dum2 = [dum2; dum(isfinite(dum))];
      end
    end
    r(it,:) = (dum2*[0.25 0.5 0.25]')';
  end
  
  r = [data(:,1:2)*[0.7,0.3]',  data(:,1:3)*[0.25 0.5 0.25]', ...
       r ,...
       data(:,end-2:end)*[0.25 0.5 0.25]', data(:,end-1:end)*[0.3,0.7]'];
else,
  r = data;
end