function a = my_lsqlin(B,x,s,a_guess)

% don't forget to addpath ~/matlab_packages/minq

boundary = 50;

epsilon = 10^-10;
au = -boundary * (s<=0) + epsilon;
ao =  boundary * (s>=0) - epsilon;

a_guess(a_guess<au)=0;
a_guess(a_guess>ao)=0;

a = minq(x'*x,-2*B'*x,2*B'*B,au,ao,0,a_guess);

a(a<au)=0;
a(a>ao)=0;

%if max(abs(a)==boundary), fprintf('my_lsqlin: boundary  reached .. '); end

if norm(B * a_guess - x) < norm(B * a - x),
  a = a_guess;
end