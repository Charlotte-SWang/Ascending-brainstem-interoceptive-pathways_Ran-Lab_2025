function weighted_rose(varargin)

% Weighted rose plot
%
% weighted_rose(theta, w) plots the histogram of radians defined in theta
% multiplied by the corresponding weight of the radians defined in array w.
% The default bin number is 20.
%
% weighted_rose(theta, w, nbin) works the same with additionally setting
% the number of bins defined in nbin.
% theta and w must be the same size
%
%
% This script is ideal for example plotting the phase amplitude coupling
% distribution
%
%
% Please cite this script as follows (or other ways described in https://blogs.mathworks.com/community/2010/12/13/citing-file-exchange-submissions/)
% if you intend to use it in a research method:  Puszta, András (2017) Weighted rose plot,
% (http://www.mathworks.com/matlabcentral/fileexchange/), MATLAB Central File Exchange. 
%
%

[cax,args,nargs] = axescheck(varargin{:});
if nargs < 2
    error(message('MATLAB:narginchk:notEnoughInputs'));
elseif nargs > 3
    error(message('MATLAB:narginchk:tooManyInputs'));
end

theta=squeeze(args{1});
w=squeeze(args{2});

if size(theta)~=size(w) 
  error('theta and w must be the same size!');
end

if nargs>2
    nbins=args{3};
else
    nbins=20;
end

angles=linspace(0,2*pi, nbins+1);
edges=linspace(0,2*pi, nbins+1);
rho=zeros(1, nbins);

for a=1:nbins
    
        rho(a)=sum(w(find(theta>=angles(a)&theta<angles(a+1))));
    
    
end

% Form radius values for histogram triangle
nn = rho(:); 
[m,n] = size(nn);
mm = 4*m;
r = zeros(mm,n);
r(2:4:mm,:) = nn;
r(3:4:mm,:) = nn;

% Form theta values for histogram triangle from triangle centers (xx)
zz = edges;

t = zeros(mm,1);
t(2:4:mm) = zz(1:m);
t(3:4:mm) = zz(2:m+1);

    h = polar(t,r);
    

