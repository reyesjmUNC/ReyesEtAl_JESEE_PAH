function h=logcolorbar(varargin)
% logcolorbar  - display a color bar, and back log transform the scale
%
% SYNTAX
%  h=logcolorbar
%
% NOTE
%  logcolorbar takes the same inputs as colorbar.m. See help colorbar
%  for more help on these inputs

h=colorbar(varargin{:});
yt=get(h,'YTick');
yt=exp(yt);
ytr=10.^floor(log10(yt))./100.*round(100*yt./(10.^floor(log10(yt))));
for k=1:length(ytr)
  str=sprintf('%.5g',ytr(k));
  ytl(k,1:length(str))=str;
end
set(h,'YTickLabel',ytl);  

