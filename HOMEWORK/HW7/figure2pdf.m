function figure2pdf(filename,openfile)
%FIGURE2PDF Exports the current MATLAB figure to a PDF file.
%
%   This function is useful for preparing high-quality graphic files
%   for research papers.
%
%   Syntax:
%      FIGURE2PDF(filename)  ........ Export figure only.
%      FIGURE2PDF(filename,true)  ... Export figure and open PDF file.
%
%   Example:
%      [x,y] = meshgrid(-2:0.1:2);
%      z = exp(-x.^2-y.^2);
%      surf(x,y,z)
%      xlabel('$x$','Interpreter','LaTeX','FontSize',11)
%      ylabel('$y$','Interpreter','LaTeX','FontSize',11)
%      zlabel('$z=f(x,y)$','Interpreter','LaTeX','FontSize',11)
%      set(gca,'TickLabelInterpreter','LaTeX','FontSize',11)
%      box on
%      figure2pdf('bell.pdf',1)
%
%   Author:
%      Ildeberto de los Santos Ruiz
%      idelossantos@ittg.edu.mx
%
% See also PRINT, SAVEFIG.

f = gcf;
a = gca;

set(a,'LooseInset',get(a,'TightInset'))
set(f,'Units','inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','auto',...
    'PaperUnits','inches',...
    'PaperPosition',[0,0,pos(3),pos(4)],...
    'PaperSize',[pos(3), pos(4)])

if nargin < 1
    filename = 'figure.pdf';
end

% if ~strcmpi(filename(end-3:end),'.pdf')
%     filename = [filename,'.pdf'];
% end
print(f,filename,'-painters','-dpdf','-r0')
if nargin > 1 && openfile == true
    open(filename)
end