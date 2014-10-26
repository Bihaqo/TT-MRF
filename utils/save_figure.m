function save_figure(figureHandler, filename)
	% Store figure in .fig, .png, .eps and tikz formats
	% 

	set(gcf, 'PaperPositionMode', 'auto');
	saveas(figureHandler, [filename, '.fig'], 'fig');
	print('-depsc2', [filename, '.eps']);
	print('-dpng', [filename, '.png']);
	try
		matlab2tikz('figurehandle', figureHandler, 'filename', [filename, '.tikz']);
	catch err
		warning('Cannot save figure in the Tikz format, make sure that matlab2tikz library is installed.');
	end
end
