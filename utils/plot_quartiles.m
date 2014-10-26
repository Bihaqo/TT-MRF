function p = plot_quartiles(x, y, spec, lineWidth)
	% Plot median and upper and lower quartiles as error bars. Returns the plot
	% handler.
	% 
    % Input:
    %   x         -- vector 1 x N with X values
    %   y         -- matrix M x N with all Y values.
    % 
    % Optional input:
    %   spec      -- character string which defines line specifications (see help plot)
    %   lineWidth -- number
	% 

	[quartiles] = quantile(y, [0.25, 0.5, 0.75]);
	errorDown = quartiles(2, :) - quartiles(1, :);
	median = quartiles(2, :);
	errorUp = quartiles(3, :) - quartiles(2, :);
	if nargin >= 3
		p = errorbar(x, median, errorDown, errorUp, spec);
	else
		p = errorbar(x, median, errorDown, errorUp);
	end

	hb = get(p, 'children');
	if nargin >= 4
		set(hb(1), 'Linewidth', lineWidth);
	end
	% Error bar line width.
	set(hb(2), 'Linewidth', 0.5);

	set(hb(1),'MarkerSize', 5);
	set(hb(1),'MarkerFaceColor', [0 0 0]);
	set(hb(1),'MarkerEdgeColor', [0 0 0]);
end
