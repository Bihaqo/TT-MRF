function plot_comparison(temperatureArray, logZ_TT, logZ_JT, logZ_BP, logZ_MF, logZ_TREEEP, logZ_AIS)
	logZ_true = logZ_JT;

	h = figure;
	hold all;
	p_arr = [];
	legends = {};

	p_arr(end + 1) = plot_method(temperatureArray, logZ_BP - logZ_true, '-s', [0 0 0], 1);
	legends{end + 1} = 'BP';

	p_arr(end + 1) = plot_method(temperatureArray, logZ_MF - logZ_true, '-v', [0 0 0], 1);
	legends{end + 1} = 'MF';

	p_arr(end + 1) = plot_method(temperatureArray, logZ_TREEEP - logZ_true, 'o-', [ 0 0 0], 1);
	legends{end + 1} = 'TREEEP';

	p_arr(end + 1) = plot_method(temperatureArray, logZ_AIS - logZ_true, '-^', [0 0 0], 1);
	legends{end + 1} = 'MCMC - AIS';

	p_arr(end + 1) = plot_method(temperatureArray, logZ_TT - logZ_true, '-', [1 0 0], 1.5);
	legends{end + 1} = 'TT';

	hLegend = legend(p_arr, legends);

	hXLabel = xlabel('temperature T');
	hYLabel = ylabel('$|\log \hat{Z} - \log Z|$', 'Interpreter','latex');
	set(gca,'YScale','log','YMinorTick','on',...
    'YColor',[0.3 0.3 0.3],...
    'XScale','log',...
    'XMinorTick','on',...
    'XColor',[0.3 0.3 0.3],...
    'TickDir','out',...
    'TickLength',[0.02 0.02],...
    'Position',[0.095360824742268 0.138888888888889 0.868556701030928 0.845238095238095],...
    'FontName','Times New Roman');
	% save_figure(h, file_name);
end


function p = plot_method(x_vals, results, spec, color, width)
	y = abs(results);

	[quartiles] = quantile(y, [0.25, 0.5, 0.75]);
	error_down = quartiles(2, :) - quartiles(1, :);
	mean_values = quartiles(2, :);
	error_up = quartiles(3, :) - quartiles(2, :);
	p = errorbar(x_vals, mean_values, error_down, error_up, spec);

	% Line color.
	set(p, 'Color', color);

	hb = get(p, 'children');
	% Line width.
	set(hb(1), 'Linewidth', width);
	% Error bar line width.
	set(hb(2), 'Linewidth', 0.5);

	set(hb(1),'MarkerSize', 5);
	set(hb(1),'MarkerFaceColor', [0 0 0]);
	set(hb(1),'MarkerEdgeColor', [0 0 0]);
end
