function Gfb_plot(figure_number, axis_vector, title_string, ...
		  xaxis_label, yaxis_label, domain, value)
% function Gfb_plot(figure_number, axis_vector, title_string, ...
%                   xaxis_labal, yaxis_label, domain, value)
%
% creates plots with labeled axes and defined axis limits
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Nov 2006

% filename : Gfb_plot.m

    figure(figure_number);
    plot(domain, value);
    axis(axis_vector);
    title(title_string);
    xlabel(xaxis_label);
    ylabel(yaxis_label);

