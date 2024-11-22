function [times, peak_pos, bulk_pos, varargout] = PeakandBulkMotion(times, x_values, y_values, label1, varargin)
% 4 input arguments plots just the first data set of peak and mean position
% values against time.  The label must be given to force future clarity.
% If 5 arguments are given the 5th argument is assumed to be L so that the
% y plotting range is set to the full traitspace.
% If 8 arguments are given the last four arguments are assumed to be the
% data and label to be included in the plots showing the peak and mean
% position.
% If 9 arguments are given the last argument is assumed to be L as with 5
% inputs.
figure;
if nargin == 8 || nargin == 9
    times2 = varargin{1};
    x_values2 = varargin{2};
    y_values2 = varargin{3};
    label2 = varargin{4};
    peak_pos = zeros(1, length(times));
    bulk_pos = zeros(1, length(times));
    peak_pos2 = zeros(1, length(times2));
    bulk_pos2 = zeros(1, length(times2));
    for i=1:length(times)
        peak_posi = max(y_values(i, :)) == y_values(i, :);
        if sum(peak_posi) > 1
            peak_pos(i) = sum(x_values(peak_posi)) / sum(peak_posi);
        else
            peak_pos(i) = x_values(peak_posi);
        end
        bulk_pos(i) = sum(x_values .* y_values(i, :)) / sum(y_values(i, :));
    end
    for i=1:length(times2)
        peak_posi2 = max(y_values2(i, :)) == y_values2(i, :);
        if sum(peak_posi2) > 1
            peak_pos2(i) = sum(x_values2(peak_posi2)) / sum(peak_posi2);
        else
            peak_pos2(i) = x_values2(peak_posi2);
        end
        bulk_pos2(i) = sum(x_values2 .* y_values2(i, :)) / sum(y_values2(i, :));
    end
    plot(times, peak_pos, 'r-', times, bulk_pos, 'b--', times2, peak_pos2, 'm-', times2, bulk_pos2, 'k--');
    xlabel('Time');
    ylabel('Position');
    if nargin == 9
        ylim([0, varargin{5}]);
    end
    legend("Peak position (" + label1 + ")", "Mean position (" + label1 + ")", "Peak position (" + label2 + ")", "Mean position (" + label2 + ")", "Location", "best", Interpreter='latex');
    varargout{1} = times2;
    varargout{2} = peak_pos2;
    varargout{3} = bulk_pos2;
elseif nargin == 4 || nargin == 5
    peak_pos = zeros(1, length(times));
    bulk_pos = zeros(1, length(times));
    for i=1:length(times)
        peak_posi = max(y_values(i, :)) == y_values(i, :);
        if sum(peak_posi) > 1
            peak_pos(i) = sum(x_values(peak_posi)) / sum(peak_posi);
        else
            peak_pos(i) = x_values(peak_posi);
        end
        bulk_pos(i) = sum(x_values .* y_values(i, :)) / sum(y_values(i, :));
    end
    plot(times, peak_pos, 'b-', times, bulk_pos, 'r--', 'LineWidth', 2);
    xlabel('Time');
    ylabel('Position');
    if nargin == 5
        ylim([0, varargin{1}]);
    end
    if strcmp(label1, "ModeMean")
        legend("Modal trait", "Mean trait", Location='best', Interpreter='latex');
    elseif ~strcmp(label1, "")
        legend("Peak position (" + label1 + ")", "Mean position (" + label1 + ")", "Location", "best", Interpreter='latex');
    else
        legend("Peak position", "Mean position", "Location", "best", Interpreter='latex');
    end
else
    error("Error : incorrect number of arguments.");
end
end