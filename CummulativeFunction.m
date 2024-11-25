function cummulative_curve = cummulative_function(object, values, varargin)
if isa(object, 'double')
    cummulative_curve = zeros(1, length(values));
    for i = 1:length(values)
        cummulative_curve(i) = nnz(object <= values(i)) / length(object);
    end
elseif isa(object, 'function_handle')
    % If cummulative_curve() is used with a function input for object then
    % that function must be like a pdf (area under curve = 1) and an
    % additional argument must be provided to cummulative_curve()
    % specifying the lower bound on the function for which
    % function(x <= lower bound) = 0.
    if nargin == 2
        error("You must provide a lower bound specifying the point for which function(x <= lower bound) = 0");
    end
    zero = cell2mat(varargin(1));
    if nargin == 4
        number = cell2mat(varargin(2));
    else
        number = 1000;
    end
    cummulative_curve = zeros(1, length(values));
    for i = 1:length(values)
        x = linspace(zero, values(i), number);
        y = object(x);
        cummulative_curve(i) = NumericalIntegrator(x, y);
    end
end
end
