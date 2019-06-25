function k = f_index(a,b)
%% find index of b where value of a matches b
%% a can be a number or string, b is a vector or cell array
%% a can be a list that gets matched to the corresponding index (returns a vector)
n = length(a);
if isnumeric(a)
    if n == 1
        k = nonzeros((1:length(b))'.*(a==b));
    else
        k = zeros(n,1);
        for i = 1:1:n
            k_i = nonzeros((1:length(b))'.*(a(i)==b));
            if ~isempty(k_i)
                k(i) = k_i;
            end
        end
    end
elseif iscell(a) && n>1
    k = [];
    for i = 1:1:length(a)
        k = [k;nonzeros((1:length(b))'.*strcmpi(a{i},b))];
    end
else
    k = nonzeros((1:length(b))'.*strcmpi(a,b));
end
end%Ends function f_index