function [prev_data,now_data] = water_dispatch_history(prev_data,now_data,dispatch,subnet)
if ~isempty(dispatch) && isfield(subnet,'Hydro')
    xi = max(1,nnz((dispatch.Timestamp-1e-6)<prev_data.Timestamp(1) & dispatch.Timestamp>0));
    xf = nnz((dispatch.Timestamp-1e-6)<=now_data.Timestamp & dispatch.Timestamp>0);
    yi = 1 + nnz(prev_data.Timestamp<dispatch.Timestamp(1));
    for n = 1:1:length(subnet.Hydro.lineNumber)
        prev_data.Hydro.OutFlow(yi:end,n) = interp1(dispatch.Timestamp(xi:xf),dispatch.LineFlows(xi:xf,subnet.Hydro.lineNumber(n)),prev_data.Timestamp(yi:end));
        now_data.Hydro.OutFlow(n) = dispatch.LineFlows(xf,subnet.Hydro.lineNumber(n));
    end
end
end%Ends function water_dispatch_history