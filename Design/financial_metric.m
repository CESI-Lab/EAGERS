function [irr,payback] = financial_metric(mC,mCbaseline,rate)
% Internal Rate of Return
Years = length(mC)/12;
npv1 = calc_npv(mCbaseline - mC, 2, (1/12:1/12:Years)');
npv2 = calc_npv(mCbaseline - mC, -1+1e-8, (1/12:1/12:Years)');
if npv1>0
    %likely no solution for irr
    irr = -1;
elseif npv2<0
    %negative rate of return
    irr = -1;
else
    lowerBound = -1+1e-8;
    upperBound = 1;
    irr = 0;
    totalSaved = sum(mCbaseline) - sum(mC);
    npv = calc_npv(mCbaseline - mC, 0, (1/12:1/12:Years)');
    while abs(npv)/abs(totalSaved) > 1e-3
        if npv > 0
            lowerBound = irr;
            irr = irr + min(.1,(upperBound-lowerBound)/2);
        else
            upperBound = irr;
            irr = irr - min(.1, (upperBound-lowerBound)/2);
        end
        npv = calc_npv(mCbaseline - mC, irr, (1/12:1/12:Years)');
    end
end
% Payback period
if irr<=0
    payback = inf;
else
    payback = round(12*Years/2)/12;%round to nearest month
    best = false;
    npv = calc_npv(mCbaseline(1:payback*12) - mC(1:payback*12), rate, (1/12:1/12:payback)');
    while ~best
        if npv<0
            payback = payback+1/12;
        else
            payback = payback - 1/12;
        end
        npv_new = calc_npv(mCbaseline(1:payback*12) - mC(1:payback*12), rate, (1/12:1/12:payback)');
        if payback == Years
            best = true;
            payback = inf;
        end
        if abs(npv_new)>npv
            best = true;
        end
    end
end
end%Ends function Financial Metric

%% External functions
% Calculate net present value
function npv = calc_npv(npb, r, t)
% npb: Net Present Benefits = Benefits(t) - Costs(t)
% r: rate of return
% t: timestep (usually year)
npv = sum(npb ./ (1+r) .^ t);
end
