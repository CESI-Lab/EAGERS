function out = eval_curve(curve,name,input)
k = f_index(name,curve.name);
x = min(curve.max_x(k),max(curve.min_x(k),input(:,1)));
switch curve.type{k}
    case 'Quadratic'
        out = curve.a0(k) + curve.a1(k)*x + curve.a2(k)*x.^2;
    case 'Cubic'
        out = curve.a0(k) + curve.a1(k)*x + curve.a2(k)*x.^2 + curve.a3(k)*x.^3;
    case 'Biquadratic'
        y = min(curve.max_y(k),max(curve.min_y(k),input(:,2)));
        out = curve.a0(k) + curve.a1(k)*x + curve.a2(k)*x.^2 + curve.b1(k)*y + curve.b2(k)*y.^2 + curve.ab(k)*x.*y;
    case 'Bicubic'
        y = min(curve.max_y(k),max(curve.min_y(k),input(:,2)));
        out = curve.a0(k) + curve.a1(k)*x + curve.a2(k)*x.^2 + curve.a3(k)*x.^3 + curve.b1(k)*y + curve.b2(k)*y.^2 + curve.b3(k)*y.^3 + curve.ab(k)*x.*y + curve.aab(k)*x.^2.*y + curve.abb(k)*x.*y.^2;
    case 'ChillerPartLoadWithLift'
        y = min(curve.max_y(k),max(curve.min_y(k),input(:,2)));
        z = input(:,3);
        out = curve.a0(k) + curve.a1(k)*x + curve.a2(k)*x.^2 + curve.a3(k)*x.^3 + curve.b1(k)*y + curve.b2(k)*y.^2 + curve.b3(k)*y.^3 + curve.ab(k)*x.*y + curve.aab(k)*x.^2.*y + curve.abb(k)*x.*y.^2 + curve.aabb(k)*x.^2.*y.^2 + curve.cb3(k)*z.*y.^3;
end
end%Ends function eval_curve