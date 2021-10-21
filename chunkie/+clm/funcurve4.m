function [r,d,d2] = funcurve4(t,icurve,cpars)
%%funcurve
% return position, first and second derivatives of four curve segments
% in the embedded eye layered media problem.
%
% Inputs:
% t - paramter values to evaluate these quantities
%
% Outputs:
% r - coordinates
% d - first derivatives w.r.t. t
% d2 - second derivatives w.r.t. t


if icurve == 1 || icurve == 2
  [r,d,d2] = clm.complexx6(t,icurve,cpars);
  %[r,d,d2] = clm.complexx1(t,icurve,cpars);
elseif icurve == 3 || icurve == 4
  [r,d,d2] = clm.circulararc(t,cpars);
  %[r,d,d2] = clm.eyeshape(t,icurve,cpars);
end

end