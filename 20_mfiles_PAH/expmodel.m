function [z] = expmodel(ar1,Cr1,x)
% this is the function for the joint double powered exponential covariance model

    z = Cr1.*( exp(-(3*x./ar1)) );

end