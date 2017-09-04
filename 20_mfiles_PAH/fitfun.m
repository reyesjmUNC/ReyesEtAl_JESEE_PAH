function error = fitfun(CPM,CPAH,CX,a1,CPM1,CPAH1,CX1,Lag)
    y1 = CPM1.*( exp(-(3*Lag./a1)) );
    y2 = CPAH1.*( exp(-(3*Lag./a1)) );
    y3 = CX1.*( exp(-(3*Lag./a1)) );
    error = sum((CPM-y1).^2 + (CPAH-y2).^2 + (CX-y3).^2);
end