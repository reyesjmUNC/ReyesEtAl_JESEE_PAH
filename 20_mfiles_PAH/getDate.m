function [ yr mo da ] = getDate(dates)
% takes dates in YYYYMMDD format and converts to YYYY MM DD
yr = floor(dates/10000);
mo = floor( (dates - yr*10000)/100 );
da = dates - yr*10000 - mo*100;

end