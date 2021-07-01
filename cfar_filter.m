function [sigbin] = cfar_filter(signal, K0, gwndw, rwndw)
    % K0 -> CFAR threshold constant.
    % gwndw -> size of guard window. 
    % rwndw -> size of reference window.
    sigbin = zeros(length(signal), 1);
    if mod(rwndw, 2) ~= 0
       rwndw = rwndw + 1;
    end
    if mod(gwndw, 2) ~= 0
       gwndw = gwndw + 1;
    end
    for i=1:length(signal)
        intv0 = max([1, i - (gwndw/2 + rwndw/2)]);
        intv1 = min( [ length(signal), i + (gwndw/2 + rwndw/2)]);
        
        intv_grd_0 = max([1, i - gwndw/2]);
        intv_grd_1 = min( [length(signal), i + gwndw/2] );
        nel_avg =  (intv1-intv0) - (intv_grd_1-intv_grd_0) ;
        averager = ( sum( abs( signal(intv0:intv1) ) ) - sum( abs( signal( intv_grd_0:intv_grd_1 ) ) ) )/nel_avg;
        sigbin(i) = signal(i) > K0*averager;
    end

end

