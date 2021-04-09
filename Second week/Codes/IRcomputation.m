function [maxIR, index]  = IRcomputation(NPort,weights,returns_p,returns_bond,returns_b)
%IR computation
%NPort : 25
%weights : portfolio 25x9
%returns : historical data 352x9

rend_p = 0.65*weights*(returns_p') + repmat(0.35*returns_bond,NPort,1) ; 
%25 x 9   9 x 350 --> 25 x 350 1st riga = 1 portafoglio, 1st colonna ultima data
rend_b = [0.60 0.40]*(returns_b');
rend_b = repmat(rend_b,25,1);
diff = (rend_p - rend_b)'; %colonne = portafogli
IR = zeros(1,NPort);
for ii = 1:NPort
    IR(ii) = mean(diff(:,ii)) / std(diff(:,ii));
end 
maxIR = max(IR);
index = find(IR == max(IR));
end 