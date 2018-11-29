%%% USES the 36x36 solid fields to create a 72x72 matrix
% essentiall just multiplying the cells

% CaCO3
worbe2_Fsed_caco3_36x36 = load('worbe2_Fsed_caco3.36x36','ascii');

[m,n]=size(worbe2_Fsed_caco3_36x36);

for i = 1:m,
    for j=1:n,
        worbe2_Fsed_caco3_72x72(2*i-1,2*j-1) = worbe2_Fsed_caco3_36x36(i,j);
        worbe2_Fsed_caco3_72x72(2*i-1,2*j) = worbe2_Fsed_caco3_36x36(i,j);
        worbe2_Fsed_caco3_72x72(2*i,2*j-1) = worbe2_Fsed_caco3_36x36(i,j);
        worbe2_Fsed_caco3_72x72(2*i,2*j) = worbe2_Fsed_caco3_36x36(i,j);
    end
end
dlmwrite('worbe2_Fsed_caco3.72x72', worbe2_Fsed_caco3_72x72,'delimiter','\t')

%detrital
worbe2_Fsed_det_36x36 = load('worbe2_Fsed_det.36x36','ascii');

[m,n]=size(worbe2_Fsed_det_36x36);

for i = 1:m,
    for j=1:n,
        worbe2_Fsed_det_72x72(2*i-1,2*j-1) = worbe2_Fsed_det_36x36(i,j);
        worbe2_Fsed_det_72x72(2*i-1,2*j) = worbe2_Fsed_det_36x36(i,j);
        worbe2_Fsed_det_72x72(2*i,2*j-1) = worbe2_Fsed_det_36x36(i,j);
        worbe2_Fsed_det_72x72(2*i,2*j) = worbe2_Fsed_det_36x36(i,j);
    end
end
dlmwrite('worbe2_Fsed_det.72x72', worbe2_Fsed_det_72x72,'delimiter','\t')

%opal
worbe2_Fsed_opal_36x36 = load('worbe2_Fsed_opal.36x36','ascii');

[m,n]=size(worbe2_Fsed_opal_36x36);

for i = 1:m,
    for j=1:n,
        worbe2_Fsed_opal_72x72(2*i-1,2*j-1) = worbe2_Fsed_opal_36x36(i,j);
        worbe2_Fsed_opal_72x72(2*i-1,2*j) = worbe2_Fsed_opal_36x36(i,j);
        worbe2_Fsed_opal_72x72(2*i,2*j-1) = worbe2_Fsed_opal_36x36(i,j);
        worbe2_Fsed_opal_72x72(2*i,2*j) = worbe2_Fsed_opal_36x36(i,j);
    end
end

dlmwrite('worbe2_Fsed_opal.72x72', worbe2_Fsed_opal_72x72,'delimiter','\t')
