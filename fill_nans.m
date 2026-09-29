function [U] = fill_nans(UwithNan)
%   Remove nan's by average of neighbors (3D). When there are several nans
%   touching, we run several iterations
U = UwithNan;
while any(isnan(U(:)))
    out = isnan(U);
    ind = find(out);
    [I,J,K] = ind2sub(size(U),ind);
    Im1 = I-1;
    Ip1 = I+1;
    Jm1 = J-1;
    Jp1 = J+1;
    Km1 = K-1;
    Kp1 = K+1;

    [top_ind,bot_ind,left_ind,right_ind,fwd_ind,bck_ind] = deal(nan(size(ind)));
    [top_val,bot_val,left_val,right_val,fwd_val,bck_val] = deal(nan(size(ind)));

    oob_top =  (Im1 < 1);
    top_ind(~oob_top) = sub2ind(size(U),Im1(~oob_top),J(~oob_top),K(~oob_top));
    top_val(~oob_top) = U(top_ind(~oob_top));

    oob_bot =  (Ip1 > size(U,1));
    bot_ind(~oob_bot) = sub2ind(size(U),Ip1(~oob_bot),J(~oob_bot),K(~oob_bot));
    bot_val(~oob_bot) = U(bot_ind(~oob_bot));

    oob_left =  (Jm1 < 1);
    left_ind(~oob_left) = sub2ind(size(U),I(~oob_left),Jm1(~oob_left),K(~oob_left));
    left_val(~oob_left) = U(left_ind(~oob_left));

    oob_right =  (Jp1 > size(U,2));
    right_ind(~oob_right) = sub2ind(size(U),I(~oob_right),Jp1(~oob_right),K(~oob_right));
    right_val(~oob_right) = U(right_ind(~oob_right));

    oob_bck =  (Km1 < 1);
    bck_ind(~oob_bck) = sub2ind(size(U),I(~oob_bck),J(~oob_bck),Km1(~oob_bck));
    bck_val(~oob_bck) = U(bck_ind(~oob_bck));

    oob_fwd =  (Kp1 > size(U,3));
    fwd_ind(~oob_fwd) = sub2ind(size(U),I(~oob_fwd),J(~oob_fwd),Kp1(~oob_fwd));
    fwd_val(~oob_fwd) = U(fwd_ind(~oob_fwd));

    mval = nanmean([top_val,bot_val,left_val,right_val,bck_val,fwd_val],2);
    U(out) = mval;
%     display('Iter!')
end
end