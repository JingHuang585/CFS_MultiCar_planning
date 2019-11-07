% This function uses two reference points on the reference trajectory to generate
% reference points for convexifing feasible set.
function x_ref = get_ref(x1_ref, x2_ref, margin)
    dir = (x2_ref - x1_ref) / norm(x2_ref - x1_ref);
    x1_ref_new = x1_ref;
    x2_ref_new = x1_ref_new + dir * margin;
    x_ref = [x1_ref_new; x2_ref_new];
end