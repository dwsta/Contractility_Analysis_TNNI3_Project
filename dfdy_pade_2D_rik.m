function dfdy = dfdy_pade_2D_rik(f,dy)
% See dfdx_pade_2D_rik
dfdy_t = dfdx_pade_2D_rik(f',dy);
dfdy = dfdy_t';
