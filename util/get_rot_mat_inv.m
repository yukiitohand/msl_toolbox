function [rot_mat] = get_rot_mat_inv(rolld,pitchd,yawd)
% create a rotation matrix for a given roll, pitch, and yaw angles
%  INPUTS:
%   rolld, pitchd, yawd: angles in degree
%  OUTPUTS:
%   rot_mat: rotation matrix

yawinv_mat = [...
    cosd(-yawd) -sind(-yawd) 0;
    sind(-yawd)  cosd(-yawd) 0;
          0              0       1
];

pitchinv_mat = [...
     cosd(-pitchd)  0  sind(-pitchd);
            0          1        0       ;
    -sind(-pitchd)  0  cosd(-pitchd);
];

rollinv_mat = [...
    1        0               0       ;
    0  cosd(-rolld) -sind(-rolld);
    0  sind(-rolld)  cosd(-rolld)
];

rot_mat = rollinv_mat * pitchinv_mat * yawinv_mat;

end