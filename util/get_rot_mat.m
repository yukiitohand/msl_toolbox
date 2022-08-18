function [rot_mat] = get_rot_mat(rolld,pitchd,yawd)
% create a rotation matrix for a given roll, pitch, and yaw angles
%  INPUTS:
%   rolld, pitchd, yawd: angles in degree
%  OUTPUTS:
%   rot_mat: rotation matrix

yaw_mat = [...
    cosd(yawd) -sind(yawd) 0;
    sind(yawd)  cosd(yawd) 0;
          0              0       1
];

pitch_mat = [...
     cosd(pitchd)  0  sind(pitchd);
            0          1        0       ;
    -sind(pitchd)  0  cosd(pitchd);
];

roll_mat = [...
    1        0               0       ;
    0  cosd(rolld) -sind(rolld);
    0  sind(rolld)  cosd(rolld)
];

rot_mat = yaw_mat * pitch_mat * roll_mat;
% rot_mat = roll_mat * pitch_mat * yaw_mat;

end