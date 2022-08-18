function [rolld,pitchd,yawd] = get_RollPitchYawd_from_rotation_matrix(rot_mat)
% [rolld,pitchd,yawd] = get_RollPitchYaw_from_rotation_matrix(rot_mat)
%  Get 
%  INPUT
%    rot_mat: 3x3 rotation matrix
%  OUTPUT
%    rolld, pitchd, yawd: roll, pitch, yaw angles in degree.

yawd   = atan2d(rot_mat(2,1),rot_mat(1,1));
pitchd = atan2d(-rot_mat(3,1),sqrt(rot_mat(3,2)^2+rot_mat(3,3)^2));
rolld  = atan2d(rot_mat(3,2),rot_mat(3,3));

end