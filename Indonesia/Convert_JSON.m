% s = struct('xsto1', xsto1)
% jsonencode(s)
% 
% txt = jsonencode(s)
% 
% 
% 
% prm.xsto3=xsto3
% prm.data=data
% 
% s = struct('prm', prm)
% txt = jsonencode(s);


% From Carel on 1 April 2024
%load('JSON_data.mat')
s = struct('data',data, 'r', r, 'p',p, 'prm',prm, 'xs',xs)
jsonencode(s)

tmp = jsonencode(s)
fid = fopen('JSON_data.json', 'w');
fprintf(fid, '%s', tmp);
fclose(fid);





% 
% %%%%%%%%%%% Below is the code to transfer xsto1.mat to sandip.json
% load('xsto1.mat')
% s = struct('xsto1', xsto1)
% jsonencode(s)
% 
% tmp = jsonencode(s)
% fid = fopen('sandip.json', 'w');
% fprintf(fid, '%s', tmp);
% fclose(fid);
