function write_row_vector(fn, vec)
% write_row_vector 将向量写为 tab 分隔文本（单行）
%
% 用法：
%   write_row_vector('M_list.txt', M_list)

fid = fopen(fn, 'w');
if fid == -1
    error('write_row_vector:CannotOpen', '无法打开文件: %s', fn);
end
fprintf(fid, '%.12e\t', vec(1:end-1));
fprintf(fid, '%.12e\n', vec(end));
fclose(fid);
end
