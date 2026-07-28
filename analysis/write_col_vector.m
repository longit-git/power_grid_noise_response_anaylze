function write_col_vector(fn, vec)
% write_col_vector 将向量写为 tab 分隔文本（单行）
%
% 用法：
%   write_col_vector('b.txt', b_vec)

fid = fopen(fn, 'w');
if fid == -1
    error('write_col_vector:CannotOpen', '无法打开文件: %s', fn);
end
fprintf(fid, '%.12e\t', vec(1:end-1));
fprintf(fid, '%.12e\n', vec(end));
fclose(fid);
end
