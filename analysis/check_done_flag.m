function tf = check_done_flag(flag_file, tag, i_setting)
% check_done_flag 检查 stage 完成标记；存在则打印跳过信息并返回 true
%
% 用法：
%   if check_done_flag(flag_file, 'HD', i_setting), return, end

tf = exist(flag_file, 'file');
if tf
    fprintf('  [Setting %d] %s 已存在，跳过\n', i_setting, tag);
end
end
