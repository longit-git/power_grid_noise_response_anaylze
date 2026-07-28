function write_done_flag(flag_file)
% write_done_flag 写 stage 完成标记（空文件）

fclose(fopen(flag_file, 'w'));
end
