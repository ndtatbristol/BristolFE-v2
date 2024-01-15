function fn_simple_text_progress_bar(n, N)
% persistent bytes_written
% if n <= 1
%     bytes_written = 0;
% end
% 
max_char_cnt = 10;
done_char = '.';
not_done_char = '-';
old_done = round((n-1) / N * max_char_cnt);
now_done = round( n    / N * max_char_cnt);
if now_done == old_done
    return
else
    fprintf(done_char);
end
% end
% for i = 1:bytes_written
%     fprintf('\b');
% end
% bytes_written = 0;
% 
% for i = 1:now_done
%     bytes_written = bytes_written + fprintf(done_char);
% end
% for i = 1:(char_cnt - now_done)
%     bytes_written = bytes_written + fprintf(not_done_char);
% end
end