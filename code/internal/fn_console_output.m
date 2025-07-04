function fn_console_output(short_comment, varargin)
if numel(varargin) < 1 || isempty(varargin{1})
    long_comment = short_comment;
else
    long_comment = varargin{1};
end
if numel(varargin) < 2
    indent = 1;
else
    indent = varargin{2};
end
comment_indent_level = fn_get_comment_indent_level;
if indent
    indent_str = repmat('  - ', 1, comment_indent_level);
else
    indent_str = '';
end
% fprintf(['Indent level %i, ', short_comment],  comment_indent_level);
switch fn_get_comment_verbosity
    case 'low'
        fprintf([indent_str, short_comment]);
    case 'high'
        fprintf([indent_str, long_comment]);
    otherwise
end
end