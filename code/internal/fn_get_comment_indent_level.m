function comment_indent_level = fn_get_comment_indent_level
global COMMENT_INDENT_LEVEL
if exist('COMMENT_INDENT_LEVEL', "var") && ~isempty(COMMENT_INDENT_LEVEL)
    comment_indent_level = COMMENT_INDENT_LEVEL;
else
    comment_indent_level = 0;
end
end