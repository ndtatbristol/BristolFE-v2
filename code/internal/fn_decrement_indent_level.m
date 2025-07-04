function fn_decrement_indent_level
global COMMENT_INDENT_LEVEL
if exist('COMMENT_INDENT_LEVEL', "var") && ~isempty(COMMENT_INDENT_LEVEL)
    COMMENT_INDENT_LEVEL = COMMENT_INDENT_LEVEL - 1;
else
    COMMENT_INDENT_LEVEL = 0;
end
end