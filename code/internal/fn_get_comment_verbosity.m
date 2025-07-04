function comment_verbosity = fn_get_comment_verbosity
global COMMENT_VERBOSITY
if exist(COMMENT_VERBOSITY,"var")
    comment_verbosity = COMMENT_VERBOSITY;
else
    comment_verbosity = 'low';
end
end