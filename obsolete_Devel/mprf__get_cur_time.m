function cur_time = mprf__get_cur_time
cur_time = datestr(now);
cur_time(cur_time == ' ' | cur_time == ':' | cur_time == '-') = '_';



end