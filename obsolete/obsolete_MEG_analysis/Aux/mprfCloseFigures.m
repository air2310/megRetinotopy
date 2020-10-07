function mprfCloseFigures(keep)

all_handles = get(0,'Children');

cl_figs = setdiff(all_handles, keep);
close(cl_figs)

end





















