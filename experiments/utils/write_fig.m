function write_fig(num,name)
sfigure(num);
fname = [name '.eps'];
print('-depsc',fname)
[s, w] = unix(['LD_LIBRARY_PATH="" epstopdf "' fname '"']);
assert(~s, w);
fprintf('Wrote %s.{eps,pdf}\n',name)

pdfname = [name '.pdf'];
[s, w] = unix(['LD_LIBRARY_PATH="" pdfcrop ' pdfname ' ' pdfname]);
assert(~s, w);
fprintf('Cropped %s\n', pdfname)
