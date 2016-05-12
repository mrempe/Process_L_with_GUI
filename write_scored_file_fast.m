function write_scored_file_fast(filename,directory,predicted_score)
%







% First copy the original file so we don't mess it up
 a = find(filename=='.');
 d = find(filename=='\');
 newfilename = strcat(directory,'\',filename(d(end)+1:a(end)-1), 'AUTOSCORED', filename(a(end):end));