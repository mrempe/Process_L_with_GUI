% script to test datetime stuff

c = find(textdata{1,1}==':');   % Find all locations of the colon in the first time stamp
s = find(textdata{1,1}=='/');
first_colon_loc  = c(1);   
second_colon_loc = c(2);
first_slash_loc  = s(1);
second_slash_loc = s(2);
hour_first_time_stamp    = str2num(textdata{1,1}(first_colon_loc-2:first_colon_loc-1));  
minute_first_time_stamp  = str2num(textdata{1,1}(first_colon_loc+1:first_colon_loc+2));
second_first_time_stamp  = str2num(textdata{1,1}(last_colon_loc+1:last_colon_loc+2));


% year
Y = 2000 + str2num(textdata{1,1}(second_slash_loc+3:second_slash_loc+4));
% month
M = str2num(textdata{1,1}(1:2));

for i=1:length(textdata)
	D(i) = str2num(textdata{i,1}(first_slash_loc+1:first_slash_loc+2));
	H(i) = str2num(textdata{i,1}(first_colon_loc-2:first_colon_loc-1)); %or try str2num(textdata{i,1}(first_colon_loc-2:first_colon_loc-1))
	m(i) = str2num(textdata{i,1}(second_colon_loc-2:second_colon_loc-1));
	s(i) = str2num(textdata{i,1}(second_colon_loc+1:second_colon_loc+2));
end



t=datetime(Y,M,D,H,m,s);  % this is a datetime vector with 