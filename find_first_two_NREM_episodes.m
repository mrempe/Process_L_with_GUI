function index2NREMend = find_first_two_NREM_episodes(state_data)

% This function just finds the first two NREM (SWS) episodes 
% longer than 1 minute each and returns the index of the 
% end of the second NREM sleep episode that is 1 minute 
% or longer.  

runs = contiguous(state_data,1); %count all contiguous runs, even of only 1
if isempty(runs)   % if no runs of NREM were found, set index2NREMend to 1 so no data is cut off
	index2NREMend=1;
else

	allruns = runs{1,2};  % go from cell array to matrix

	counter = 0;  % counter for the number of NREM episodes 1 min or longer
	i=0;

	while counter < 2 & i<length(state_data)-1
		i=i+1;
		if(allruns(i,2)-allruns(i,1))>=5
			counter = counter+1;
		end
	end

	index2NREMend = allruns(i,2);

end
