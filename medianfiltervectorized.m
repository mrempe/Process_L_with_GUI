function mY=medianfiltervectorized(y,MedianWidth)
% function mY=mymedianfilter(y,MedianWidth)
% Performs a median-based filter operation that replaces 
% each element of y with the median of the MedianWidth 
% adjacent points. (MedianWidth points to the left and MedianWidth
% points to the right).  MedianWidth must be a positive integer. 

% I've made improvements on a code from Tom O'Haver.
% 1) The code is cleaner and easier to read. (and faster)
% 2) I fixed an error in the original version that made the 
%    computation biased to the right.
% 2) I hope to vectorize the code too. (done in this version!)



if(MedianWidth<=0)
  error('MedianWidth must be a positive integer')
end

MW=round(MedianWidth);
total_width = 2*MedianWidth+1;  %MedianWidth is the number of 
                                % points on either side. total_width 
                                % is the total number of points in sliding
                                % window. 
ly=length(y);
S=zeros(ly-total_width+1,total_width);
for i=1:total_width
S(:,i) = y(i:ly-total_width+i);
end
mY=median(S,2);
mY=[y(1:MW); mY; y(end-MW+1:end)];



