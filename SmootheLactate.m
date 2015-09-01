function[numchanged,SmoothedData]=SmootheLactate(InputVector);
%this function was made to work with a 1-dimensional matrix (a vector called InputVector) that contains
%only the lactate data of interest. 

SmoothedData(:,1)=InputVector;                  %Establish column 1 of matrix; this is the raw data.
SmoothedData(1:10,2)=InputVector(1:10);         %For first ten values of smoothed data, must output raw data because they lack ten previous values to compare to.

MatrixforMeanSEM=zeros(size(SmoothedData,1)-10,11); %create a matrix with ten fewer rows than original vector has rows.  We will subject to smoothing every value after the first ten.

for matrixcolumns=1:11  %create each column from the original vector shifted one point to the left per row. This aligns consecutive values in columns of ten.   
    MatrixforMeanSEM(:,matrixcolumns)=SmoothedData(matrixcolumns:length(SmoothedData)-(11-matrixcolumns),1); %  matrix row 1 now has values 1-11 in a row, row 2 has values 2-12 in a row, etc
end

FlipMatrix=MatrixforMeanSEM';  % Each column now represents one input row (ten consecutive values for calculating mean + SEM) from row 11 through to last row.

numchanged=abs(FlipMatrix(11,:)-mean(FlipMatrix(1:10,:)))./((std(FlipMatrix(1:10,:))/sqrt(1))*10)>1; %Does absolute value of ((value minus mean of previous five)/five*SEM of previous ten) exceed 1?
FlipMatrix(11,find(numchanged>0))=mean(FlipMatrix(1:10,find(numchanged>0))); %if yes, reset means in flipped Matrix
SmoothedData(11:end,2)=FlipMatrix(11,:)';  %now replace values in row 2 with flipped matrix values.
numchanged=sum(numchanged); %sums the values in numchanged to tally total number of changes.
return
