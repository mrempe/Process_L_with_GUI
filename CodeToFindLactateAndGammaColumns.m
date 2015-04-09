tLines1=fgetl(fid);  % get the first line
tLines2=fgetl(fid);  % get the second line
if ~isempty(strfind(tLines2,'[nA]')) %if textdata{2,3}=='[nA]' %nA is nanoamps, the units for lactate.
    LactateYesOrNo=2; %this tells the DetectArtifact function that there is a Lactate column in the data
    RawLactate=subData(1:end,1);
    [SmoothedLactatePercent(FileCounter),SmoothedLactate]=SmootheLactate(RawLactate);
    Lactate=SmoothedLactate(:,2);
end
if  ~isempty(strfind(tLines1,'Gamma'))
    ColumnsHeads = textscan(tLines1,'%s','delimiter', sprintf('\t'));
    HeadChars=char(ColumnsHeads{1,1});
    for i=1:length(HeadChars)
        Gamma(i)=~isempty(strfind(HeadChars(i,:),'Gamma'));
        EEG1(i)=~isempty(strfind(HeadChars(i,:),'EEG 1'));
        EEG2(i)=~isempty(strfind(HeadChars(i,:),'EEG 2'));
    end
    
    GammaEEG1Col=intersect(find(Gamma>0),find(EEG1>0))-2; %Find Gamma columns for EEG1.  Subtract two to account for the fact that the first two columns are time stamp and state.
    GammaEEG2Col=intersect(find(Gamma>0),find(EEG2>0))-2; %Find Gamma columns for EEG1.  Subtract two to account for the fact that the first two columns are time stamp and state.
    
    GammaEEG1=subData(IntervalStart:IntervalStop,GammaEEG1Col(1:end));
    GammaEEG2=subData(IntervalStart:IntervalStop,GammaEEG2Col(1:end));
end %if clause that will find and extract Gamma values from the file.