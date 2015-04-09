function TimeStampMatrix=create_TimeStampMatrix_from_textdata(textdata)



% for some reason it seems that this only works of the options containing fewer charaters are tried first. 
% if the input had 10:00:00 only, (without a date) and I put sscanf(textdata{i,1},'%f:%f:%f') at the end, the second
% try statement would be executed and TimeStampMatrix would have only one row.  

for i=1:length(textdata)
  try
    TimeStampMatrix(:,i) = sscanf(textdata{i,1},'%f:%f:%f');
  catch exception1
    try
      TimeStampMatrix(:,i) = sscanf(textdata{i,1},'"%f:%f:%f"');
    catch exception2
    try
      TimeStampMatrix(:,i) = sscanf(textdata{i,1},'"%f/%f/%f,%f:%f:%f"');
    catch exception3
      try 
        TimeStampMatrix(:,i) = sscanf(textdata{i,1},'%f/%f/%f,%f:%f:%f');
      catch exception4 
        try   
          TimeStampMatrix(:,i) = sscanf(textdata{i,1},'%f/%f/%f %f:%f:%f');  
        catch exception5
          try
            TimeStampMatrix(:,i) = sscanf(textdata{i,1},'"""%f/%f/%f,%f:%f:%f"""');
          catch exception6
            
            end
          end
        end
       end  
     end
  end
end