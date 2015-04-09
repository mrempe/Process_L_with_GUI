function display_progress(step,total_steps)

if(step==ceil(total_steps/4))
   disp('1/4 of the way through')
end

if(step==ceil(total_steps/2))
   disp('1/2 of the way through')
end

if(step==ceil(total_steps*.75))
   disp('3/4 of the way through')
end
