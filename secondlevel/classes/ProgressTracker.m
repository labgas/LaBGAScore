classdef ProgressTracker < handle
   properties
       count = 0
       total = 0
       lastPrint
       startTime
   end
   methods
       function obj = ProgressTracker(total)
           obj.total     = total;
           obj.lastPrint = tic;
           obj.startTime = tic;
       end
       
       function update(obj)
           % Increment count
           obj.count = obj.count + 1;

           pct = obj.count / obj.total;

           % Only print when enough time passed or at the end
           if toc(obj.lastPrint) > 0.1 || obj.count == obj.total
               
               % Compute elapsed time
               elapsed = toc(obj.startTime);

               % Compute ETA
               if pct > 0
                   remaining = elapsed * (1/pct - 1);
               else
                   remaining = NaN;
               end

               eta_str = obj.formatTime(remaining);

               % Print progress line
               fprintf(1, '\rProgress: %d/%d (%.1f%%) | ETA: %s', ...
                   obj.count, obj.total, pct*100, eta_str);

               drawnow limitrate;  % IMPORTANT for printing during parfor
               obj.lastPrint = tic;

               if obj.count == obj.total
                   fprintf('\n'); % final newline
               end
           end
       end

       function str = formatTime(~, seconds)
           if isnan(seconds)
               str = '--:--:--';
               return;
           end

           hrs = floor(seconds/3600);
           mins = floor(mod(seconds,3600)/60);
           secs = floor(mod(seconds,60));

           str = sprintf('%02d:%02d:%02d', hrs, mins, secs);
       end
   end
end


