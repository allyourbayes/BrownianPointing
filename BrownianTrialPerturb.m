function [] = BrownianTrialPerturb(initials,sinPercent,tfreq,stepSize,updateRate,numConv)


%BrownianPointing - 
% ___________________________________________________________________
% 
% Try to follow a dot with the mouse.
% _______________________________________________________________________
%

% HISTORY
% 4/23/13   lkc Wrote it based on MouseTraceDemo 
% 4/25/13   lkc made ver 3, which does real rotation, and needs a new
%               naming convention
% 4/26/13   lkc Going to implement some radness
%               - sinusoidal varying of parameters
%               - treating each frame as a trial + using lag in a revco way
%               the goal being to get a quantitative measure of learning
%               and/or sensorimotor adaptation.
% 4/28/13   lkc Implemented sinusoidal variation
% 4/29/13   lkc Origin of transform can now mooove wif cursor, such that
%               each frame's movement is like a little mini trial starting
%               at the origin
%   7/29    Tested handling actual cursor goes off-screen, place back in
%           center and handle drawn cursor smoothly

% NOTES
%   4/28    - need to implement sinusoidal variation
%           & *really* need to code a file saving & analysis pipeline
%   4/29    - I think a cool condition would be to have the target jump a
%           relatively long way, and sit there until until the subject's
%           responds "stops" by come criteria.
%   4/29    - Yummy!.
%   7/29    - Cursor path flipped for y but not x when looking at raw
%           traces because x,y = (0,0) starts in top left by default. Logic
%           of analyses ought to be the same, though.



%% set up initial parameters
trialTime = 32*8; %20;      % trial duration in s
%stepSize = 8;      % ave or max stepsize in pixels - determines target speed
frameRate = 60;     % will determine this "live" in future versions
totFrames = trialTime*frameRate;    % yup
%tfreq = 0.25; %0.2;  % temp. freq. of parameter variation
%updateRate = 3;

% these "maximums" are the extreme vals of the parameters in the dynamic
% case, and are simply the parameters in the static case
aMax = 0.6;
bMax = 0.6;
%thMax = 5*pi/12; 
thMax = pi/4;

% flags for stuff that will become arguments or gui checkboxes at some
% point...
mouseVisibleFlag = 1;   % visual feedback of whey yo finger at
dynoUpdateFlag = 1;     % dynamic updating of transform
scalerFlag = 0;        % scaler or rotation transformation?
origAtCurs = 1;         % do transform around cursor or screen center?

% dot colorz & sizes
dotSizes = [30, 10]; % first one is target
%dotColors = [rTarg, rMaus; gTarg, gMaus; bTarg, bMaus];
if mouseVisibleFlag
    dotColors = [20, 200; 200, 200; 20, 200];
else
    dotColors = [20, 0; 200, 0; 20, 0];
end

% D-fine the tform
if scalerFlag   % Kat's scaler rotation
    a = 1;
    b = aMax;
    c = bMax;
    d = 1;
    if dynoUpdateFlag
        b = aMax.*sin(2.*pi.*tfreq.*(1/frameRate).*(1:totFrames));
        c = bMax.*(pi./6).*sin(2.*pi.*tfreq.*(1/frameRate).*(1:totFrames));
    else
        b = repmat(b, 1, totFrames);
        c = repmat(c, 1, totFrames);
%         m = [a, b; c, d];
    end
    
else            % a rotation transformation
    th = thMax; %pi./2;
    if dynoUpdateFlag
%         th = linspace(0, pi./2, totFrames);
        th = thMax.*sin(2.*pi.*tfreq.*(1/frameRate).*(1:totFrames));
    else
        th = repmat(th, 1, totFrames);
    end
    
    m = [cos(th(1)), -sin(th(1)); sin(th(1)), cos(th(1))]; % le rotation
    
end % end rotation type selection

% some analysis params
maxLag =  .7;    % specify max xcorr lag in secs
frontPorch = 1;     % # of initial secs to drop from the xcorr

%% generate target sequence
% generate position jitter directly
% targCoords = [0;0];
% for i = 1:totFrames-1
%     % need to work in step size here
%     temp = stepSize*randn(2,1);        
%     targCoords(:,i+1) = targCoords(:,i) + temp;
% end

% or, to drop out some of the high freq. jitter, let's generate random
% velocities, smooove them, and then integrate them to get positions.

%{

for i = 1:totFrames
    temp(i,:) = randn(2,1);
    latestCoords = stepSize*sum(temp, 2);
    if temp(1,end) > 0.95*theRect(RectRight)
       temp(1,end) = randn(1,1)-0.5;
    elseif temp(1,end) < 0.05*theRect(RectRight)
       temp(1,end) = randn(1,1)+0.5;
    end
        if temp(2,end) > 0.95*theRect(RectBottom)
       temp(2,end) = randn(1,1)-0.5;
    elseif temp(2,end) < 0.05*theRect(RectBottom)
       temp(2,end) = randn(1,1)+0.5;
    end
end

%temp = randn(2,totFrames);
temp(1,:) = conv(temp(1,:), [.25,.25,.25,.25], 'same');
temp(2,:) = conv(temp(2,:), [.25,.25,.25,.25], 'same');
temp(:,1) = [0; 0];
targCoords = stepSize*cumsum(temp, 2);


respCoords = zeros(size(targCoords));   % response coordinate array
tFormRespCoords = respCoords;

%}

%% stimulus presentation
try
    
    % Open up a window on the screen and clear it.
    whichScreen = max(Screen('Screens'));
    [theWindow,theRect] = Screen(whichScreen,'OpenWindow',0,[],[],2);

    % Move the cursor to the center of the screen
    centerX = theRect(RectRight)/2;
    centerY = theRect(RectBottom)/2;
    
    
    for i = 1:totFrames
    temp(:,i) = randn(2,1);
    latestCoords = stepSize*sum(temp, 2);
    if latestCoords(1,1) + centerX > 0.95*theRect(RectRight)
       temp(1,end) = randn(1,1)-0.5;
    elseif latestCoords(1,1) + centerX < 0.05*theRect(RectRight)
       temp(1,end) = randn(1,1)+0.5;
    end
        if latestCoords(2,1) + centerY > 0.95*theRect(RectBottom)
       temp(2,end) = randn(1,1)-0.5;
    elseif latestCoords(2,1) + centerY < 0.05*theRect(RectBottom)
       temp(2,end) = randn(1,1)+0.5;
    end
    end

%numConv = 15;
%temp = randn(2,totFrames);
temp(1,:) = conv(temp(1,:), repmat(1/numConv,1,numConv), 'same');
temp(2,:) = conv(temp(2,:), repmat(1/numConv,1,numConv), 'same');
temp(:,1) = [0; 0];
targCoords = stepSize*cumsum(temp, 2);


respCoords = zeros(size(targCoords));   % response coordinate array
tFormRespCoords = respCoords;
    
    
    targCoords(1,:) = targCoords(1,:) + centerX;
    targCoords(2,:) = targCoords(2,:) + centerY;
    SetMouse(centerX,centerX);
    HideCursor();

    % Wait for a click
    Screen(theWindow,'FillRect',0);
    Screen(theWindow,'DrawText','Klick to start trial',50,50,255);
    Screen('Flip', theWindow);
    while(1)
        [x,y,buttons] = GetMouse(whichScreen);
        if buttons(1)
            break
        end
    end
    
    % short delay before the stimulus begins (make random)
    pause(0.5);

    % Loop and track the mouse, drawing the contour
    SetMouse(centerX,centerY);
    x = centerX;
    y = centerY;
    % DRIFT CONCEPT ADDED 07/29
    drift_x = 0;
    drift_y = 0;
    tFormedCoords = [centerX; centerY];
    txOld = x;
    tyOld = y;
    
    % Query the frame duration
    ifi = Screen('GetFlipInterval', theWindow);

    % Flips on this number of frames (1 = flips every frame)
    waitframes = 1;
    
    Screen('DrawDots', theWindow, [centerX, centerX; centerY, centerY], dotSizes, dotColors);
    vbl = Screen('Flip', theWindow);
    %Screen('Flip', theWindow);
    sampleTime = 1/frameRate; %0.01;
    startTime = GetSecs;
    nextTime = startTime+sampleTime;
    
    pause(0.5); % ready... set...
    repeatCoords = [targCoords(1,1), tFormedCoords(1,1) ; targCoords(2,1), tFormedCoords(2,1)];
    
    for i = 1:totFrames
        xOld = x;
        yOld = y;
        [x,y,buttons] = GetMouse(whichScreen);
        % NEW PART ADDED TO DEAL WITH LEAVING SCREEN  07/29
        if x > 0.95*theRect(RectRight)
            drift_x = drift_x + (0.95-0.5)*theRect(RectRight);
            xOld = x - xOld + centerX;
            x = centerX;
            SetMouse(centerX,round(y),whichScreen);
        elseif x < 0.05*theRect(RectRight)
            drift_x = drift_x - (0.95-0.5)*theRect(RectRight);
            xOld = x - xOld + centerX;
            x = centerX;
            SetMouse(centerX,round(y),whichScreen);
        end
        if y > 0.95*theRect(RectBottom)
            drift_y = drift_y + (0.95-0.5)*theRect(RectBottom);
            yOld = y - yOld + centerY;
            y = centerY;
            SetMouse(round(x),centerY,whichScreen);
        elseif y < 0.05*theRect(RectBottom)
            drift_y = drift_y - (0.95-0.5)*theRect(RectBottom);
            yOld = y - yOld + centerY;
            y = centerY;
            SetMouse(round(x),centerY,whichScreen);
        end
        % END NEW PART  07/29
        
        if scalerFlag   % scaler tranformation
            m = [a, b(i); c(i), d];
        else    % rotation
            %sinPercent = 1;
            m = [cos(pi/4*round(th(i)/thMax)), -sin(pi/4*round(th(i)/thMax)); sin(pi/4*round(th(i)/thMax)), cos(pi/4*round(th(i)/thMax))]; % le rotation
            m = [cos(th(i)), -sin(th(i)); sin(th(i)), cos(th(i))]; % le rotation
            %m = [cos(thMax*(th(i)>0) - thMax*(th(i)<=0)), -sin(thMax*(th(i)>0) - thMax*(th(i)<=0));...
            %    sin(thMax*(th(i)>0) - thMax*(th(i)<=0)), cos(thMax*(th(i)>0) - thMax*(th(i)<=0))];
            m = [cos(sinPercent*th(i)+ (1-sinPercent)*thMax*(th(i)>0) - (1-sinPercent)*thMax*(th(i)<=0)), -sin(sinPercent*th(i)+ (1-sinPercent)*thMax*thMax*(th(i)>0) - (1-sinPercent)*thMax*(th(i)<=0));...
                sin(sinPercent*th(i)+ (1-sinPercent)*thMax*thMax*(th(i)>0) - (1-sinPercent)*thMax*(th(i)<=0)), cos(sinPercent*th(i)+ (1-sinPercent)*thMax*thMax*(th(i)>0) - (1-sinPercent)*thMax*(th(i)<=0))];
            %m = [cos(-thMax*(th(i)>0)), -sin(-thMax*(th(i)>0)); sin(-thMax*(th(i)>0)), cos(-thMax*(th(i)>0))]; % le rotation
            %thMax*(th(1:totFrames)>0) - (thMax*th(1:totFrames)<=0)
            %m = [cos(0) -sin(0); sin(0) cos(0)];
            %m = [1+thMax*sin(th(i)) 0; 0 1+thMax*sin(th(i))];
        end
        
        % this is a bit subtle - get mo vector from mouse movement, but
        % apply transformed version of it to cursor
        if origAtCurs
%             tFormedCoords = m*[x-xOld; y-yOld];   % transform coordinates!
%             tFormedCoords = tFormedCoords + [xOld; yOld];
            deltaCoords = m*[x-xOld; y-yOld];   % transform coordinates!
            tFormedCoords = tFormedCoords + deltaCoords;

        else
            tFormedCoords = m*[x-centerX; y-centerY];   % transform coordinates!
            tFormedCoords = tFormedCoords + [centerX; centerY];
        end
        
        % target is drawn first, so cursor will be visible on top of it
        tempCoords = [targCoords(1,i), tFormedCoords(1,1) ; targCoords(2,i), tFormedCoords(2,1)];
        
                %line below makes this the 'continuous cursor' condition,
        %commenting it out will make it refresh with target instead
        %Additionally, using 2 as column number makes the cursor continuous
        %but using 1 as a column number makes the target continuous instead
        repeatCoords(:,1) = tempCoords(:,1);
        if mod(i,round(frameRate/updateRate)) == 0
           repeatCoords = tempCoords; 
        end
        Screen('DrawDots', theWindow, repeatCoords, dotSizes, dotColors);
        
        %Screen('DrawDots', theWindow, tempCoords, dotSizes, dotColors);
        %Screen('Flip', theWindow, GetSecs()+sampleTime);
        vbl = Screen('Flip', theWindow, vbl + (waitframes - 0.5)*ifi);  % , GetSecs()+sampleTime
        
        % DRIFT_X and Y ADDED TO X AND Y 07/29
        respCoords(:,i) = [drift_x + x; drift_y + y];               % where the mouse or finger actually is
        tFormRespCoords(:,i) = tFormedCoords;   % where the cursor is on screen
        txOld = tFormedCoords(1);
        tyOld = tFormedCoords(2);
        
        
    end % end trial for loop

    % Close up
    ShowCursor(0);
    Screen(theWindow,'Close');

catch
    Screen('CloseAll');
    ShowCursor;
    psychrethrow(psychlasterror);
end %try..catch..


%keyboard; 




%% data analysis!
% let's look at x and y correlation peaks and latencies as our index of
% adaptation

% convert to samples for lag and front bumper
maxLag = round(maxLag*frameRate);
frontPorch = round(frontPorch*frameRate);

% compute transformed target coordinates just for the hell of it
tFormedTargCoords = m*[targCoords(1,:)-centerX; targCoords(2,:)-centerY];   % transform coordinates!
tFormedTargCoords = tFormedTargCoords + repmat([centerX; centerY], 1, totFrames);

% center and scale coords to screen width
% respCoords = where the mouse is on the desk (or the finger on the pad)
normRespCoords(1,:) = (respCoords(1,:) - centerX)./ centerX;
normRespCoords(2,:) = (respCoords(2,:) - centerY)./ centerY;

normTargCoords(1,:) = (targCoords(1,:) - centerX)./ centerX;
normTargCoords(2,:) = (targCoords(2,:) - centerY)./ centerY;

% tFormRespCoords = where the cursor is on the screen
normTformRespCoords(1,:) = (tFormRespCoords(1,:) - centerX)./ centerX;
normTformRespCoords(2,:) = (tFormRespCoords(2,:) - centerY)./ centerY;

normTformTargCoords(1,:) = (tFormedTargCoords(1,:) - centerX)./ centerX;
normTformTargCoords(2,:) = (tFormedTargCoords(2,:) - centerY)./ centerY;

% compute velocities (from raw data, tho it shouldn't make a diff -ha!
respVel = diff(respCoords, 1, 2);
targVel = diff(targCoords, 1, 2);

% compute positional and velocity cross correlations
% velocity
[xCorX, xLagsVel] = xcorr(respVel(1,frontPorch:end), targVel(1,frontPorch:end), maxLag, 'unbiased');
[xCorY, yLagsVel] = xcorr(respVel(2,frontPorch:end), targVel(2,frontPorch:end), maxLag, 'unbiased');

% find peak lag and height
[xMax, xLagAtMax] = max(xCorX);
[yMax, yLagAtMax] = max(xCorY);
xRelLag = xLagAtMax-maxLag-1; % lag relative to zero
yRelLag = yLagAtMax-maxLag-1; % lag relative to zero
xMaxLagSecs = xRelLag./frameRate;
yMaxLagSecs = yRelLag./frameRate;


% cursor error
mausErr = normRespCoords(:, (xRelLag+1):end) - normTargCoords(:, 1:(end-xRelLag));
% error re. maus on screen
cursErr = normTformRespCoords(:, (xRelLag+1):end) - normTargCoords(:, 1:(end-xRelLag));

pythagCursErr = sqrt(cursErr(1,:).^2+cursErr(2,:).^2);

cursErrNoLag = normTformRespCoords - normTargCoords;
pythagCursErrNoLag = sqrt(cursErrNoLag(1,:).^2+cursErrNoLag(2,:).^2);

theLen = length(pythagCursErr);
    [fAx, myInds] = fftaxes(theLen, frameRate);
    theMid = round(theLen./2);
    % determine how much of the spectrum to plot in a non-magic way
    stimFreqInd =  find(fAx>3*tfreq, 1);
    stimFreqInd =  find(fAx>1.5, 1);
    notMagicNumber = stimFreqInd - theMid;
    absCursFT = fftshift(abs(fft(killDC(pythagCursErr))));
    
    pythagRespVel = sqrt(respVel(1,:).^2+respVel(2,:).^2);
    
    theLen = length(pythagRespVel);
    [fAx2, myInds2] = fftaxes(theLen, frameRate);
    theMid2 = round(theLen./2);
    % determine how much of the spectrum to plot in a non-magic way
    stimFreqInd =  find(fAx2>3*tfreq, 1);
    stimFreqInd =  find(fAx>1.5, 1);
    notMagicNumber2 = stimFreqInd - theMid2;
    absRespVelFT = fftshift(abs(fft(killDC(pythagRespVel))));
    


%% plotting!!

% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 16)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 16)

    
    
   
    window = 1/tfreq*frameRate/2;
    angles_to_test = [0.1 0.2 0.28 0.35 0.4 0.5];
    num_angles = numel(angles_to_test);
    
    cursErr_snippetAve4_angles_cond = zeros(num_angles,16,4,window);
cursErr_snippet_freqAve4_angles_cond = zeros(num_angles,16,4,window);
VelXErr_snippet_freqAve4_angles_cond = zeros(num_angles,16,4,window);
VelYErr_snippet_freqAve4_angles_cond = zeros(num_angles,16,4,window);

[cursErr_snippetAve_angles_cond] = zeros(num_angles,16,10,window);
[cursErr_snippet_freqAve_angles_cond] = zeros(num_angles,16,10,window);
[cursErr_init_angles_cond] = zeros(num_angles,16,10);
[cursErr_init_freq_angles_cond] = zeros(num_angles,16,10);
    
    for i = 1:numel(angles_to_test)
    [all_inds,sum_inds] = IndsForAngles(0.0625,angles_to_test(i),60,window,trialTime);
    [cursErrMaxAve, CursErrLagAtMaxAve, cursErr_snippetAve, cursErr_freqMaxAve, CursErr_freqLagAtMaxAve, cursErr_snippet_freqAve,...
    cursErr_init,cursErr_freq_init,cursErr_snippetAve4,cursErr_snippet_freqAve4,VelXErr_snippet_freqAve4,VelYErr_snippet_freqAve4] = windowedErrsTrialAve(respVel,targVel,normRespCoords,normTargCoords,normTformRespCoords,all_inds,60,window,tfreq);
    
    [cursErr_snippetAve_angles_cond] = shortened_snippet_placement(1,i,cursErr_snippetAve,cursErr_snippetAve_angles_cond,1,window);
[cursErr_snippet_freqAve_angles_cond] = shortened_snippet_placement(1,i,cursErr_snippet_freqAve,cursErr_snippet_freqAve_angles_cond,1,window);
[cursErr_init_angles_cond] = shortened_snippet_placement(1,i,cursErr_init,cursErr_init_angles_cond,1,window);
[cursErr_init_freq_angles_cond] = shortened_snippet_placement(1,i,cursErr_freq_init,cursErr_init_freq_angles_cond,1,window);

[cursErr_snippetAve4_angles_cond] = shortened_snippet_placement(1,i,cursErr_snippetAve4,cursErr_snippetAve4_angles_cond,4,window);
[cursErr_snippet_freqAve4_angles_cond] = shortened_snippet_placement(1,i,cursErr_snippet_freqAve4,cursErr_snippet_freqAve4_angles_cond,4,window);
[VelXErr_snippet_freqAve4_angles_cond] = shortened_snippet_placement(1,i,VelXErr_snippet_freqAve4,VelXErr_snippet_freqAve4_angles_cond,4,window);
[VelYErr_snippet_freqAve4_angles_cond] = shortened_snippet_placement(1,i,VelYErr_snippet_freqAve4,VelYErr_snippet_freqAve4_angles_cond,4,window);
    
    end
    
    cursErr_snippetAve4_angles_cond = squeeze(cursErr_snippetAve4_angles_cond(:,1,:,:));
    
    
    
      for j = 1:num_angles
        for k = 1:4
          [cursErr_freq4MaxCond(j,k), cursErr_freq4LagAtMaxCond(j,k)] = max(squeeze(cursErr_snippetAve4_angles_cond(j,k,:)));    
        end
      end
      

    
      
      for i = 1:numel(angles_to_test)
[angle_inds,sum_inds] = IndsForAngles(tfreq,angles_to_test(i),frameRate,window,trialTime);
    all_first_inds(i,:) = angle_inds(1:4);
      end
      
      
      WhereMaxLag = cursErr_freq4LagAtMaxCond(:,1:4) + all_first_inds;


angle_at_Max = th(WhereMaxLag);


%other_th = pi/6*round(th(1:totFrames)/thMax);
other_th = thMax*(th(1:totFrames)>0) - thMax*(th(1:totFrames)<=0);
numFramesToPlot = 6000;


positive_inds = find(other_th > 0);
negative_inds = find(other_th <= 0);
verid_inds = find(other_th == 0);

positive_inds(positive_inds>numel(pythagCursErr)) = [];
negative_inds(negative_inds>numel(pythagCursErr)) = [];
verid_inds(verid_inds>numel(pythagCursErr)) = [];


[R,P]=corrcoef(cursErr(:,positive_inds)');



[R,P]=corrcoef(cursErr(:,negative_inds)');



[R,P]=corrcoef(cursErr(:,verid_inds)');




something = diff(positive_inds);
positive_change = find(something~=1); %length of positive section
positive_change_vals = something(positive_change); %size of jump to next positive section

something2 = diff(negative_inds);
negative_change = find(something2~=1); %length of positive section
negative_change_vals = something2(negative_change); %size of jump to next positive section

numDivisions = floor(positive_change(1)/(frameRate/2));

%this is only true for the current design
for i = 2:(numel(positive_change)-1)
windowAroundPositiveChangeErr(i,:) = pythagCursErr(((i-1)*2*positive_change-window/2+positive_inds(1)):((i-1)*2*positive_change+window/2+positive_inds(1)));
end
for i = 2:(numel(negative_change)-1)
windowAroundNegativeChangeErr(i,:) = pythagCursErr(((i-1)*2*negative_change-window/2+negative_inds(1)):((i-1)*2*negative_change+window/2+negative_inds(1)));
end



clean_divide = floor(numel(positive_inds)/positive_change(1));
remainder = rem(numel(positive_inds),positive_change(1));
another_idea = reshape(positive_inds(1:(end-remainder)),positive_change(1),clean_divide);
another_idea_trunc = another_idea(1:(numDivisions*frameRate/2),:);
another_idea_trunc = reshape(another_idea_trunc,frameRate/2,numel(another_idea_trunc)/(frameRate/2));
error_x = cursErr(1,another_idea_trunc);
error_x = reshape(error_x,(numDivisions*frameRate/2),numel(error_x)/(numDivisions*frameRate/2));
error_y = cursErr(2,another_idea_trunc);
error_y = reshape(error_y,(numDivisions*frameRate/2),numel(error_y)/(numDivisions*frameRate/2));
error_x_pos_ave = mean(error_x,2);
error_y_pos_ave = mean(error_y,2);
error_x_pos_abs_ave = mean(abs(error_x),2);
error_y_pos_abs_ave = mean(abs(error_y),2);


positive_divisions = [];
for i = 1:(size(another_idea_trunc,2)/numDivisions)
   positive_divisions((end+1):(end+frameRate/2),:) = another_idea_trunc(:,((i-1)*numDivisions+1):(i*numDivisions)); 
end


for i = 1:numDivisions
    [b,bint,r,rint,stats] = regress(cursErr(2,positive_divisions(:,i))',[cursErr(1,positive_divisions(:,i))' ones(size(cursErr(1,positive_divisions(:,i))'))]);
    slope_pos(i) = b(1);
    intercept_pos(i) = b(2);
end


clean_divide = floor(numel(negative_inds)/positive_change(1));
remainder = rem(numel(negative_inds),positive_change(1));
another_idea = reshape(negative_inds(1:end-remainder),positive_change(1),(numel(negative_inds)-remainder)/positive_change(1));
another_idea_trunc = another_idea(1:(numDivisions*frameRate/2),:);
another_idea_trunc = reshape(another_idea_trunc,frameRate/2,numel(another_idea_trunc)/(frameRate/2));
error_x = cursErr(1,another_idea_trunc);
error_x = reshape(error_x,(numDivisions*frameRate/2),numel(error_x)/(numDivisions*frameRate/2));
error_y = cursErr(2,another_idea_trunc);
error_y = reshape(error_y,(numDivisions*frameRate/2),numel(error_y)/(numDivisions*frameRate/2));
error_x_neg_ave = mean(error_x,2);
error_y_neg_ave = mean(error_y,2);
error_x_neg_abs_ave = mean(abs(error_x),2);
error_y_neg_abs_ave = mean(abs(error_y),2);


negative_divisions = [];
for i = 1:(size(another_idea_trunc,2)/numDivisions)
   negative_divisions((end+1):(end+frameRate/2),:) = another_idea_trunc(:,((i-1)*numDivisions+1):(i*numDivisions)); 
end


for i = 1:numDivisions
    [b,bint,r,rint,stats] = regress(cursErr(2,negative_divisions(:,i))',[cursErr(1,negative_divisions(:,i))' ones(size(cursErr(1,negative_divisions(:,i))'))]);
    slope_neg(i) = b(1);
    intercept_neg(i) = b(2);
end
%all_sections

mean_NegChangeErr = mean(windowAroundNegativeChangeErr,1);
mean_PosChangeErr = mean(windowAroundPositiveChangeErr,1);


save([initials '_' num2str(tfreq) 'pi_4Sine' num2str(sinPercent) 'ThetaStep' num2str(stepSize) 'Conv' num2str(numConv) '_FBTarg' num2str(updateRate) '.mat']);