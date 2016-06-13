%pass in all necessary figure pointers and plot positions to create plots
%in the master script - analysisFigs.m

function [] = brownianPointingSim(verid_delta_scalar,rot_delta_scalar,master_bDuFH, master_bDuFV, master_xcor,master_xCorrPos,nRows,nColumns,plot_pos)



% BrownianPointing - 
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

% NOTES
%   4/28    - need to implement sinusoidal variation
%           & *really* need to code a file saving & analysis pipeline
%   4/29    - I think a cool condition would be to have the target jump a
%           relatively long way, and sit there until until the subject's
%           responds "stops" by come criteria.
%   4/29    - Yummy!.
%   6/17    Using the old BrownianPointing to set up graphs for Methods
%               section

%   6/19    Using brownianpointingMethods to simulate learning and set up
%               graphs for hypotheses


% here are the input parameters:
%verid_delta_scalar - magnitude (usually between 0 and 1) for veridical and
      %hybrid models
%rot_delta_scalar - magnitude (usually between 0 and 1) for rotation and
      %hybrid models
%master_bDuFH - pointer to the horizontal FFT plot grid 
%master_bDuFV - pointer to the vertical FFT plot grid 
%master_xcor - pointer to the xcorr velocity plot grid
%master_xCorrPos - pointer to xcorr position plot grid
%nRows - number of rows in the plot grids (depending on # models)
%nColumns - number of columns in the plot grids (depending on # models)
%plot_pos - designated plot position for this particular run

% There are no outputs because all graphs have been created elsewhere


%% set up initial parameters
trialTime = 20;      % trial duration in s
stepSize = 5;      % ave or max stepsize in pixels - determines target speed
frameRate = 60;     % will determine this "live" in future versions
frameRateDiscrete = 12;
stepSizeDiscrete = stepSize*frameRate./frameRateDiscrete;
totFrames = trialTime*frameRate;    % yup
totFramesDiscrete = trialTime*frameRateDiscrete;
tfreq = 0.5;  % temp. freq. of parameter variation
trunc = 600; % graphing param for methods

% these "maximums" are the extreme vals of the parameters in the dynamic
% case, and are simply the parameters in the static case
aMax = 0.6;
bMax = 0.6;
thMax = 5*pi/12; 

% flags for stuff that will become arguments or gui checkboxes at some
% point...
mouseVisibleFlag = 1;   % visual feedback of whey yo finger at
dynoUpdateFlag = 1;     % dynamic updating of transform
scalerFlag = 0;        % scaler or rotation transformation?
origAtCurs = 1;         % do transform around cursor or screen center?
dynamic_mag_flag = 0;  % 0 means magnitude is constant, 1 magnitude changes

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
        th2 = thMax.*sin(2.*pi.*2*tfreq.*(1/frameRate).*(1:totFrames));
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
temp = randn(2,totFrames);
temp(1,:) = conv(temp(1,:), [.25,.25,.25,.25], 'same');
temp(2,:) = conv(temp(2,:), [.25,.25,.25,.25], 'same');
temp(:,1) = [0; 0];
targCoords = stepSize*cumsum(temp, 2);
temp2 = randn(2,totFramesDiscrete);
temp2(1,:) = conv(temp2(1,:), [.25,.25,.25,.25], 'same');
temp2(2,:) = conv(temp2(2,:), [.25,.25,.25,.25], 'same');
temp2(:,1) = [0; 0];
targCoords2 = stepSizeDiscrete*cumsum(temp2, 2);


respCoords = zeros(size(targCoords));   % response coordinate array
tFormRespCoords = respCoords;

%% stimulus presentation
try
    
    % Open up a window on the screen and clear it.
    whichScreen = max(Screen('Screens'));
    [theWindow,theRect] = Screen(whichScreen,'OpenWindow',0,[],[],2);

    % Move the cursor to the center of the screen
    centerX = theRect(RectRight)/2;
    centerY = theRect(RectBottom)/2;
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
    tFormedCoords = [centerX; centerY];
    txOld = x;
    tyOld = y;
    Screen('DrawDots', theWindow, [centerX, centerX; centerY, centerY], dotSizes, dotColors);
    Screen('Flip', theWindow);
    sampleTime = 0.01;
    startTime = GetSecs;
    nextTime = startTime+sampleTime;
    
    pause(0.5); % ready... set...
    
    % some of these may not be necessary
    opp_m_super_old = [cos(0), -sin(0); sin(0), cos(0)];
    opp_m_old = [cos(0), -sin(0); sin(0), cos(0)];
    opposite_m = [cos(0), -sin(0); sin(0), cos(0)];
    
    for i = 1:totFrames
        %dynamic mag is the 'stiffness' option
        if dynamic_mag_flag
            dynamic_mag = 1./thMax.*(thMax-abs(th(i))); %stiffness for large perturbations
        else
            dynamic_mag = 1; % all is equal
        end
        xOld = x;
        yOld = y;
        opp_m_super_old = opp_m_old;
        opp_m_old = opposite_m;
        
        
        % next line commented out because we create simulated input
        %[x,y,buttons] = GetMouse(whichScreen);
        if i < 3 % I use this because tFormRespCoords don't exist at t = 1
          verid_deltaCoords = [(targCoords(1,i)-tFormedCoords(1,1)) (targCoords(2,i)-tFormedCoords(2,1))];
          opp_deltaCoords = opposite_m*[0;0];
          %opp_tFormedCoords = tFormedCoords;
        else
            % need to compare to tFormResCoords(i-1) if using targCoords(i)
          
            %this creates the change necessary for the opposite
            %perturbation in the 'rotation' and 'hybrid' models
            opp_deltaCoords = opposite_m*[(targCoords(1,i)-tFormRespCoords(1,i-1)); (targCoords(2,i)-tFormRespCoords(2,i-1))];
            verid_deltaCoords = [(targCoords(1,i)-tFormedCoords(1,1)) (targCoords(2,i)-tFormedCoords(2,1))];
        end
            % the following junk code is sticking around just in case
            % + opp_m_old*[0.3*(targCoords(1,i)-tFormRespCoords(1,i-1)); 0.3*(targCoords(2,i)-tFormRespCoords(2,i-1))]...
                %+opp_m_super_old*[0.3*(targCoords(1,i)-tFormRespCoords(1,i-1)); 0.3*(targCoords(2,i)-tFormRespCoords(2,i-1))];%[x-xOld; y-yOld];   % transform coordinates!
            %opp_tFormedCoords = tFormedCoords + opp_deltaCoords;

            %x and y created by adding the combination of rotation and
            %veridical target-tracking residuals to the current position,
            %with the option of providing 'stiffness' during the larger
            %perturbations or keeping the magnitude even across.
            x = tFormedCoords(1) + dynamic_mag*rot_delta_scalar*opp_deltaCoords(1) + dynamic_mag*verid_delta_scalar*verid_deltaCoords(1);
            y = tFormedCoords(2) + dynamic_mag*rot_delta_scalar*opp_deltaCoords(2) + dynamic_mag*verid_delta_scalar*verid_deltaCoords(2); 

        %end
        
        if scalerFlag   % scaler tranformation
            m = [a, b(i); c(i), d];
        else    % rotation
            m = [cos(th(i)), -sin(th(i)); sin(th(i)), cos(th(i))]; % le rotation
            opposite_m = [cos(-th(i)), -sin(-th(i)); sin(-th(i)), cos(-th(i))]; % le rotation
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
        Screen('DrawDots', theWindow, tempCoords, dotSizes, dotColors);
        Screen('Flip', theWindow, GetSecs()+sampleTime);
        
        respCoords(:,i) = [x; y];               % where the mouse or finger actually is
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

% added for the methods figs
normTargCoords2(1,:) = (targCoords2(1,:) - centerX)./ centerX;
normTargCoords2(2,:) = (targCoords2(2,:) - centerY)./ centerY;

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

% (normalized) position
[xCorXPos, xLagsPos] = xcorr(normTformRespCoords(1,frontPorch:end), normTargCoords(1,frontPorch:end), maxLag, 'unbiased');
[xCorYPos, yLagsPos] = xcorr(normTformRespCoords(2,frontPorch:end), normTargCoords(2,frontPorch:end), maxLag, 'unbiased');

% find peak lag and height
[xMax, xLagAtMax] = max(xCorX);
[yMax, yLagAtMax] = max(xCorY);
xRelLag = xLagAtMax-maxLag-1; % lag relative to zero
yRelLag = yLagAtMax-maxLag-1; % lag relative to zero
xMaxLagSecs = xRelLag./frameRate;
yMaxLagSecs = yRelLag./frameRate;

if xRelLag < 0 
    xRelLag = 0; 
    fprintf('Warning: x lag less than zero; S has ESP.');
end
if yRelLag < 0 
    yRelLag = 0; 
    fprintf('Warning: y lag less than zero; S has ESP.');
end

% compute running diffs between mouse/cursor and target position after
% subtracting the response lag.  The idea is to treat each target movement
% as a "trial", and look at the associated response error on each of these
% trials.  rad.
% cursor error
mausErr = normRespCoords(:, (xRelLag+1):end) - normTargCoords(:, 1:(end-xRelLag));
% error re. maus on screen
cursErr = normTformRespCoords(:, (xRelLag+1):end) - normTargCoords(:, 1:(end-xRelLag));

%% spatial cross correlations
% okay, now lets see if we can do some sort of correlation on the spatial
% distributions of the target vs. cursor and target vs. mouse data...
% The general plan is to make spatial arrays that are all zeros except for 
% where the target, etc. went.

% make matrices of the right size (based on max excursion during that
% trial?)  We'll want to be working in pixel space.
% tFormRespCoords, respCoords, and targCoords are the relevant data

% find overall max and min of the relevant data
maxCoords = ceil(max([tFormRespCoords, respCoords, targCoords], [], 2));
minCoords = floor(min([tFormRespCoords, respCoords, targCoords], [], 2));
xArrSize = maxCoords(1)-minCoords(1)+1;
yArrSize = maxCoords(2)-minCoords(2)+1;

% init arrays
targMat = zeros(yArrSize, xArrSize);
respMat = zeros(yArrSize, xArrSize);
cursMat = zeros(yArrSize, xArrSize);

% set pixels on a path to 1
ttemp = ceil(targCoords) - repmat(minCoords, 1, totFrames)+1;
tinds = sub2ind(size(targMat), ttemp(2,:), ttemp(1,:));
targMat(round(tinds)) = 1;

rtemp = ceil(respCoords) - repmat(minCoords, 1, totFrames)+1;
rinds = sub2ind(size(respMat), rtemp(2,:), rtemp(1,:));
respMat(round(rinds)) = 1;

ctemp = ceil(tFormRespCoords) - repmat(minCoords, 1, totFrames)+1;
cinds = sub2ind(size(cursMat), ctemp(2,:), ctemp(1,:));
cursMat(round(cinds)) = 1;

% blur a bit
mykern = ggaus(10, 1);
targMat = conv2(targMat, mykern, 'same');
respMat = conv2(respMat, mykern, 'same');
cursMat = conv2(cursMat, mykern, 'same');

% compute cross correlations
targAndResp = xcorr2(targMat, respMat);
targAndCurs = xcorr2(targMat, cursMat);


%% plotting!!

% ONLY PLOT XCORRS AND FFTS

%{

% raw data traces
rawFig = figure;
    plot(normTargCoords(1,:),normTargCoords(2,:), 'go:');
    hold on;
    plot(normTformRespCoords(1,:),normTformRespCoords(2,:), 'co:');
    plot(normRespCoords(1,:),normRespCoords(2,:), 'ko-');
    legend('target on screen', 'cursor on screen', 'finger/mouse location');
    % mark the beginning and end points
    plot(normTargCoords(1,1),normTargCoords(2,1), 'b*', 'MarkerSize', 20, ...
                                                        'MarkerFaceColor', 'b');
    plot(normTargCoords(1,end),normTargCoords(2,end), 'rh', 'MarkerSize', 20, ...
                                                        'MarkerFaceColor', 'r');

    hold off;
    axis equal;
    
    % This following plot is for the methods section only, raw data traces
    tracesMethodsFig = figure;
    subplot(1,2,1);
    plot(normTargCoords(1,:),normTargCoords(2,:), 'go:');
    hold on;
    legend('continuous target trace');
    % mark the beginning and end points
    plot(normTargCoords(1,1),normTargCoords(2,1), 'b*', 'MarkerSize', 20, ...
                                                        'MarkerFaceColor', 'b');
    plot(normTargCoords(1,end),normTargCoords(2,end), 'rh', 'MarkerSize', 20, ...
                                                        'MarkerFaceColor', 'r');

    hold off;
    axis equal;
    subplot(1,2,2);
    
    plot(normTargCoords2(1,:),normTargCoords2(2,:), 'go:');
    hold on;
    legend('discrete target trace');
    % mark the beginning and end points
    plot(normTargCoords2(1,1),normTargCoords2(2,1), 'b*', 'MarkerSize', 20, ...
                                                        'MarkerFaceColor', 'b');
    plot(normTargCoords2(1,end),normTargCoords2(2,end), 'rh', 'MarkerSize', 20, ...
                                                        'MarkerFaceColor', 'r');

    hold off;
    axis equal;

% plots of the positions and velocities as a function of time
% The instances of 1:end below can be replace by frontPorch:end to omit
% plotting the first n=frontPorch samples
tracesFig = figure;
subplot(2,2,1); % H pos
    plot(normRespCoords(1,1:end), 'k')
    hold on
    plot(normTformRespCoords(1,1:end), 'c')
    plot(normTargCoords(1,1:end), 'g')
    legend('mouse', 'cursor', 'target');
    title('Horizontal Position');
    hold off;
subplot(2,2,2); % V pos
    plot(normRespCoords(2,1:end), 'k')
    hold on
    plot(normTformRespCoords(2,1:end), 'c')
    plot(normTargCoords(2,1:end), 'g')
    legend('mouse', 'cursor', 'target');
    title('Vertical Position');
    hold off;
subplot(2,2,3); % H vel
    plot(targVel(1,1:end), 'g')
    hold on
   plot(respVel(1,1:end), 'k')
    legend('response', 'target');
    title('Horizontal Velocity');
    hold off;
subplot(2,2,4); % V vel
    plot(targVel(2,1:end), 'g')
    hold on
    plot(respVel(2,1:end), 'k')
    legend('response', 'target');
    title('Vertical Velocity');
    hold off;

%}

% plots of the cross-correlations
%xcorFig = figure;
figure(master_xcor);
subplot(nRows,nColumns,plot_pos);
%{
subplot(1, 2, 1);   % position
    plot(xLagsPos./frameRate, xCorXPos);
    hold on; 
    plot(yLagsPos./frameRate, xCorYPos, 'k'); 
    legend('horizontal position', 'vertical position');
    xlabel('lag (s)');
    hold off;
subplot(1, 2, 2);   % velocity
%}
    plot(xLagsVel./frameRate, xCorX);
    hold on; 
    plot(yLagsVel./frameRate, xCorY, 'k'); 
    myylim = ylim;
    line([xMaxLagSecs, xMaxLagSecs], [myylim(1), myylim(2)], ...
                                        'Color', 'b', ...
                                        'LineStyle', ':');
    line([yMaxLagSecs, yMaxLagSecs], [myylim(1), myylim(2)], ...
                                        'Color', 'k', ...
                                        'LineStyle', ':');
    legend('horizontal velocity', 'vertical velocity', 'Location', 'NorthWest');
    text(xMaxLagSecs, xMax, [num2str(xMaxLagSecs, '%1.3f'), '  ', num2str(xMax, '%1.2f')]);
    text(yMaxLagSecs, yMax, [num2str(yMaxLagSecs, '%1.3f'), '  ', num2str(yMax, '%1.2f')]);
    line([-1 1],[0 0]);
    xlabel('lag (s)');
    hold off;

    
    %{
    
% plot RT (lag) shifted errors
lagErrFig = figure;
subplot(2,2,1); % H pos, lag shifted
    plot(normRespCoords(1,(xRelLag+1):end), 'k')
    hold on
    plot(normTformRespCoords(1,(xRelLag+1):end), 'c')
    plot(normTargCoords(1,1:(end-xRelLag)), 'g')
    legend('mouse', 'cursor', 'target');
    title('Horizontal Position');
    hold off;
subplot(2,2,2); % V pos, lag shifted
    plot(normRespCoords(2,(xRelLag+1):end), 'k')
    hold on
    plot(normTformRespCoords(2,(xRelLag+1):end), 'c')
    plot(normTargCoords(2,1:(end-xRelLag)), 'g')
    legend('mouse', 'cursor', 'target');
    title('Vertical Position');
    hold off;
subplot(2,2,3); % H pos error
    plot(mausErr(1,1:end), 'k')
    hold on
    plot(cursErr(1,1:end), 'c')
    mylim = xlim;
    line([mylim(1), mylim(2)], [0, 0], ...
        'Color', 'g', 'LineStyle', ':');
    legend('mouse error', 'cursor error');
    title('Horizontal Position');
    hold off;
subplot(2,2,4); % V pos error
    plot(mausErr(2,1:end), 'k')
    hold on
    plot(cursErr(2,1:end), 'c')
    mylim = xlim;
    line([mylim(1), mylim(2)], [0, 0], ...
        'Color', 'g', 'LineStyle', ':');
    legend('mouse error', 'cursor error');
    title('Vertical Position');
    hold off;

    
    
% NEED TO FIX FOR SCALER TRANSFORMS
% plot error against mag. of tranform
errWithTformFig = figure;
subplot(1,2,1) % H err
    plot(mausErr(1,1:end), 'k')
    hold on
    plot(cursErr(1,1:end), 'c')
    mylim = xlim;
    line([mylim(1), mylim(2)], [0, 0], ...
        'Color', 'g', 'LineStyle', ':');
    mylim = ylim;
    scFac = 0.5*(mylim(2)./thMax);
    plot(scFac.*th, 'r:');
    legend('mouse error', 'cursor error', 'tForm');
    title('Horizontal Position');
    hold off;
subplot(1,2,2); % V pos error
    plot(mausErr(2,1:end), 'k')
    hold on
    plot(cursErr(2,1:end), 'c')
    mylim = xlim;
    line([mylim(1), mylim(2)], [0, 0], ...
        'Color', 'g', 'LineStyle', ':');
    mylim = ylim;
    scFac = 0.5*(mylim(2)./thMax);
    plot(scFac.*th, 'r:');
    legend('mouse error', 'cursor error', 'tForm');
    title('Vertical Position');
    hold off;
    
    % another plot for the methods section
    th_disc = NaN(size(th));
    th2_disc = NaN(size(th2));
    th_disc(1:(frameRate/frameRateDiscrete):end) = th(1:(frameRate/frameRateDiscrete):end);
    th2_disc(1:(frameRate/frameRateDiscrete):end) = th2(1:(frameRate/frameRateDiscrete):end);
    frequencyFig = figure;
    subplot(2,2,1);
    mylim = ylim;
    scFac = 0.5*(mylim(2)./thMax);
    plot(scFac.*th(1:trunc), 'ro:');
    hold on;
    plot(-scFac.*th(1:trunc), 'bo:');
    legend('perturbation','correct adaptive response');
    xlabel('Trial');
    ylabel('Perturbed Response');
    hold off;
    
        subplot(2,2,2);
    mylim = ylim;
    scFac = 0.5*(mylim(2)./thMax);
    plot(scFac.*th2(1:trunc), 'ro:');
        hold on;
        plot(-scFac.*th2(1:trunc), 'bo:');
        legend('perturbation','correct adaptive response');
            xlabel('Trial');
    ylabel('Perturbed Response');
    hold off;
    
        subplot(2,2,3);
    mylim = ylim;
    scFac = 0.5*(mylim(2)./thMax);
    plot(scFac.*th_disc(1:trunc), 'ro:');
        hold on;
        plot(-scFac.*th_disc(1:trunc), 'bo:');
        legend('perturbation','correct adaptive response');
            xlabel('Trial');
    ylabel('Perturbed Response');
    hold off;
    
        subplot(2,2,4);
    mylim = ylim;
    scFac = 0.5*(mylim(2)./thMax);
    plot(scFac.*th2_disc(1:trunc), 'ro:');
        hold on;
        plot(-scFac.*th2_disc(1:trunc), 'bo:');
        legend('perturbation','correct adaptive response');
            xlabel('Trial');
    ylabel('Perturbed Response');
    hold off;
    
  %}  
    
    
% DITTO
% Fourier plots of above - see if the error is frequency matched (or
% related) to the magnitude of the transform
%baronDuFourierH = figure;
figure(master_bDuFH);
subplot(nRows,nColumns,plot_pos);
%subplot(1,2,1) % H err
    magicNumber = 40;
    theLen = length(cursErr(1,:));
    [fAx, myInds] = fftaxes(theLen, frameRate);
    theMid = round(theLen./2);
    hCursFT = normalize(fftshift(abs(fft(killDC(cursErr(1,:))))));
    plot(fAx(theMid-magicNumber:theMid+magicNumber), hCursFT(theMid-magicNumber:theMid+magicNumber));
    [peaks,locs] = findpeaks(hCursFT(theMid-magicNumber:theMid+magicNumber),'sortstr','descend');
    hold on
    
        transFT = normalize(fftshift(abs(fft(killDC(th(1:length(cursErr(1,:))))))));
    plot(fAx(theMid-magicNumber:theMid+magicNumber), transFT(theMid-magicNumber:theMid+magicNumber), 'r');
    
    line([fAx(theMid-magicNumber+locs(1)-1), fAx(theMid-magicNumber+locs(1)-1)], [myylim(1), myylim(2)], ...
                                        'Color', 'b', ...
                                        'LineStyle', ':');
    
    text(fAx(theMid-magicNumber+locs(1)-1), peaks(1), num2str(fAx(theMid-magicNumber+locs(1)-1), '%1.3f'));
    
        line([fAx(theMid-magicNumber+locs(3)-1), fAx(theMid-magicNumber+locs(3)-1)], [myylim(1), myylim(2)], ...
                                        'Color', 'b', ...
                                        'LineStyle', ':');
    
    text(fAx(theMid-magicNumber+locs(3)-1), peaks(3), num2str(fAx(theMid-magicNumber+locs(3)-1), '%1.3f'));
    
    
        line([fAx(theMid-magicNumber+locs(5)-1), fAx(theMid-magicNumber+locs(5)-1)], [myylim(1), myylim(2)], ...
                                        'Color', 'b', ...
                                        'LineStyle', ':');
    
    text(fAx(theMid-magicNumber+locs(5)-1), peaks(5), num2str(fAx(theMid-magicNumber+locs(5)-1), '%1.3f'));
    
            line([fAx(theMid-magicNumber+locs(7)-1), fAx(theMid-magicNumber+locs(7)-1)], [myylim(1), myylim(2)], ...
                                        'Color', 'b', ...
                                        'LineStyle', ':');
    
    text(fAx(theMid-magicNumber+locs(7)-1), peaks(7), num2str(fAx(theMid-magicNumber+locs(7)-1), '%1.3f'));
    
            line([fAx(theMid-magicNumber+locs(9)-1), fAx(theMid-magicNumber+locs(9)-1)], [myylim(1), myylim(2)], ...
                                        'Color', 'b', ...
                                        'LineStyle', ':');
    
    text(fAx(theMid-magicNumber+locs(9)-1), peaks(9), num2str(fAx(theMid-magicNumber+locs(9)-1), '%1.3f'));
                                    

    xlabel('Hz');
    ylabel('Signal Magnitude');
    legend('cursor error', 'perturbation function');
    title('Horizontal Frequency');
    axis([-2 2 0 1.2])
    hold off;
%subplot(1,2,2); % V pos error
%baronDuFourierV = figure;
figure(master_bDuFV);

baronDuFourierV = subplot(nRows,nColumns,plot_pos);
    theLen = length(cursErr(1,:));
    [fAx, myInds] = fftaxes(theLen, frameRate);
    theMid = round(theLen./2);
    vCursFT = normalize(fftshift(abs(fft(killDC(cursErr(2,:))))));
    plot(fAx(theMid-magicNumber:theMid+magicNumber), vCursFT(theMid-magicNumber:theMid+magicNumber));
    [peaks,locs] = findpeaks(vCursFT(theMid-magicNumber:theMid+magicNumber),'sortstr','descend');
    hold on
    
        transFT = normalize(fftshift(abs(fft(killDC(th(1:length(cursErr(2,:))))))));
    plot(fAx(theMid-magicNumber:theMid+magicNumber), transFT(theMid-magicNumber:theMid+magicNumber), 'r');
    
        line([fAx(theMid-magicNumber+locs(1)-1), fAx(theMid-magicNumber+locs(1)-1)], [myylim(1), myylim(2)], ...
                                        'Color', 'b', ...
                                        'LineStyle', ':');
    
    text(fAx(theMid-magicNumber+locs(1)-1), peaks(1), num2str(fAx(theMid-magicNumber+locs(1)-1), '%1.3f'));
    
        line([fAx(theMid-magicNumber+locs(3)-1), fAx(theMid-magicNumber+locs(3)-1)], [myylim(1), myylim(2)], ...
                                        'Color', 'b', ...
                                        'LineStyle', ':');
    
    text(fAx(theMid-magicNumber+locs(3)-1), peaks(3), num2str(fAx(theMid-magicNumber+locs(3)-1), '%1.3f'));
    
    
        line([fAx(theMid-magicNumber+locs(5)-1), fAx(theMid-magicNumber+locs(5)-1)], [myylim(1), myylim(2)], ...
                                        'Color', 'b', ...
                                        'LineStyle', ':');
    
    text(fAx(theMid-magicNumber+locs(5)-1), peaks(5), num2str(fAx(theMid-magicNumber+locs(5)-1), '%1.3f'));
    
            line([fAx(theMid-magicNumber+locs(7)-1), fAx(theMid-magicNumber+locs(7)-1)], [myylim(1), myylim(2)], ...
                                        'Color', 'b', ...
                                        'LineStyle', ':');
    
    text(fAx(theMid-magicNumber+locs(7)-1), peaks(7), num2str(fAx(theMid-magicNumber+locs(7)-1), '%1.3f'));
    
            line([fAx(theMid-magicNumber+locs(9)-1), fAx(theMid-magicNumber+locs(9)-1)], [myylim(1), myylim(2)], ...
                                        'Color', 'b', ...
                                        'LineStyle', ':');
    
    text(fAx(theMid-magicNumber+locs(9)-1), peaks(9), num2str(fAx(theMid-magicNumber+locs(9)-1), '%1.3f'));
    

        xlabel('Hz');
    ylabel('Signal Magnitude');
    legend('cursor error', 'perturbation function');
    title('Vertical frequency');
    axis([-2 2 0 1.2])
    hold off;


    figure(master_xCorrPos);

subplot(nRows,nColumns,plot_pos);
    
        imagesc(targAndCurs);
    title('target and Cursor');
    
    
 %{   
    
% show image version for 2D cross-correlations
imFig = figure;
subplot(1, 3, 1);
    imagesc(targMat);
    title('target');
subplot(1, 3, 2);
    imagesc(respMat);
    title('mouse/finger')
subplot(1, 3, 3);
    imagesc(cursMat);
    title('cursor');
    colormap jet
    
TwoDXcorrFig = figure;
subplot(1,2,1);
    imagesc(targAndResp);
    title('target and mouse');
subplot(1,2,2);
    imagesc(targAndCurs);
    title('target and Cursor');


    
figure(TwoDXcorrFig);
figure(imFig);
figure(errWithTformFig);
figure(lagErrFig);
figure(tracesFig);
figure(rawFig);
figure(xcorFig);

%}