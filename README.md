# BrownianPointing
my grad school project on visuomotor coordination (think hand-eye coordination)

Background Motivation

I wanted to see the flexibility of the mapping between what the eye sees and what the hand does when normal rules of hand-eye coordination break down. To study this, Larry Cormack and I developed a new experimental paradigm.

In our studies, participants continuously tracked a randomly moving target on screen with a computer mouse while cursor position was transformed through multiple condition-dependent rotation angles. The target changed position in a Brownian manner and never went offscreen. At each timepoint, the change in mouse position was multiplied by a rotation transformation matrix, resulting in a constant angular rotation being applied to the mouse movements at 60 Hz with the previous timepoint’s cursor position as the rotation origin.


BrownianPointing.m

The original file from which the rest descend. Consists of a trial where you use your cursor to track an on-screen square moving with Brownian motion. The cursor may undergo perturbations or transformations, based upon specific built-in settings. After the trial, the visual window closes and basic analyses are plotted in new graphical windows.

 Larry Cormack adapted this first m-file from MouseTraceDemo.m in the Psychtoolbox library for MATLAB. At each timepoint, the change in mouse position was multiplied by a rotation transformation matrix, resulting in a constant angular rotation being applied to the mouse movements at 60 Hz with the previous timepoint’s cursor position as the rotation origin.

He created traces of the mouse, cursor, and target, sampled at 60 hz. Here, he included cross-correlations of the x and y position of mouse-target, cursor-target. From that, we get relative lag between the target and the human response, from which we can calculate mouse-target and cursor-target spatial error. Lastly, we can look at error in Fourier space.

My contribution to the original document: updated measurement strategy for x and y-position to handle times when the cursor went off-screen.


BrownianTrialPerturb.m and BrownianTrialPerturbWrapper.m

I changed BrownianPointing so that I could input different parameters and run multiple trials across different conditions. I also set up a way to save these to a file with names that reflected the condition. BrownianTrialPerturbWrapper is just one example of a wrapper function I could use to run multiple conditions on the same person.

The problem space was pretty unconstrained at this point and I was experimenting with multiple parameters:
tfreq - frequency of change in rotation angle (e.g. the frequency of a sinusoidal change in angle with an amplitude of 45 degrees)
stepSize - represents target speed. Especially relevant with a low target update rate.
sinPercent - there were two original conditions here (1) pure sinusoidal change in rotation angle and (2) alternation between two rotation angles. sinPercent represented the weighted average of these two conditions.
updateRate - how often we would present the updated target position
numConv - how smooth the target motion would be


BrownianPointingSim.m and AnalysisFigs.m

I created a simulation of a person's cursor activity based upon real lag curves (similar to an impulse response) that I got from myself. Earlier, I created simulated data both for a person who adapted to the cursor rotation and for a person who moved their mouse without noticing the rotation. I ran the positions through the simulator to calculate error, lags, etc akin to what exist for the empirical data. This simulated data, once run through AnalysisFigs.m helped to assess whether the empirical data followed the adaptation or no-adaptation behavior.
