
initials = 'KLS';
sinPercent = [0 0.25 0.5 0.75 1];
sinPercentShuffle = sinPercent(randperm(numel(sinPercent)));
tfreq = 0.25;
stepSize = 5;
updateRate = 5;
numConv = 8;
for i = sinPercent
BrownianTrialPerturb(initials,i,tfreq,stepSize,updateRate,numConv);
end