% QuestTest.m
%
% Created by hacking a copy of QuestDemo.m, in order to
% do some simulations for Beau Watson.
%
% 128 repetitions of each condition
% Quest assumes LogWeibull psychometric function with parameters:
% beta = {1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5} (in separate conditions)
% gamma = 0.03
% delta = .01
% sample range = {-30,30} dB
% sample step size 0.5 dB
% Prior is always Gaussian, sd=12 dB, centered on guess.
% Initial guess drawn from Normal distribution with mean 0 and sd=6 dB.
% Quest uses mode of posterior density.
% Simulated observer has LogWeibull psychometric function with parameters:
% {beta=3.5,gamma=.03,delta=.01}
%
%
% By commenting and uncommenting five lines below, you can use
% this file to implement three QUEST-related procedures for measuring
% threshold.
%
% QuestMode: In the original algorithm of Watson & Pelli (1983)
% each trial and the final estimate are at the MODE of the posterior pdf.
%
% QuestMean: In the improved algorithm of King-Smith et al. (1994).
% each trial and the final estimate are at the MEAN of the posterior pdf.
%
% QuestQuantile & QuestMean: In the ideal algorithm of Pelli (1987)
% each trial is at the best QUANTILE, and the final estimate is at 
% the MEAN of the posterior pdf.
%
% King-Smith, P. E., Grigsby, S. S., Vingrys, A. J., Benes, S. C., and Supowit, A.
% (1994) Efficient and unbiased modifications of the QUEST threshold method: theory, 
% simulations, experimental evaluation and practical implementation. 
% Vision Res, 34 (7), 885-912.
%
% Pelli, D. G. (1987) The ideal psychometric procedure. Investigative Ophthalmology 
% & Visual Science, 28 (Suppl), 366.
%
% Watson, A. B. and Pelli, D. G. (1983) QUEST: a Bayesian adaptive psychometric 
% method. Percept Psychophys, 33 (2), 113-20.
%
% Copyright (c) 1996-1997 Denis G. Pelli
%
% 6/24/97  dgp	wrote it, based on QuestDemo.m

fprintf('beta	mean	sd		reps\n')
reps=64;
for beta = [1.5:.5:5.5]
	estimate=1:reps;
	for i=1:reps
		gamma = 0.03;
		delta = .01;
		tGuess=randn*6/20;
		tGuessSd=12/20;
		range = 3;
		grain=0.5/20;
		pThreshold=0.6370;
		tActual=0;
		trialsDesired=64;
		wrongRight=char('wrong','right');
		q=QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma,grain,range);
		qObserver=QuestCreate(tActual,tGuessSd,pThreshold,3.5,delta,gamma,grain,range);
		%fprintf('Quest estimate of thresh after initialization is %g\n',QuestMean(q));
		%fprintf('Initial value was %g\n',tGuess);
		%fprintf('Quest estimate of sd after initialization is %g\n',QuestSd(q));
		%fprintf('Initial value was %g\n',tGuessSd);
		
		for k=1:trialsDesired
			tTest=QuestMode(q);	% 0.03 s
		 	response=QuestSimulate(qObserver,tTest,tActual); % 0.10 s
			q=QuestUpdate(q,tTest,response); % 0.10 s
		end
		t=QuestMode(q);
		%fprintf('%4.1f\t%+6.3f\n',beta,t);
		estimate(i)=t;
	end
	fprintf('%4.1f\t%+6.3f\t%5.3f\t%d\n',beta,mean(estimate),std(estimate),length(estimate));
end
	
		
		
