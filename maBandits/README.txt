        maBandits: Matlab files for multi-armed bandits
    
      authors: Olivier Cappé, Aurélien Garivier, Emilie Kaufmann
     ............................................................


These are written in object oriented Matlab and require the Stats
toolbox (at least for policyBayesUCB, policyCPUCB, policyThompson,
gamePoisson and plotResults).


Files included in the current release:

    demo.m                Script to run several policies several times and plot
                            the result; this is the one you should try first
    demo.png              Expected graphical output of demo (fourth plot)
    
    experiment.m          Runs a given policy several times
    plotResults.m         Plots the results (saved in mat files) for several
                            algorithms in a given scenario
    
    Game.m                Generic game class
    gameBernoulli.m       Bernoulli bandit environment
    gameExp.m             (possibly bounded) exponential bandit environment
    gamePoisson.m         (possibly bounded) Poisson bandit environment
    
    
    Policy.m              Generic policy class
    policyBayesUCB.m      Policies, names should be explicit
    policyCPUCB.m      
    policyDMED.m
    policyGittins.m
    policyKLUCB.m
    policyKLUCBexp.m
    policyKLUCBpoisson.m
    policyKLempUCB.m
    policyMOSS.m
    policyThompson.m
    policyUCB.m
    policyUCBtuned.m
    policyUCBV.m
    
    data/                 Sub-directory that contains a mat file needed by
                             Gittins (see below)
    results/              Empty directory for the results


# HOWTO: 

Run demo: for each of the four scenarios, wait for a few minutes and obtain
a bunch of .mat files (one for each policy) in the results/ subdirectory as
well a plot that should look like the one in the file demo.png included in the
main directory (which corresponds to the fourth scenario: bounded Poisson
rewards). See experiment.m for the format of the saved .mat files and
plotResults.m for the explanation of the plots.


# NOTES:

About the rewards: Some policies (policyUCB, policyKLempUCB, policyKLUCB...)
require the rewards to be bounded in [0,1] (policyGittins even requires the
rewards to be binary). This is why, in gameExp.m and gamePoisson.m, the
rewards are normalized when the parameter bound is set (that is, when the arms
are bounded Exponential or bounded Poisson distributed). But the policies
policyKLUCBpoisson and policyKLUCBexp, for obvious reasons, do not work with
[0,1]-valued rewards: the first requires integers, the second requires
positive reals. Thus, they accept a third parameter 'bound' to be set when
called with bounded distributions, so that they can 'de-normalize' the
rewards.
    
About Gittins: Computing Gittins indices require much memory and time as we
did no use any approximation of the indices. To alleviate this issue
policyGittins uses a pre-computed table of indices gathered during preliminar
experiments. This file is saved in Matlab's 5.0 matfile format and should be
compatible with any current version of matlab. With the provided data file,
trying to use policyGittins with an horizon larger than 500 or with very
different parameters will be much slower. To update the table with your own
computations, consider using the method saveGI in policyGittins.m.


    --
    
    $Id: README.txt,v 1.18 2012-07-06 16:29:42 cappe Exp $
