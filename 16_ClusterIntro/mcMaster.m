%This script is meant to be run on the Monte Carlo experiments on the PSU
%cluster. 

numSim = 1000;

for M = [250]
    for beta = [20 50 80]
        for delta = [0 80]
            tag = sprintf('b%dd%dM%d', beta, delta, M);
            fprintf('Starting Monte Carlo: %s\n', tag);
            
            %We just want to show how to submit a job.
            %runMonte(beta/100, delta/100, M, numSim, tag);
            pause(10)
            
            %Clean up data files (can re-create them using seed):
            %Should no longer be necessary.
            %unix(sprintf('rm *Data*%s*.mat', tag));
        end
    end
end

