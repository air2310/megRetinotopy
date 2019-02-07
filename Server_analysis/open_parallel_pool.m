function mpool = open_parallel_pool(n_cores, sys)
% Function to start a parallel pool
if ~exist('sys','var') || isempty(sys)
    sys = 'default';
end

switch sys
    case 'Spinoza_grid_matlab2016'
        configCluster;
        ClusterInfo.setEmailAddress('a.edadan@uu.nl');
        c = parcluster;
        mpool = c.parpool(n_cores);

    case 'Spinoza_grid_matlab2018'
        configCluster;
        c = parcluster('p6444');
        c.AdditionalProperties.EmailAddress = 'a.edadan@uu.nl';
        mpool = c.parpool(n_cores);

    case 'NYU_HPC'    
        
    otherwise
        mpool = parpool(n_cores);
        pctRunOnAll warning('off','all');
               
end

end