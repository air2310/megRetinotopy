function beta_out = mprfRecomputeBetas(stim, sigma, X, Y, mresp)

beta_out = nan(size(sigma));

fprintf('Recomputing betas:\n')
figure;

for n = 1:1000:length(sigma)
    
    if mod(n-1,10000) == 0
        f_done = n/length(sigma)*100;
        fprintf('%0.2f percent\n',f_done)
    end
    
    end_idx = n+999;
    if end_idx > length(sigma)
        end_idx = length(sigma);
    end
    
    
    
    cur_sigma = sigma(n:end_idx);
    cur_x0 = X(n:end_idx);
    cur_y0 = Y(n:end_idx);
    
    RF = rfGaussian2d(stim.X,stim.Y,cur_sigma,cur_sigma,false, cur_x0,cur_y0);
    
    
    tmp = stim.im_conv' * RF;
    
    beta_out(n:end_idx) = mresp(n:end_idx) ./max(tmp);
    
    plot(nanmean(tmp,2)); hold on;
    
    
end

fprintf('Done\n');

end









































