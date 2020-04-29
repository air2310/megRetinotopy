function stim = mprfLoadRMStimulus_02(vw)

if notDefined('vw')
    error('Need volume view');
end


rmparams = viewGet(vw,'rmparams');

stim.im = rmparams.stim.images_unconvolved;
stim.im_conv = rmparams.analysis.allstimimages';

stim.window = rmparams.stim.stimwindow;
stim.X = rmparams.analysis.X;
stim.Y = rmparams.analysis.Y;

[stim.full_x, stim.full_y] = meshgrid(unique(stim.X), unique(stim.Y));
[~, stim.full_im] = rmStimulusMatrix(rmparams,[],[],false,false);





return


% 
% stim.full_im = zeros([size(stim.full_x) size(stim.im,2)]);
% 
% tmp1 = zeros(numel(stim.full_x),1);
% for n = 1:size(stim.full_im,3)
%    tmp2 = tmp1;
%    tmp2(stim.window) = stim.im(:,n);
%    stim.full_im(:,:,n) = reshape(tmp2,size(stim.full_x));
% end
% 

% OLD CODE FROM rmStimulusMatrix
% modelCoords = [stim.Y(:) stim.X(:)];
% imageCoords = [Y(:) X(:)];
% 
% images = zeros([numel(Y) size(stim.im,2)]);
% images2 = zeros([size(Y) size(stim.im,2)]);
% [~,Ia,Ib] = intersectCols(imageCoords',modelCoords');
% 
% for n = 1:size(stim.im,2)
%     
%     images(Ia,n) = stim.im(Ib,n);
%     images2(:,:,n) = reshape(images(:,n),size(Y));
%     
%     
% end
% 
% 


end