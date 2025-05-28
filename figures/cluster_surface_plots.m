FAPS = spm_read_vols(spm_vol('FAPS.nii'));

%% INS (sagittal) R
figure;surf(flipud(squeeze(FAPS(60,:,:))')); % 60th slice in x dimension -> sagittal slice

figure;surf(flipud(squeeze(FAPS(60,50:70,10:30))')) % bounding box in y- and z-dimensions
colormap(canlabCmap)
caxis([-.0002 .0002])
colorbar;
set(gca,'visible','off')

%% INS (sagittal) L
figure;surf(flipud(squeeze(FAPS(15,:,:))'))
caxis([-.0002 .0002])

figure;surf(flipud(squeeze(FAPS(15,50:70,10:30))'))
colormap(canlabCmap)
caxis([-.0002 .0002])
colorbar;
set(gca,'visible','off')
%% AMY/HTH (coronal)
figure;surf(flipud(squeeze(FAPS(:,50,:))'))

figure;surf(flipud(squeeze(FAPS(15:65,50,10:30))'))
colormap(canlabCmap)
caxis([-.0003 .0003])
colorbar;
set(gca,'visible','off')


%% PCC (sagittal)

figure;surf(flipud(squeeze(FAPS(39,30:50,30:50))'))
colormap(canlabCmap)
caxis([-.0002 .0002])
colorbar;
set(gca,'visible','off')


%% HTH (sagittal)

figure;surf(flipud(squeeze(FAPS(39,40:60,15:30))'))
colormap(canlabCmap)
caxis([-.0002 .0002])
colorbar;
set(gca,'visible','off')