function andy_mk_video(fname, pp);
   dd = eidors_readdata(fname);
   pp = set_parameters(pp,dd);
   dd = preproc_data(dd,pp);
   [imdl,ROIs,pp] = mk_imdl(pp); pp.ROIs = ROIs;
   img= inv_solve(imdl, dd(:,pp.hi),dd);
   switch pp.show
      case 'static'; pp.fname = fname;
          static_show(img,pp);
      case 'dynamic'; pp.fname = fname;
%         for FS = [105:5:250]; pp.FS = FS;
          for FS = [10:5:5000]; pp.FS = FS;
          clf;
          dynamic_show(img,pp);
          drawnow;
              if isfield(pp,'video')
                  print('-dpng', '-r75',sprintf( ...
                    '%s-im%04d.png',pp.video,FS))
              end
          end
      otherwise;     error('huh?');
   end


function  pp = set_parameters(pp,dd);
   if ~isfield(pp,'Nel'); 
       switch size(dd,1);
          case  208; pp.Nel = 16;
          case 1024; pp.Nel = 32;
          otherwise; error 'huh?';
       end
   end
   if ~isfield(pp,'FR'); 
       pp.FR = 50;
   end
   if ~isfield(pp,'LPF');
       pp.LPF = 1.2;
   end

function dd = preproc_data(di,pp)
    di = real(freq_filt(di,@(f) f<pp.LPF,pp.FR,2));
    
    [~, msel] = mk_stim_patterns(pp.Nel,1,[0,1],[0,1],{'no_meas_current'},1);
    dd = zeros(pp.Nel^2,size(di,2));
    dd(msel, :) = di; 

function pp = get_boundary(fmdl,pp);
    fmdl2 = mdl2d_from3d(fmdl);
    b = fmdl2.boundary; bs = size(b);
    pp.boundaryx = reshape(fmdl2.nodes(b,1),bs)';
    pp.boundaryy = reshape(fmdl2.nodes(b,2),bs)';

function [img,vh,vi] = set_elem_background(fmdl)
    img = mk_image(fmdl,1);
    vh = fwd_solve(img);
%   img.elem_data(vertcat(fmdl.mat_idx{2:3}))= 2;
    m_frac = elem_select(fmdl, @(x,y,z) ...
           (x.^2)/0.6 + (y-0.25).^2/0.55 < 1);
    img.elem_data = img.elem_data + m_frac*2;
    vi = fwd_solve(img);
 

function [imdl,ROIs,pp] = mk_imdl(pp)
%   fmdl.electrode = fmdl.electrode(fliplr([9:16,1:8]));
    switch 1;
       case 0;
          fmdl= mk_library_model('adult_male_32el_lungs');
          fmdl.electrode(2:2:32)=[];
          evol = get_elem_volume(fmdl);
          pp.lmag   = sum(evol(fmdl.mat_idx{3})) / ...
                      sum(evol(fmdl.mat_idx{2}));
       case 1;
          fmdl = mk_library_model({'adult_male','boundary'},...
                  [16 1.015 0.5],[0.05],0.08);
          pp.lmag = 1.2;
       case 2;
          fmdl = mk_library_model({'adult_male','boundary'},...
                  [16 1.000 0.5],[0.05],0.08);
          pp.lmag = 1.5;
    end
    pp = get_boundary(fmdl,pp);

    fmdl = mdl_normalize(fmdl,1);
    [fmdl.stimulation,fmdl.meas_select] = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current_next1'},1);
    
    [img,vh,vi] = set_elem_background(fmdl);
    opt.imgsz = [100 100];
    opt.noise_figure = 0.8;
    opt.square_pixels = 1;
    imdl=mk_GREIT_model(img, 0.25, [], opt);
    rmdl = imdl.rec_model;
    xy = rmdl.coarse2fine'*interp_mesh(rmdl);
    ROI = getfield(inv_solve(imdl,vh,vh.meas*2),'elem_data')<-10;
    %ROIl=ROI.*(xy(:,1)<0);
    %ROIr=ROI.*(xy(:,1)>0);
    ROIl=(xy(:,1)<0);
    ROIr=(xy(:,1)>0);
    ROIs = [ROIl(:), ROIr(:), ROI(:)];
    ROI = getfield(inv_solve(imdl,vh,vi),'elem_data')>500;
    ROIs = [ROIs, ROI(:)];

function cmap = tidal_colourmap(pp)
%   sd  = -scale_data*100;
    sd  = -linspace(-100,100,255);
    LIM = 30;
    BG  = sd< LIM;
    p20 = (sd>LIM) .* (sd<= 100) .* (sd -LIM) / (100-LIM);
    P20 = (p20>0);
    p10 = (sd> 10) .* (sd<= LIM) .* (sd - 10) / (LIM-10);
    n20 = (sd<-20) .* (sd>=-100) .* (sd + 20) /-(100-20);
    n10 = 0; n20=0;
    red = 255*p20 + 170*n20;
    grn = 255*p20 ;
    blu = 192*P20 + 192*p10 + 170*n20;
    bgval = min(255,pp.bgval);
    red = red + bgval*BG.*(1-p10);
    grn = grn + bgval*BG.*(1-p10);
    blu = blu + bgval*BG.*(1-p10);

    red = red/255;
    grn = grn/255;
    blu = blu/255;
    cmap = [red(:),grn(:),blu(:)];

function dynamic_show(img,pp);
    tt = (1:size(img.elem_data,2))/50;
    ROIs = pp.ROIs(:,1:2)/sum(pp.ROIs(:,1:2),[1,2]);
    sig = -img.elem_data.' * ROIs;
    sig(:,2) = sig(:,2) * pp.lmag;
    vsig = sig + [pp.volofs,0];
    sig(:,2) = sig(:,2) +  pp.volsigofs;
    flow = [0,0;diff(sig)] + [pp.flowofs,0];
    grey = 0.7*[1,1,1]; blu=[2,2,8]/10; red=[8,2,2]/10;
    colorrdr=[(4*grey+0*blu);(4*grey+0*red);
              (3*grey+1*blu);(3*grey+1*red);
              (2*grey+2*blu);(2*grey+2*red);
              (1*grey+3*blu);(1*grey+3*red);
              (0*grey+4*blu);(0*grey+4*red)]/4;
    lw = 'LineWidth';

%   subplot(3,2,[2,4]);
%   subplot(2,2,2);
    axes('Position',[0.52,0.44,0.4,0.5]);
    mp100 = max(5,pp.FS- 50);
    mp200 = max(4,pp.FS-150);
    mp300 = max(3,pp.FS-250);
    mp400 = max(2,pp.FS-350);
    idx4 = 1:mp400;
    idx3 = mp400:mp300;
    idx2 = mp300:mp200;
    idx1 = mp200:mp100;
    idx0 = mp100:pp.FS;
    set(gca,'Colororder',colorrdr,'NextPlot','ReplaceChildren')
   
    hh= plot(-sig(idx4,:),-flow(idx4,:),  ...
             -sig(idx3,:),-flow(idx3,:),  ...
             -sig(idx2,:),-flow(idx2,:),  ...
             -sig(idx1,:),-flow(idx1,:),  ...
             -sig(idx0,:),-flow(idx0,:), lw, 2 );
    set(hh(1:2),lw,1);
    set(hh(9:10),lw,3);
    minsig = min(-sig(:));  maxsig = max(-sig(:));
    minflow= min(-flow(:)); maxflow= max(-flow(:));
    axis([minsig, maxsig, minflow, maxflow]);
    lw = {'LineWidth',1,'Color',[0,0,0],'LineStyle','--'};
    line([minsig;maxsig],-[1;1]*[0,pp.flowofs],lw{:});
    line([0,0],[minflow;maxflow],lw{:});
    axis off

    %subplot(2,1,2);
    axes('position',[0.1000, 0.1000, 0.8, 0.3 ]);
    idx1 = 1:pp.FS;

    idx2 = pp.FS+1:length(tt);
    colorrdr=[grey;grey;blu;red];
    set(gca,'Colororder',colorrdr,'NextPlot','ReplaceChildren')
    plot(tt(idx2),vsig(idx2,:), ...
         tt(idx1),vsig(idx1,:),'LineWidth',2)
    ylim(pp.yl); xlim([0,100])
    box off; set(gca,'YTickLabel',[]);
    vert_lines(tt(pp.FS),pp.yl);


    img.calc_colours.backgnd = [1,1,1];
    img.calc_colours.ref_level = 0;
    img.calc_colours.clim = pp.clim;
    img.calc_colours.greylev = 0.05;
    img.calc_colours.cmap_type = tidal_colourmap(pp);
 
    img.get_img_data.frame_select = pp.FS;
    axes('Position',[0.1,0.28,0.4,0.8]);

    img.elem_data = img.elem_data .* (pp.ROIs(:,3)*0.3+0.7);
    img.fwd_model = rmfield(img.fwd_model,'electrode');
    hh=show_fem(img); set(hh,'EdgeColor','none');
    if pp.bgval<=255
    line(pp.boundaryx, pp.boundaryy, ...
             'Color',60/255*[1,1,1],'LineWidth',4);
    end
    xlim(mean(xlim)+1.02*(xlim-mean(xlim)))
    ylim(mean(ylim)+1.02*(ylim-mean(ylim)))
    axis off

    kids = get(gcf,'Children');
    set(gcf,'Children',kids([3,2,1]));

function static_show(img,pp);
    subplot(2,1,1);
    tt = (1:size(img.elem_data,2))/50;
    sig = -(pp.ROIs/sum(pp.ROIs,[1,2]))'*img.elem_data;
    plot(tt,sig,'LineWidth',2)
    vert_lines(tt(pp.FS),pp.yl);
    box off; title(pp.fname); ylim(pp.yl);

    img.show_slices.img_cols=min(12,length(pp.FS));
    img.calc_colours.backgnd = [1,1,1];
    img.calc_colours.ref_level = 0;
    img.calc_colours.greylev = 0.05;
    img.calc_colours.cmap_type = tidal_colourmap;
 
    img.get_img_data.frame_select = pp.FS;

    subplot(2,1,2); show_slices(img);

function hh=vert_lines(loc,yl)
    hh= line([1;1]*loc,yl'*ones(size(loc)), ...
           'Color',[0,0,0]);

