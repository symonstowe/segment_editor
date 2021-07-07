    fn = '2432623_xueliu_01.eit'; hi=183;
    pp.FR= 20;
    yl=3*[-4,47];
    FS = hi+[1:36]*10;
    pp.clim =   250;
    pp.video = 'Jon-bg1';
    pp.flowofs = 0.8;
    pp.volofs = 40;
    pp.volsigofs = 0;
    pp.bgval = 128;%pp.video = 'Jon-bg3';

set(gcf,'PaperPosition',[1,1,13.2,8]);


pp.yl = yl;
pp.FS = FS;
pp.hi = hi;
pp.LPF= 0.8;
pp.show = 'dynamic';
andy_mk_video(fn, pp)
