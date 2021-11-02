function cMap=vef_colormap(clim,cfun,varargin)

colormap(cfun(varargin{:}));
caxis(clim);
