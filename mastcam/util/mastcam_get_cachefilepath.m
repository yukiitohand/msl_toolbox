function [cache_imgpath,cache_hdrpath] = mastcam_get_cachefilepath(basename_cache,dirpath_cache)
cache_imgpath = joinPath(dirpath_cache,[basename_cache '.img']);
cache_hdrpath = joinPath(dirpath_cache,[basename_cache '.hdr']);
end