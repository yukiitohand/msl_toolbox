function [basename_cache] = mastcam_get_basename_cache(basename_cache_com,data_code,cache_vr)
basename_cache = sprintf('%s_%s_%s',basename_cache_com,data_code,cache_vr);
end