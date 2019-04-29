function median_itst = get_median_intensities(itst)
itst = reshape(itst, 3, size(itst, 1) / 3, 12);
itst = permute(itst, [1, 3, 2]);

extract_from = [0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0;
                1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0;
                0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0];
            
median_itst = squeeze(median(sum(itst .* extract_from), 2));