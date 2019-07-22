function mySaveTiffStack16bit(img, fn)


% pref the TIFF object
t = Tiff(fn, 'w');
% http://www.mathworks.com/help/matlab/import_export/exporting-to-images.html
tags.ImageLength   = size(img, 1);
tags.ImageWidth    = size(img, 2);
tags.Photometric   = Tiff.Photometric.MinIsBlack;
tags.BitsPerSample = 16;
tags.SampleFormat  = Tiff.SampleFormat.Int;
tags.RowsPerStrip  = 16;
tags.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tags.SamplesPerPixel = 1;
% unspec               = Tiff.ExtraSamples.Unspecified;
% tags.ExtraSamples    = [ unspec, unspec, unspec, unspec, unspec];
t.setTag(tags);
%%
disp(t);
%%

for f = 1 : size(img, 3)
    t.setTag(tags)
    t.write(img(:, :, f), 'WriteMode', 'append');
    t.writeDirectory();  
end
t.close();
%%
% wdata = imread('output.tif');
% flag  = isequal(img,wdata);
% disp(flag);