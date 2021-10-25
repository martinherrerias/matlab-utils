function [I]=geotiffread(filename)
% GEOTIFF_READ: read geotiff using imread and assign map info from infinfo.
%
% output:
% I.z, image data
% I.x, x coordinate in map
% I.y, y coordinate in map
% I.info, misc. info
%
% imshow(I.z, 'xdata', I.x, 'ydata', I.y);
% shows image with map coordinate

% Version by Yushin Ahn, ahn.74@osu.edu
% Glacier Dynamics Laboratory, 
% Byrd Polar Resear Center, Ohio State University 
% Referenced enviread.m (Ian Howat)
% 
% Copyright (c) 2010, Yushin Ahn
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the Ohio State University nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

    Tinfo = imfinfo(filename);
    info.samples = Tinfo(1).Width;
    info.lines = Tinfo(1).Height;
    info.imsize = Tinfo(1).Offset;
    info.bands = Tinfo(1).SamplesPerPixel;

    sub = [1, info.samples, 1, info.lines];

    info.map_info.dx = Tinfo(1).ModelPixelScaleTag(1);
    info.map_info.dy = Tinfo(1).ModelPixelScaleTag(2);
    info.map_info.mapx = Tinfo(1).ModelTiepointTag(4);
    info.map_info.mapy = Tinfo(1).ModelTiepointTag(5);

    xm = info.map_info.mapx;
    ym = info.map_info.mapy;
    xx = xm + ((0:info.samples-1).*info.map_info.dx);
    yy = ym - ((0:info.lines  -1).*info.map_info.dy);

    I.x = xx(sub(1):sub(2));
    I.y = yy(sub(3):sub(4));
    I.z = imread(filename);
    I.info = info;
end





