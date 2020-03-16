function map = bw
% create red-white-blue colorbar
%   written by Tomoyasu Horikawa 20160310


%% color configuration : from blue to to white then to red

%% color interpolation
T = [8,54,105
    20,77,138
    35,104,173
    51,126,184
    71,149,196
    104,171,207
    145,196,221
    175,212,230
    208,228,239
    226,237,243
    245,245,245
    249,234,225
    252,220,200
    248,193,166
    244,165,130
    230,132,105
    214,97,78
    196,60,60
    177,23,42
    138,11,36
    103,0,31]./255;
x = round(linspace(0,255,size(T,1)));
map = interp1(x/255,T,linspace(0,1,1238));
map=map(end/2:-1:1,:);
% plot(map)

%% debug
if 0
r=normpdf([-2:4/64:2]-0.5)';
g=normpdf([-2:4/64:2])';
b=normpdf([-2:4/64:2]+0.5)';
map=rwb;
map=map./max(map(:));
imagesc(randn(10))
colormap mem
subplot(211)
plot([r,g,b])
subplot(212)
plot(rwb)

r=normpdf([-1:2/64:1]-0.5)';
g=normpdf([-1:2/64:1])';
b=normpdf([-1:2/64:1]+0.5)';
map=[r,g,b];
map=map./max(map(:));
end

%% debug
if 0
red_top     = [];
white_middle= [1 1 1];
blue_bottom = [0.1 0.3 0.8];
T = [0.8 0 0.1
    1 1 1
    0.1 0.3 0.8];
x = [0
    127
    255];
T = [0,   0,   0          %// dark
    101, 67,  33         %// brown
    255, 105, 180        %// pink
    255, 255, 255        %// white
    255, 255, 255]./255; %// white again  -> note that this means values between 161 and 255 will be indistinguishable
x = [0
    50
    120
    160
    255];
T = [8,54,105
    20,77,138
    35,104,173
    51,126,184
    71,149,196
    104,171,207
    145,196,221
    175,212,230
    208,228,239
    226,237,243
    245,245,245
    249,234,225
    252,220,200
    248,193,166
    244,165,130
    230,132,105
    214,97,78
    196,60,60
    177,23,42
    138,11,36
    103,0,31]./255;
x = round(linspace(0,255,size(T,1)));
map = interp1(x/255,T,linspace(0,1,255));
I = linspace(0,1,255);
imagesc(I(ones(1,10),:)')
colormap(map)

end