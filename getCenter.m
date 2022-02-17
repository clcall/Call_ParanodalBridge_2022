function center = getCenter(an)
%cuprs
if contains(an,'a609_1_r1c1')
    center = [145.85 153.60 31];
elseif contains(an,'a609_1_r1c6')
    center = [110.98 91.33 37];
elseif contains(an,'a609_1_r1c3')
    center = [138.03 130.77 41];
elseif contains(an,'a609_1_r1c4')
    center = [204.80 104.34 27];
elseif contains(an,'a008_0_r1c1')
    center = [141.14 133.36 29];
elseif contains(an,'a008_0_r1c3')
    center = [134.14 138.81 28];
elseif contains(an,'a008_0_r1c8')
    center = [142.57 135.83 28];
elseif contains(an,'a010_0_r1c6')
    center = [189.27 137.98 32];
elseif contains(an,'a010_0_r1c7')
    center = [129.52 173.65 27];
elseif contains(an,'a010_0_r1c8')
    center = [139.07 139.20 19];
elseif contains(an,'a010_0_r2c5')
    center = [204.8 146.6 6];
%ctrls
elseif contains(an,'a411_2_r1c1')
    center = [134.53 137.38 29];
elseif contains(an,'a411_2_r1c4')
    center = [135.44 142.44 16];
elseif contains(an,'a411_2_r2c3')
    center = [137.90 137.51 35];
elseif contains(an,'a614_4_r2c1')
    center = [126.01 130.89 56];
elseif contains(an,'a614_4_r2c2')
    center = [136.51 132.60 39];
elseif contains(an,'a909_4_r1c5')
    center = [135.83 136.73 39];
elseif contains(an,'a909_4_r2c1')
    center = [133.49 134.40 26];
elseif contains(an,'a909_4_r2c2')
    center = [133.88 131.16 24];
elseif contains(an,'a909_4_r2c4')
    center = [138.94 125.97 28];
elseif contains(an,'a909_4_r2c7')
    center = [136.60 137.51 12];
elseif contains(an,'a129_1_r1c1')
    center = [112.6,228.9,22];
elseif contains(an,'a129_3_r1c2')
    center = [121.5 142.5 30];
elseif contains(an,'a129_3_r1c3')
    center = [122.1 120.2 33];
elseif contains(an,'a129_3_r2c4')
    center = [97.8 114.1 33];
elseif contains(an,'a129_6_r1c6')
    center = [122.3 104.9 28];
elseif contains(an,'a129_6_r1c7')
    center = [129.3 124.6 55];
elseif contains(an,'a129_6_r1c8')
    center = [114.7 90.0 32];
%bza stable
elseif contains(an,'a106_5_r1c1')
    center = [368./4.8177 704./4.8177 40];
elseif contains(an,'a106_6_r1c1')
    center = [534./4.8177 530./4.8177 64];
%bsln only
elseif contains(an, 'a105_1_r1c1')
    center = [115 111 25];
elseif contains(an, 'a105_1_r1c2')
    center = [106 134 34];
elseif contains(an, 'a106_5_r1c2')
    center = [195 159 49];
elseif contains(an, 'a106_5_r1c3')
    center = [126 123 33];
elseif contains(an, 'a226_2_r1c1')
    center = [114 130 40];
else
    fprintf(['No match for ' an '. Make sure soma centroid is logged in getCenter.m.\n'])
    return
end
end