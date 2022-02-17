function [ctrl, cupr, cupenrich, cupbza, cupsg, ctrlbza, ctrlenrich, ctrlU50, cuprU50, cuprFC, cuprFX] = getFigColors
ctrl = [103 187 247];
cupr = [255 150 38];
cupenrich = [115 247 103];
cupbza = [235 103 247];
ctrlU50 = [103 150 230];
cuprU50 = [28 204 186];
cuprFC = 'b';
cuprFX = [125 125 125];

% cupsg = mean([cupr;cupbza]).*0.7;
cupsg = [0 0 0];

ctrlbza = mean([ctrl;cupbza]);
ctrlenrich = mean([cupenrich;ctrl]);

ctrl = (ctrl./255) - .1;
cupr = (cupr./255);
cupenrich = (cupenrich./255);
ctrlenrich = ctrlenrich./255;
cupbza = (cupbza./255);
cupsg = cupsg./255;
ctrlbza = ctrlbza./255;
ctrlU50 = ctrlU50./255;
cuprU50 = cuprU50./255;
cuprFX = cuprFX./255;