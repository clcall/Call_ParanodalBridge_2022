function timeline = gettimeline(name)
switch name
    case 'fish01_cellB'
        dec16 = 77;
        dec16b = 80: 1.47 :102.75-6;
        dec17 = 102.75: 1.83 :126.75-3;
        dec18 = 126.75: 1.8 :151.5-8;
        dec19 = 151.5: 1.8 :151.5+(1.8*8);
        timeline = [dec16,dec16b,dec17,dec18,dec19]; 
    case 'fish05_cellG'
        dec16 = 78;
        dec16b = 81: 1.46 :103.5-7;
        dec17 = 103.5: 1.83 :127.25-4;
        dec18 = 127.25: 1.8 :151.5-12;
        timeline = [dec16,dec16b,dec17,dec18];
    case 'fish07_cellE'
        dec16 = 78.5;
        dec16b = 80.5: 1.46 :103-6;
        dec17 = 103: 1.83 :127-4;
        dec18 = 127: 1.8 :151.5-8;
        dec19 = 151.5: 1.8 :151.5+(1.8*9);
        timeline = [dec16,dec16b,dec17,dec18,dec19]; 
    case 'fish08_cellO'
        dec17 = 103: 1.83 :127-4;
        dec18 = 127: 1.8 :151.5-8;
        dec19 = 151.5: 1.8 :151.5+(1.8*9);
        timeline = [dec17,dec18,dec19]; 
    case 'fishA_cell02'
        nov20 = 79: 1.8 :79+(1.8*9);
        nov21 = 101: 1.16 :101+(1.16*1);
        nov21b = 105.5: 1.98 :105.5+(1.98*8);
        nov22 = 126: 1.09 :126+(1.09*1);
        nov22b = 130.75: 1.71 :130.75+(1.71*6);
        nov23 = 149: 2.05 :149+(2.05*11);
        nov25 = 197;
        timeline = [nov20,nov21,nov21b,nov22,nov22b,nov23,nov25];
    case 'fishC_cell01'
        nov20 = 79.5: 1.8 :79.5+(1.8*9);
        nov21 = 101.5: 1.16 :101.5+(1.16*1);
        nov21b = 106: 1.98 :106+(1.98*8);
        nov22 = 127.5;
        nov22b = 129.5: 1.71 :129.5+(1.71*6);
        nov23 = 150: 2.05 :150+(2.05*10);
        nov25 = 199;
        timeline = [nov20,nov21,nov21b,nov22,nov22b,nov23,nov25];
    case 'fishD_cell01'
        nov20 = 80: 1.8 :80+(1.8*9);
        nov21 = 102: 1.16 :102+(1.16*1);
        nov21b = 106.5: 1.98 :106.5+(1.98*8);
        nov22 = 126.25;
        nov22b = 130: 1.71 :130+(1.71*6);
        nov23 = 150.5: 2.05 :150.5+(2.05*8);
        timeline = [nov20,nov21,nov21b,nov22,nov22b,nov23];
end
