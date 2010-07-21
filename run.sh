#! /bin/bash

#TESTING PARAMETER SENSITIVITY
#./profile.o ../images/411-b-14.9-h.tiff ../images/411-b-14.9-v.tiff 1571 687 1547 2290 3037 814 ../images/sensitivity2

./cpdeflec ../images/411-b-14.9-h.tiff ../images/411-b-14.9-v.tiff 1571 687 1547 2290 3037 814 ../images/411-b-14.9
#./profile.o ../images/411-t-14.9-h.tiff ../images/411-t-14.9-v.tiff 1351 679 1325 2304 2824 820 ../images/411-t-14.9
#./profile.o ../images/411-l-14.9-h.tiff ../images/411-l-14.9-v.tiff 1462 668 1436 2281 2923 806 ../images/411-l-14.9
#./profile.o ../images/411-r-14.9-h.tiff ../images/411-r-14.9-v.tiff 1389 669 1366 2288 2856 807 ../images/411-r-14.9
#./profile.o ../images/411-b-14.4-left-h.tiff ../images/411-b-14.4-left-v.tiff 95 522 96 2277 1672 735 ../images/411-b-14.4-left
#./profile.o ../images/411-b-14.4-righ-1.2-h.tiff ../images/411-b-14.4-righ-1.2-v.tiff 2338 493 2350 2165 3961 543 ../images/411-b-14.4-righ-1.2
#./profile.o ../images/411-b-14.4-inf-h.tiff ../images/411-b-14.4-inf-v.tiff 1390 484 1417 2166 2924 596 ../images/411-b-14.4-inf
#./profile.o ../images/411-b-14.4-1-h.tiff ../images/411-b-14.4-1-v.tiff 1386 474 1413 2173 2934 588 ../images/411-b-14.4-1
#./profile.o ../images/411-b-14.4-0.7-h.tiff ../images/411-b-14.4-0.7-v.tiff 1383 469 1410 2177 2940 583 ../images/411-b-14.4-0.7
#./profile.o ../images/411-b-14.0-ora-h.tiff ../images/411-b-14.0-ora-v.tiff 1119 695 1130 2376 2625 820 ../images/411-b-14.0-ora
# Need to swap around pricipal point...
#./profile.o ../images/411-b-14.4-rot-1-h.tiff ../images/411-b-14.4-rot-1-v.tiff 530 1659 734 3257 1983 1579 ../images/411-b-14.4-rot-1

#./profile.o ../images/411-b-14.0-set-h.tiff ../images/411-b-14.0-set-v.tiff 1119 695 1129 2376 2625 820 ../images/411-b-14.0-set

# This panel needs to have its first point shifted.
#./profile.o ../images/037-t-15.4-h.tiff ../images/037-t-15.4-v.tiff 1567 629 1545 2279 3055 769 ../images/037-t-15.4

#./profile.o ../images/015-b-15.4-h.tiff ../images/015-b-15.4-v.tiff 1452 777 1434 2175 2727 885 ../images/015-b-15.4

#./profile.o ../images/018-b-15.0-h.tiff ../images/018-b-15.0-v.tiff 1525 786 1516 2171 2791 882 ../images/018-b-15.0
# Need to move first corner.
#./profile.o ../images/018-r-15.3-h.tiff ../images/018-r-15.3-v.tiff 1504 775 1489 2162 2770 879 ../images/018-r-15.3

#./profile.o ../images/019-b-14.9-h.tiff ../images/019-b-14.9-v.tiff 1367 798 1360 2192 2634 896 ../images/019-b-14.9

#./profile.o ../images/035-r-x-h.tiff ../images/035-r-x-v.tiff 1133 672 1121 2343 2627 818 ../images/035-r-x
#./profile.o ../images/035-r-14.8-h.tiff ../images/035-r-14.8-v.tiff 1335 652 1318 2266 2792 790 ../images/035-r-14.8
#./profile.o ../images/035-l-x-h.tiff ../images/035-l-x-v.tiff 1044 703 1029 2308 2478 848 ../images/035-l-x

#./profile.o ../images/373-b-15.5-h.tiff ../images/373-b-15.5-v.tiff 1421 717 1403 2262 2821 841 ../images/373-b-15.5
#./profile.o ../images/373-b-15.5-nef-h.tiff ../images/373-b-15.5-nef-v.tiff 1421 717 1403 2262 2821 841 ../images/373-b-15.5-nef
#./profile.o ../images/373-t-15.6-h.tiff ../images/373-t-15.6-v.tiff 1412 706 1391 2252 2811 835 ../images/373-t-15.6
#./profile.o ../images/373-t-14.4-b-h.tiff ../images/373-t-14.4-b-v.tiff 1447 740 1421 2277 2836 870 ../images/373-t-14.4-b
#./profile.o ../images/373-t-14.4-a-h.tiff ../images/373-t-14.4-a-v.tiff 1383 717 1354 2259 2771 855 ../images/373-t-14.4-a
#./profile.o ../images/373-r-15.5-h.tiff ../images/373-r-15.5-v.tiff 1413 707 1393 2252 2813 834 ../images/373-r-15.5
#./profile.o ../images/373-l-14.4-h.tiff ../images/373-l-14.4-v.tiff 1590 727 1570 2289 3008 846 ../images/373-l-14.4

#./profile.o ../images/247-t-14.4-h.tiff ../images/247-t-14.4-v.tiff 1425 658 1426 2245 2872 770 ../images/247-t-14.4

#./profile.o ../images/040407-b-14.4-h.tiff ../images/040407-b-14.4-v.tiff 1437 606 1435 2196 2891 725 ../images/040407-b-14.4

#Need to move 1st and 2nd corner.
#./profile.o ../images/210906-b-14.4-h.tiff ../images/210906-b-14.4-v.tiff 1476 597 1472 2193 2958 711 ../images/210906-b-14.4
#Need to move 1st and 2nd corner.
#./profile.o ../images/210906-b-14.4-dirt-h.tiff ../images/210906-b-14.4-dirt-v.tiff 1476 597 1471 2193 2957 711 ../images/210906-b-14.4-dirt
