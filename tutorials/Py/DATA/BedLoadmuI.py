from numpy import *
yMuI = array([ 0.      ,  0.003125,  0.00625 ,  0.009375,  0.0125  ,  0.015625,
        0.01875 ,  0.021875,  0.025   ,  0.028125,  0.03125 ,  0.034375,
        0.0375  ,  0.040625,  0.04375 ,  0.046875,  0.05    ,  0.053125,
        0.05625 ,  0.059375,  0.0625  ,  0.065625,  0.06875 ,  0.071875,
        0.075   ,  0.078125,  0.08125 ,  0.084375,  0.0875  ,  0.090625,
        0.09375 ,  0.096875,  0.1     ,  0.103125,  0.10625 ,  0.109375,
        0.1125  ,  0.115625,  0.11875 ,  0.121875,  0.125   ,  0.128125,
        0.13125 ,  0.134375,  0.1375  ,  0.140625,  0.14375 ,  0.146875,
        0.15    ,  0.153125,  0.15625 ,  0.159375,  0.1625  ,  0.165625,
        0.16875 ,  0.171875,  0.175   ,  0.178125,  0.18125 ,  0.184375,
        0.1875  ,  0.190625,  0.19375 ,  0.196875,  0.2     ,  0.203125,
        0.20625 ,  0.209375,  0.2125  ,  0.215625,  0.21875 ,  0.221875,
        0.225   ,  0.228125,  0.23125 ,  0.234375,  0.2375  ,  0.240625,
        0.24375 ,  0.246875,  0.25    ,  0.253125,  0.25625 ,  0.259375,
        0.2625  ,  0.265625,  0.26875 ,  0.271875,  0.275   ,  0.278125,
        0.28125 ,  0.284375,  0.2875  ,  0.290625,  0.29375 ,  0.296875,
        0.3     ,  0.303125,  0.30625 ,  0.309375,  0.3125  ,  0.315625,
        0.31875 ,  0.321875,  0.325   ,  0.328125,  0.33125 ,  0.334375,
        0.3375  ,  0.340625,  0.34375 ,  0.346875,  0.35    ,  0.353125,
        0.35625 ,  0.359375,  0.3625  ,  0.365625,  0.36875 ,  0.371875,
        0.375   ,  0.378125,  0.38125 ,  0.384375,  0.3875  ,  0.390625,
        0.39375 ,  0.396875,  0.4     ,  0.403125,  0.40625 ,  0.409375,
        0.4125  ,  0.415625,  0.41875 ,  0.421875,  0.425   ,  0.428125,
        0.43125 ,  0.434375,  0.4375  ,  0.440625,  0.44375 ,  0.446875,
        0.45    ,  0.453125,  0.45625 ,  0.459375,  0.4625  ,  0.465625,
        0.46875 ,  0.471875,  0.475   ,  0.478125,  0.48125 ,  0.484375,
        0.4875  ,  0.490625,  0.49375 ,  0.496875,  0.5     ,  0.503125,
        0.50625 ,  0.509375,  0.5125  ,  0.515625,  0.51875 ,  0.521875,
        0.525   ,  0.528125,  0.53125 ,  0.534375,  0.5375  ,  0.540625,
        0.54375 ,  0.546875,  0.55    ,  0.553125,  0.55625 ,  0.559375,
        0.5625  ,  0.565625,  0.56875 ,  0.571875,  0.575   ,  0.578125,
        0.58125 ,  0.584375,  0.5875  ,  0.590625,  0.59375 ,  0.596875,
        0.6     ,  0.603125,  0.60625 ,  0.609375,  0.6125  ,  0.615625,
        0.61875 ,  0.621875,  0.625   ,  0.628125,  0.63125 ,  0.634375,
        0.6375  ,  0.640625,  0.64375 ,  0.646875,  0.65    ,  0.653125,
        0.65625 ,  0.659375,  0.6625  ,  0.665625,  0.66875 ,  0.671875,
        0.675   ,  0.678125,  0.68125 ,  0.684375,  0.6875  ,  0.690625,
        0.69375 ,  0.696875,  0.7     ,  0.703125,  0.70625 ,  0.709375,
        0.7125  ,  0.715625,  0.71875 ,  0.721875,  0.725   ,  0.728125,
        0.73125 ,  0.734375,  0.7375  ,  0.740625,  0.74375 ,  0.746875,
        0.75    ,  0.753125,  0.75625 ,  0.759375,  0.7625  ,  0.765625,
        0.76875 ,  0.771875,  0.775   ,  0.778125,  0.78125 ,  0.784375,
        0.7875  ,  0.790625,  0.79375 ,  0.796875,  0.8     ,  0.803125,
        0.80625 ,  0.809375,  0.8125  ,  0.815625,  0.81875 ,  0.821875,
        0.825   ,  0.828125,  0.83125 ,  0.834375,  0.8375  ,  0.840625,
        0.84375 ,  0.846875,  0.85    ,  0.853125,  0.85625 ,  0.859375,
        0.8625  ,  0.865625,  0.86875 ,  0.871875,  0.875   ,  0.878125,
        0.88125 ,  0.884375,  0.8875  ,  0.890625,  0.89375 ,  0.896875,
        0.9     ,  0.903125,  0.90625 ,  0.909375,  0.9125  ,  0.915625,
        0.91875 ,  0.921875,  0.925   ,  0.928125,  0.93125 ,  0.934375,
        0.9375  ,  0.940625,  0.94375 ,  0.946875,  0.95    ,  0.953125,
        0.95625 ,  0.959375,  0.9625  ,  0.965625,  0.96875 ,  0.971875,
        0.975   ,  0.978125,  0.98125 ,  0.984375,  0.9875  ,  0.990625,
        0.99375 ,  0.996875,  1.      ])
UMuI = array([  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   1.00000000e-08,   5.00000000e-08,
         9.00000000e-08,   1.50000000e-07,   2.30000000e-07,
         3.10000000e-07,   4.20000000e-07,   5.30000000e-07,
         6.60000000e-07,   8.00000000e-07,   9.60000000e-07,
         1.13000000e-06,   1.32000000e-06,   1.52000000e-06,
         1.73000000e-06,   1.97000000e-06,   2.22000000e-06,
         2.48000000e-06,   2.77000000e-06,   3.07000000e-06,
         3.39000000e-06,   3.72000000e-06,   4.08000000e-06,
         4.46000000e-06,   4.86000000e-06,   5.28000000e-06,
         5.72000000e-06,   6.19000000e-06,   6.69000000e-06,
         7.21000000e-06,   7.76000000e-06,   8.35000000e-06,
         8.97000000e-06,   9.63000000e-06,   1.03300000e-05,
         1.10700000e-05,   1.18700000e-05,   1.27400000e-05,
         1.36700000e-05,   1.46900000e-05,   1.58100000e-05,
         1.70500000e-05,   1.84600000e-05,   2.00700000e-05,
         2.19600000e-05,   2.42500000e-05,   2.70700000e-05,
         3.06600000e-05,   3.52600000e-05,   4.11100000e-05,
         4.83900000e-05,   5.72600000e-05,   6.77900000e-05,
         8.00500000e-05,   9.40900000e-05,   1.09920000e-04,
         1.27580000e-04,   1.47080000e-04,   1.68420000e-04,
         1.91620000e-04,   2.16690000e-04,   2.81290000e-04,
         3.45060000e-04,   4.07990000e-04,   4.70100000e-04,
         5.31380000e-04,   5.91830000e-04,   6.51440000e-04,
         7.10230000e-04,   7.68190000e-04,   8.25320000e-04,
         8.81620000e-04,   9.37090000e-04,   9.91740000e-04,
         1.04555000e-03,   1.09853000e-03,   1.15068000e-03,
         1.20200000e-03,   1.25250000e-03,   1.30216000e-03,
         1.35100000e-03,   1.39900000e-03,   1.44618000e-03,
         1.49252000e-03,   1.53804000e-03,   1.58272000e-03,
         1.62658000e-03,   1.66961000e-03,   1.71181000e-03,
         1.75317000e-03,   1.79371000e-03,   1.83342000e-03,
         1.87230000e-03,   1.91035000e-03,   1.94757000e-03,
         1.98396000e-03,   2.01952000e-03,   2.05426000e-03,
         2.08816000e-03,   2.12123000e-03,   2.15347000e-03,
         2.18489000e-03,   2.21547000e-03,   2.24522000e-03,
         2.27415000e-03,   2.30224000e-03,   2.32951000e-03,
         2.35595000e-03,   2.38155000e-03,   2.40633000e-03,
         2.43028000e-03,   2.45339000e-03,   2.47568000e-03,
         2.49714000e-03,   2.51777000e-03,   2.53757000e-03,
         2.55654000e-03,   2.57468000e-03,   2.59199000e-03,
         2.60847000e-03,   2.62412000e-03,   2.63895000e-03,
         2.65294000e-03,   2.66610000e-03,   2.67844000e-03,
         2.68994000e-03,   2.70061000e-03,   2.71046000e-03,
         2.71947000e-03,   2.72766000e-03,   2.73502000e-03,
         2.74154000e-03,   2.74724000e-03,   2.75211000e-03,
         2.75614000e-03,   2.75935000e-03,   2.76173000e-03,
         2.76328000e-03,   2.76400000e-03,   2.76389000e-03,
         2.76295000e-03,   2.76118000e-03,   2.75858000e-03,
         2.75515000e-03,   2.75090000e-03,   2.74581000e-03,
         2.73989000e-03,   2.73315000e-03,   2.72557000e-03,
         2.71716000e-03,   2.70793000e-03,   2.69786000e-03,
         2.68697000e-03,   2.67525000e-03,   2.66269000e-03,
         2.64931000e-03,   2.63510000e-03,   2.62006000e-03,
         2.60418000e-03,   2.58748000e-03,   2.56995000e-03,
         2.55159000e-03,   2.53240000e-03,   2.51238000e-03,
         2.49153000e-03,   2.46985000e-03,   2.44735000e-03,
         2.42401000e-03,   2.39984000e-03,   2.37484000e-03,
         2.34902000e-03,   2.32236000e-03,   2.29488000e-03,
         2.26656000e-03,   2.23742000e-03,   2.20744000e-03,
         2.17664000e-03,   2.14500000e-03,   2.11254000e-03,
         2.07925000e-03,   2.04513000e-03,   2.01017000e-03,
         1.97439000e-03,   1.93778000e-03,   1.90034000e-03,
         1.86207000e-03,   1.82297000e-03,   1.78304000e-03,
         1.74229000e-03,   1.70070000e-03,   1.65828000e-03,
         1.61503000e-03,   1.57096000e-03,   1.52605000e-03,
         1.48031000e-03,   1.43375000e-03,   1.38635000e-03,
         1.33813000e-03,   1.28907000e-03,   1.23919000e-03,
         1.18848000e-03,   1.13693000e-03,   1.08456000e-03,
         1.03136000e-03,   9.77330000e-04,   9.22470000e-04,
         8.66770000e-04,   8.10250000e-04,   7.52900000e-04,
         6.94730000e-04,   6.35720000e-04,   5.75880000e-04,
         5.15210000e-04,   4.53710000e-04,   3.91380000e-04,
         3.28230000e-04,   2.64240000e-04,   1.99430000e-04,
         1.33780000e-04,   6.73000000e-05,   0.00000000e+00])
PpMuI = array([ 0.2990625,  0.2971875,  0.2953125,  0.2934375,  0.2915625,
        0.2896875,  0.2878125,  0.2859375,  0.2840625,  0.2821875,
        0.2803125,  0.2784375,  0.2765625,  0.2746875,  0.2728125,
        0.2709375,  0.2690625,  0.2671875,  0.2653125,  0.2634375,
        0.2615625,  0.2596875,  0.2578125,  0.2559375,  0.2540625,
        0.2521875,  0.2503125,  0.2484375,  0.2465625,  0.2446875,
        0.2428125,  0.2409375,  0.2390625,  0.2371875,  0.2353125,
        0.2334375,  0.2315625,  0.2296875,  0.2278125,  0.2259375,
        0.2240625,  0.2221875,  0.2203125,  0.2184375,  0.2165625,
        0.2146875,  0.2128125,  0.2109375,  0.2090625,  0.2071875,
        0.2053125,  0.2034375,  0.2015625,  0.1996875,  0.1978125,
        0.1959375,  0.1940625,  0.1921875,  0.1903125,  0.1884375,
        0.1865625,  0.1846875,  0.1828125,  0.1809375,  0.1790625,
        0.1771875,  0.1753125,  0.1734375,  0.1715625,  0.1696875,
        0.1678125,  0.1659375,  0.1640625,  0.1621875,  0.1603125,
        0.1584375,  0.1565625,  0.1546875,  0.1528125,  0.1509375,
        0.1490625,  0.1471875,  0.1453125,  0.1434375,  0.1415625,
        0.1396875,  0.1378125,  0.1359375,  0.1340625,  0.1321875,
        0.1303125,  0.1284375,  0.1265625,  0.1246875,  0.1228125,
        0.1209375,  0.1190625,  0.1171875,  0.1153125,  0.1134375,
        0.1115625,  0.1096875,  0.1078125,  0.1059375,  0.1040625,
        0.1021875,  0.1003125,  0.0984375,  0.0965625,  0.0946875,
        0.0928125,  0.0909375,  0.0890625,  0.0871875,  0.0853125,
        0.0834375,  0.0815625,  0.0796875,  0.0778125,  0.0759375,
        0.0740625,  0.0721875,  0.0703125,  0.0684375,  0.0665625,
        0.0646875,  0.0628125,  0.0609375,  0.0590625,  0.0571875,
        0.0553125,  0.0534375,  0.0515625,  0.0496875,  0.0478125,
        0.0459375,  0.0440625,  0.0421875,  0.0403125,  0.0384375,
        0.0365625,  0.0346875,  0.0328125,  0.0309375,  0.0290625,
        0.0271875,  0.0253125,  0.0234375,  0.0215625,  0.0196875,
        0.0178125,  0.0159375,  0.0140625,  0.0121875,  0.0103125,
        0.0084375,  0.0065625,  0.0046875,  0.0028125,  0.0009375,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,  0.       ])
