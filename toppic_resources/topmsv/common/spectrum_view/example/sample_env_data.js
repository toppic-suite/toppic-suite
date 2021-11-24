"use strict";
let sampleEnvelopes = [
    {
        "displayColor_": "red",
        "displayLevel_": 7,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 198.82540980333332,
                "monoMz_": 198.82540980333332,
                "intensity_": 2247.436393022484
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 199.15874313666666,
                "monoMz_": 199.15874313666666,
                "intensity_": 706.0636069775157
            }
        ],
        "monoMass_": 593.4544,
        "charge_": 3,
        "intensity_": 2793.16
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 199.23765424777778,
                "monoMz_": 199.23765424777778,
                "intensity_": 1978.441545389315
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 199.34876535888887,
                "monoMz_": 199.34876535888887,
                "intensity_": 1887.7984546106852
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 199.45987646999998,
                "monoMz_": 199.45987646999998,
                "intensity_": 984.7411259943406
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 199.5709875811111,
                "monoMz_": 199.5709875811111,
                "intensity_": 366.7288907416424
            }
        ],
        "monoMass_": 1784.0734,
        "charge_": 9,
        "intensity_": 3737.81
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 3,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 356.76062647000003,
                "monoMz_": 356.76062647000003,
                "intensity_": 5072.519040037839
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 357.26062647000003,
                "monoMz_": 357.26062647000003,
                "intensity_": 1894.9809599621613
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 357.76062647000003,
                "monoMz_": 357.76062647000003,
                "intensity_": 438.3042976543544
            }
        ],
        "monoMass_": 711.5067,
        "charge_": 2,
        "intensity_": 7359.43
    },
    {
        "displayColor_": "red",
        "displayLevel_": 6,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 365.76587647,
                "monoMz_": 365.76587647,
                "intensity_": 15296.814960750255
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 366.26587647,
                "monoMz_": 366.26587647,
                "intensity_": 5881.757363921583
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 366.76587647,
                "monoMz_": 366.76587647,
                "intensity_": 1384.2450392497474
            }
        ],
        "monoMass_": 729.5172,
        "charge_": 2,
        "intensity_": 23177.39
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 4,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 366.23447647,
                "monoMz_": 366.23447647,
                "intensity_": 1319.9666795240848
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 366.5678098033333,
                "monoMz_": 366.5678098033333,
                "intensity_": 766.5333204759153
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 366.90114313666663,
                "monoMz_": 366.90114313666663,
                "intensity_": 256.7243453990118
            }
        ],
        "monoMass_": 1095.6816,
        "charge_": 3,
        "intensity_": 2302.72
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 3,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 404.57770980333333,
                "monoMz_": 404.57770980333333,
                "intensity_": 1074.3867697622882
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 404.91104313666665,
                "monoMz_": 404.91104313666665,
                "intensity_": 687.7532302377119
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 405.24437647,
                "monoMz_": 405.24437647,
                "intensity_": 252.02067732131005
            }
        ],
        "monoMass_": 1210.7113,
        "charge_": 3,
        "intensity_": 1786.02
    },
    {
        "displayColor_": "red",
        "displayLevel_": 5,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 413.90834313666664,
                "monoMz_": 413.90834313666664,
                "intensity_": 4780.008965504166
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 414.24167647,
                "monoMz_": 414.24167647,
                "intensity_": 3181.814832091344
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 414.57500980333333,
                "monoMz_": 414.57500980333333,
                "intensity_": 1200.2910344958339
            }
        ],
        "monoMass_": 1238.7032,
        "charge_": 3,
        "intensity_": 11233.24
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 3,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 455.91840980333336,
                "monoMz_": 455.91840980333336,
                "intensity_": 2172.8756087722345
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 456.25174313666673,
                "monoMz_": 456.25174313666673,
                "intensity_": 1575.7243912277663
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 456.58507647000005,
                "monoMz_": 456.58507647000005,
                "intensity_": 643.854267101711
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 456.91840980333336,
                "monoMz_": 456.91840980333336,
                "intensity_": 190.71813769454656
            }
        ],
        "monoMass_": 1364.7334,
        "charge_": 3,
        "intensity_": 4776.39
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 5,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 461.92107647000006,
                "monoMz_": 461.92107647000006,
                "intensity_": 4215.704203933385
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 462.2544098033334,
                "monoMz_": 462.2544098033334,
                "intensity_": 3103.22214802853
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 462.58774313666675,
                "monoMz_": 462.58774313666675,
                "intensity_": 1282.5957960666142
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 462.92107647000006,
                "monoMz_": 462.92107647000006,
                "intensity_": 383.6798396365533
            }
        ],
        "monoMass_": 1382.7414,
        "charge_": 3,
        "intensity_": 8455.83
    },
    {
        "displayColor_": "red",
        "displayLevel_": 4,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 466.9277431366667,
                "monoMz_": 466.9277431366667,
                "intensity_": 2201.9276981094117
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 467.0944098033334,
                "monoMz_": 467.0944098033334,
                "intensity_": 1827.3980946801719
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 466.76107647000003,
                "monoMz_": 466.76107647000003,
                "intensity_": 1458.4428253536353
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 467.26107647000003,
                "monoMz_": 467.26107647000003,
                "intensity_": 1085.5923018905887
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 467.4277431366667,
                "monoMz_": 467.4277431366667,
                "intensity_": 511.29456959253133
            }
        ],
        "monoMass_": 2794.5228,
        "charge_": 6,
        "intensity_": 6131.37
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 470.59127647,
                "monoMz_": 470.59127647,
                "intensity_": 1176.353277765729
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 470.92460980333334,
                "monoMz_": 470.92460980333334,
                "intensity_": 883.2167222342713
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 471.25794313666665,
                "monoMz_": 471.25794313666665,
                "intensity_": 370.6751995504845
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 471.59127647,
                "monoMz_": 471.59127647,
                "intensity_": 112.36018324602246
            }
        ],
        "monoMass_": 1408.752,
        "charge_": 3,
        "intensity_": 2186.59
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 480.2537193271429,
                "monoMz_": 480.2537193271429,
                "intensity_": 1087.1251804807491
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 480.39657647,
                "monoMz_": 480.39657647,
                "intensity_": 1970.847256658574
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 480.5394336128572,
                "monoMz_": 480.5394336128572,
                "intensity_": 1923.5990711171794
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 480.68229075571435,
                "monoMz_": 480.68229075571435,
                "intensity_": 1327.9226676281714
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 480.82514789857146,
                "monoMz_": 480.82514789857146,
                "intensity_": 721.6527433414257
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 480.9680050414286,
                "monoMz_": 480.9680050414286,
                "intensity_": 326.6872668027465
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 481.11086218428574,
                "monoMz_": 481.11086218428574,
                "intensity_": 127.54271012659237
            }
        ],
        "monoMass_": 3354.7251,
        "charge_": 7,
        "intensity_": 6025.54
    },
    {
        "displayColor_": "red",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 480.27197647,
                "monoMz_": 480.27197647,
                "intensity_": 794.601536647581
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 480.4386431366667,
                "monoMz_": 480.4386431366667,
                "intensity_": 1229.118463352419
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 480.6053098033334,
                "monoMz_": 480.6053098033334,
                "intensity_": 1042.108138318087
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 480.77197647,
                "monoMz_": 480.77197647,
                "intensity_": 631.4885562871173
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 480.9386431366667,
                "monoMz_": 480.9386431366667,
                "intensity_": 303.128366557609
            }
        ],
        "monoMass_": 2875.5882,
        "charge_": 6,
        "intensity_": 3408.01
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 4,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 483.43902647,
                "monoMz_": 483.43902647,
                "intensity_": 3655.9682952621606
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 483.6056931366667,
                "monoMz_": 483.6056931366667,
                "intensity_": 3117.838278752273
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 483.2723598033333,
                "monoMz_": 483.2723598033333,
                "intensity_": 2346.7554829804076
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 483.7723598033333,
                "monoMz_": 483.7723598033333,
                "intensity_": 1899.0271429027123
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 483.93902647,
                "monoMz_": 483.93902647,
                "intensity_": 915.8598714579066
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 484.1056931366667,
                "monoMz_": 484.1056931366667,
                "intensity_": 369.21170473783945
            }
        ],
        "monoMass_": 2893.5905,
        "charge_": 6,
        "intensity_": 10203.04
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 486.47922647,
                "monoMz_": 486.47922647,
                "intensity_": 856.3240196080646
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 486.97922647,
                "monoMz_": 486.97922647,
                "intensity_": 443.1815220724283
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 487.47922647,
                "monoMz_": 487.47922647,
                "intensity_": 133.5812035969757
            }
        ],
        "monoMass_": 970.9439,
        "charge_": 2,
        "intensity_": 2411.68
    },
    {
        "displayColor_": "red",
        "displayLevel_": 3,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 486.48028647,
                "monoMz_": 486.48028647,
                "intensity_": 586.3263814751936
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 486.58028647,
                "monoMz_": 486.58028647,
                "intensity_": 1531.5963340421365
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 486.68028647,
                "monoMz_": 486.68028647,
                "intensity_": 2096.1305250876144
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 486.78028647,
                "monoMz_": 486.78028647,
                "intensity_": 1991.735532139224
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 486.88028647,
                "monoMz_": 486.88028647,
                "intensity_": 1470.8957688057444
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 486.98028647,
                "monoMz_": 486.98028647,
                "intensity_": 896.8680768443131
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 487.08028647,
                "monoMz_": 487.08028647,
                "intensity_": 468.7394749123854
            }
        ],
        "monoMass_": 4854.7301,
        "charge_": 10,
        "intensity_": 5359.37
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 1,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 488.77622647000004,
                "monoMz_": 488.77622647000004,
                "intensity_": 1972.6174139955892
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 489.02622647000004,
                "monoMz_": 489.02622647000004,
                "intensity_": 2050.2263756868842
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 489.27622647000004,
                "monoMz_": 489.27622647000004,
                "intensity_": 1156.5571577256333
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 489.52622647000004,
                "monoMz_": 489.52622647000004,
                "intensity_": 463.8517236748597
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 489.77622647000004,
                "monoMz_": 489.77622647000004,
                "intensity_": 147.15333743156182
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 490.02622647000004,
                "monoMz_": 490.02622647000004,
                "intensity_": 39.07362431311578
            }
        ],
        "monoMass_": 1951.0758,
        "charge_": 4,
        "intensity_": 3373.24
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 5,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 492.78122647000004,
                "monoMz_": 492.78122647000004,
                "intensity_": 5976.597326723034
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 493.03122647000004,
                "monoMz_": 493.03122647000004,
                "intensity_": 6301.175076711756
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 493.28122647000004,
                "monoMz_": 493.28122647000004,
                "intensity_": 3609.6338137518987
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 493.53122647000004,
                "monoMz_": 493.53122647000004,
                "intensity_": 1471.0424028275124
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 493.78122647000004,
                "monoMz_": 493.78122647000004,
                "intensity_": 474.4249232882432
            }
        ],
        "monoMass_": 1967.0958,
        "charge_": 4,
        "intensity_": 12218.46
    },
    {
        "displayColor_": "red",
        "displayLevel_": 1,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 497.28767647000006,
                "monoMz_": 497.28767647000006,
                "intensity_": 1932.1280193307684
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 497.6210098033334,
                "monoMz_": 497.6210098033334,
                "intensity_": 1543.376693646466
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 497.9543431366667,
                "monoMz_": 497.9543431366667,
                "intensity_": 684.1669373919092
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 498.28767647000006,
                "monoMz_": 498.28767647000006,
                "intensity_": 218.27198066923168
            }
        ],
        "monoMass_": 1488.8412,
        "charge_": 3,
        "intensity_": 3088.8
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 3,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 499.68371647000004,
                "monoMz_": 499.68371647000004,
                "intensity_": 2610.405868149625
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 499.78371647,
                "monoMz_": 499.78371647,
                "intensity_": 2541.1014979465363
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 499.88371647,
                "monoMz_": 499.88371647,
                "intensity_": 1921.048163620799
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 499.58371647,
                "monoMz_": 499.58371647,
                "intensity_": 1859.8935327461907
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 499.98371647000005,
                "monoMz_": 499.98371647000005,
                "intensity_": 1198.394131850375
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 499.48371647000005,
                "monoMz_": 499.48371647000005,
                "intensity_": 693.3343707243065
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 500.08371647,
                "monoMz_": 500.08371647,
                "intensity_": 640.5148702029346
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 500.18371647000004,
                "monoMz_": 500.18371647000004,
                "intensity_": 300.8654086125095
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 500.28371647,
                "monoMz_": 500.28371647,
                "intensity_": 126.48839102383077
            }
        ],
        "monoMass_": 4984.7644,
        "charge_": 10,
        "intensity_": 10559.51
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 4,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 505.27890980333336,
                "monoMz_": 505.27890980333336,
                "intensity_": 6133.372661082358
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 505.44557647000005,
                "monoMz_": 505.44557647000005,
                "intensity_": 5440.873487380698
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 505.11224313666673,
                "monoMz_": 505.11224313666673,
                "intensity_": 3766.903433730912
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 505.61224313666673,
                "monoMz_": 505.61224313666673,
                "intensity_": 3438.0327889821224
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 505.77890980333336,
                "monoMz_": 505.77890980333336,
                "intensity_": 1717.6273389176424
            }
        ],
        "monoMass_": 3024.6298,
        "charge_": 6,
        "intensity_": 12599.81
    },
    {
        "displayColor_": "red",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 505.27925147,
                "monoMz_": 505.27925147,
                "intensity_": 480.28799062905966
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 505.52925147,
                "monoMz_": 505.52925147,
                "intensity_": 516.9265504226884
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 505.77925147,
                "monoMz_": 505.77925147,
                "intensity_": 301.26344957731163
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 506.02925147,
                "monoMz_": 506.02925147,
                "intensity_": 124.65298653163065
            }
        ],
        "monoMass_": 2017.0879,
        "charge_": 4,
        "intensity_": 6256.65
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 506.29087647000006,
                "monoMz_": 506.29087647000006,
                "intensity_": 415.5526890808573
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 506.39087647,
                "monoMz_": 506.39087647,
                "intensity_": 1130.0859124745157
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 506.49087647000005,
                "monoMz_": 506.49087647000005,
                "intensity_": 1606.7989449941908
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 506.59087647000007,
                "monoMz_": 506.59087647000007,
                "intensity_": 1583.6973377744418
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 506.69087647000003,
                "monoMz_": 506.69087647000003,
                "intensity_": 1211.735674020216
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 506.79087647000006,
                "monoMz_": 506.79087647000006,
                "intensity_": 764.8110550058091
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 506.89087647,
                "monoMz_": 506.89087647,
                "intensity_": 413.49129457834545
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 506.99087647000005,
                "monoMz_": 506.99087647000005,
                "intensity_": 196.43325986416832
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 507.09087647000007,
                "monoMz_": 507.09087647000007,
                "intensity_": 83.51091676368722
            }
        ],
        "monoMass_": 5052.836,
        "charge_": 10,
        "intensity_": 5470.44
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 508.08040980333334,
                "monoMz_": 508.08040980333334,
                "intensity_": 458.24311634169453
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 508.1915209144444,
                "monoMz_": 508.1915209144444,
                "intensity_": 1124.0858123801117
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 508.30263202555557,
                "monoMz_": 508.30263202555557,
                "intensity_": 1450.114187619888
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 508.41374313666665,
                "monoMz_": 508.41374313666665,
                "intensity_": 1302.4015541537683
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 508.5248542477778,
                "monoMz_": 508.5248542477778,
                "intensity_": 910.995549325669
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 508.6359653588889,
                "monoMz_": 508.6359653588889,
                "intensity_": 526.9249906190703
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 508.74707647,
                "monoMz_": 508.74707647,
                "intensity_": 261.53753398060826
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 508.8581875811111,
                "monoMz_": 508.8581875811111,
                "intensity_": 114.21905605182555
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 508.96929869222225,
                "monoMz_": 508.96929869222225,
                "intensity_": 44.681585401470294
            }
        ],
        "monoMass_": 4563.6582,
        "charge_": 9,
        "intensity_": 3702.48
    },
    {
        "displayColor_": "red",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 540.19697647,
                "monoMz_": 540.19697647,
                "intensity_": 488.4586336216338
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 540.3080875811111,
                "monoMz_": 540.3080875811111,
                "intensity_": 1275.9471110678262
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 540.4191986922222,
                "monoMz_": 540.4191986922222,
                "intensity_": 1746.2510378619425
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 540.5303098033334,
                "monoMz_": 540.5303098033334,
                "intensity_": 1659.2813274351088
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 540.6414209144444,
                "monoMz_": 540.6414209144444,
                "intensity_": 1225.378492475515
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 540.7525320255555,
                "monoMz_": 540.7525320255555,
                "intensity_": 747.1656899558599
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 540.8636431366666,
                "monoMz_": 540.8636431366666,
                "intensity_": 390.49896213805755
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 540.9747542477777,
                "monoMz_": 540.9747542477777,
                "intensity_": 179.42182837456608
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 541.0858653588889,
                "monoMz_": 541.0858653588889,
                "intensity_": 73.80098687517147
            }
        ],
        "monoMass_": 4852.7073,
        "charge_": 9,
        "intensity_": 5142.61
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 4,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 543.3122542477778,
                "monoMz_": 543.3122542477778,
                "intensity_": 1059.1052189162376
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 543.4233653588889,
                "monoMz_": 543.4233653588889,
                "intensity_": 2789.6154683617237
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 543.53447647,
                "monoMz_": 543.53447647,
                "intensity_": 3846.614781083762
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 543.6455875811112,
                "monoMz_": 543.6455875811112,
                "intensity_": 3680.4223664592546
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 543.7566986922222,
                "monoMz_": 543.7566986922222,
                "intensity_": 2735.629153096761
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 543.8678098033333,
                "monoMz_": 543.8678098033333,
                "intensity_": 1678.2576419285053
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 543.9789209144444,
                "monoMz_": 543.9789209144444,
                "intensity_": 882.2506805092273
            }
        ],
        "monoMass_": 4880.7448,
        "charge_": 9,
        "intensity_": 14599.78
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 543.8622478985715,
                "monoMz_": 543.8622478985715,
                "intensity_": 2439.958993500933
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 543.7193907557142,
                "monoMz_": 543.7193907557142,
                "intensity_": 2241.6028199193543
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 544.0051050414286,
                "monoMz_": 544.0051050414286,
                "intensity_": 1866.194242720251
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 544.1479621842857,
                "monoMz_": 544.1479621842857,
                "intensity_": 1118.8846513731523
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 543.5765336128571,
                "monoMz_": 543.5765336128571,
                "intensity_": 1097.4720592379867
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 544.2908193271428,
                "monoMz_": 544.2908193271428,
                "intensity_": 557.2833645705618
            }
        ],
        "monoMass_": 3797.9848,
        "charge_": 7,
        "intensity_": 10257.47
    },
    {
        "displayColor_": "red",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 551.87397647,
                "monoMz_": 551.87397647,
                "intensity_": 516.2428874031338
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 551.9850875811111,
                "monoMz_": 551.9850875811111,
                "intensity_": 1378.940536680728
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 552.0961986922222,
                "monoMz_": 552.0961986922222,
                "intensity_": 1926.8346897247059
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 552.2073098033333,
                "monoMz_": 552.2073098033333,
                "intensity_": 1867.1810052949493
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 552.3184209144443,
                "monoMz_": 552.3184209144443,
                "intensity_": 1405.0343293677336
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 552.4295320255555,
                "monoMz_": 552.4295320255555,
                "intensity_": 872.3519966854191
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 552.5406431366666,
                "monoMz_": 552.5406431366666,
                "intensity_": 464.00726604030734
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 552.6517542477777,
                "monoMz_": 552.6517542477777,
                "intensity_": 216.88595780143018
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 552.7628653588888,
                "monoMz_": 552.7628653588888,
                "intensity_": 90.72585917295125
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 552.87397647,
                "monoMz_": 552.87397647,
                "intensity_": 34.42531027529403
            }
        ],
        "monoMass_": 4957.8003,
        "charge_": 9,
        "intensity_": 4615.72
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 6,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 555.0922320255555,
                "monoMz_": 555.0922320255555,
                "intensity_": 32509.252754471
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 555.2033431366667,
                "monoMz_": 555.2033431366667,
                "intensity_": 31646.155825594305
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 555.3144542477778,
                "monoMz_": 555.3144542477778,
                "intensity_": 23924.187830963478
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 554.9811209144444,
                "monoMz_": 554.9811209144444,
                "intensity_": 23162.585439371298
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 555.4255653588889,
                "monoMz_": 555.4255653588889,
                "intensity_": 14924.45991144455
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 554.8700098033333,
                "monoMz_": 554.8700098033333,
                "intensity_": 8634.589194061153
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 555.53667647,
                "monoMz_": 555.53667647,
                "intensity_": 7976.790146884108
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 555.6477875811112,
                "monoMz_": 555.6477875811112,
                "intensity_": 3746.8922871348054
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 555.7588986922223,
                "monoMz_": 555.7588986922223,
                "intensity_": 1575.250471381632
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 555.8700098033333,
                "monoMz_": 555.8700098033333,
                "intensity_": 600.7872455290013
            }
        ],
        "monoMass_": 4984.7646,
        "charge_": 9,
        "intensity_": 135674.42
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 556.30881647,
                "monoMz_": 556.30881647,
                "intensity_": 2379.310736425231
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 556.5088164699999,
                "monoMz_": 556.5088164699999,
                "intensity_": 3556.3554637628936
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 556.70881647,
                "monoMz_": 556.70881647,
                "intensity_": 2922.592151004438
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 556.90881647,
                "monoMz_": 556.90881647,
                "intensity_": 1719.4154659732637
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 557.10881647,
                "monoMz_": 557.10881647,
                "intensity_": 801.9445362371065
            }
        ],
        "monoMass_": 2776.5077,
        "charge_": 5,
        "intensity_": 9759.65
    },
    {
        "displayColor_": "red",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 556.8026764699999,
                "monoMz_": 556.8026764699999,
                "intensity_": 2432.0237594390096
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 557.8026764699999,
                "monoMz_": 557.8026764699999,
                "intensity_": 701.7220541309647
            }
        ],
        "monoMass_": 555.7954,
        "charge_": 1,
        "intensity_": 1581.36
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 557.2941764699999,
                "monoMz_": 557.2941764699999,
                "intensity_": 3670.7396490241135
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 557.4608431366667,
                "monoMz_": 557.4608431366667,
                "intensity_": 6601.142664972139
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 557.6275098033333,
                "monoMz_": 557.6275098033333,
                "intensity_": 6398.717436945089
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 557.7941764699999,
                "monoMz_": 557.7941764699999,
                "intensity_": 4390.2147569475155
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 557.9608431366667,
                "monoMz_": 557.9608431366667,
                "intensity_": 2372.4100107521704
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 558.1275098033333,
                "monoMz_": 558.1275098033333,
                "intensity_": 1068.2981675742724
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 558.2941764699999,
                "monoMz_": 558.2941764699999,
                "intensity_": 414.97733502786105
            }
        ],
        "monoMass_": 3337.7214,
        "charge_": 6,
        "intensity_": 20303.89
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 557.7992209144444,
                "monoMz_": 557.7992209144444,
                "intensity_": 55.65418643002552
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 557.8547764699999,
                "monoMz_": 557.8547764699999,
                "intensity_": 301.2868626030358
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 557.9103320255555,
                "monoMz_": 557.9103320255555,
                "intensity_": 836.6982429914208
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 557.965887581111,
                "monoMz_": 557.965887581111,
                "intensity_": 1586.246150603477
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 558.0214431366666,
                "monoMz_": 558.0214431366666,
                "intensity_": 2305.697417399463
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 558.0769986922221,
                "monoMz_": 558.0769986922221,
                "intensity_": 2736.769626085633
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 558.1325542477776,
                "monoMz_": 558.1325542477776,
                "intensity_": 2759.4458135699742
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 558.1881098033333,
                "monoMz_": 558.1881098033333,
                "intensity_": 2428.0717750500085
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 558.2436653588888,
                "monoMz_": 558.2436653588888,
                "intensity_": 1901.241056985668
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 558.2992209144444,
                "monoMz_": 558.2992209144444,
                "intensity_": 1344.489750748163
            },
            {
                "displayLevel_": -1,
                "peakId_": "10",
                "pos_": 558.3547764699999,
                "monoMz_": 558.3547764699999,
                "intensity_": 868.6094125911205
            },
            {
                "displayLevel_": -1,
                "peakId_": "11",
                "pos_": 558.4103320255555,
                "monoMz_": 558.4103320255555,
                "intensity_": 517.4169514547207
            }
        ],
        "monoMass_": 10022.255,
        "charge_": 18,
        "intensity_": 4746.9
    },
    {
        "displayColor_": "red",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 557.9278389699999,
                "monoMz_": 557.9278389699999,
                "intensity_": 1410.5302142987425
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 558.0528389699999,
                "monoMz_": 558.0528389699999,
                "intensity_": 3392.2289585681424
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 558.1778389699999,
                "monoMz_": 558.1778389699999,
                "intensity_": 4296.270741127661
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 558.3028389699999,
                "monoMz_": 558.3028389699999,
                "intensity_": 3792.0680543108374
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 558.4278389699999,
                "monoMz_": 558.4278389699999,
                "intensity_": 2608.623459626972
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 558.5528389699999,
                "monoMz_": 558.5528389699999,
                "intensity_": 1484.7292588723399
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 558.6778389699999,
                "monoMz_": 558.6778389699999,
                "intensity_": 725.463508486294
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 558.8028389699999,
                "monoMz_": 558.8028389699999,
                "intensity_": 311.98740497064927
            }
        ],
        "monoMass_": 4455.3645,
        "charge_": 8,
        "intensity_": 13759.93
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 560.1118564699999,
                "monoMz_": 560.1118564699999,
                "intensity_": 12266.641335953955
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 560.31185647,
                "monoMz_": 560.31185647,
                "intensity_": 10180.18758049769
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 559.91185647,
                "monoMz_": 559.91185647,
                "intensity_": 8124.787686248284
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 560.51185647,
                "monoMz_": 560.51185647,
                "intensity_": 6047.687858142749
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 560.7118564699999,
                "monoMz_": 560.7118564699999,
                "intensity_": 2848.35288079513
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 560.91185647,
                "monoMz_": 560.91185647,
                "intensity_": 1122.1686640460446
            }
        ],
        "monoMass_": 2794.5229,
        "charge_": 5,
        "intensity_": 35088.79
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 3,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 560.4620598033333,
                "monoMz_": 560.4620598033333,
                "intensity_": 25181.859717740655
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 560.62872647,
                "monoMz_": 560.62872647,
                "intensity_": 24578.16139652301
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 560.7953931366667,
                "monoMz_": 560.7953931366667,
                "intensity_": 16967.099920728953
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 560.2953931366667,
                "monoMz_": 560.2953931366667,
                "intensity_": 13890.388358610511
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 560.9620598033333,
                "monoMz_": 560.9620598033333,
                "intensity_": 9220.68317894747
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 561.12872647,
                "monoMz_": 561.12872647,
                "intensity_": 4174.1402822593445
            }
        ],
        "monoMass_": 3355.7287,
        "charge_": 6,
        "intensity_": 83840.49
    },
    {
        "displayColor_": "red",
        "displayLevel_": 3,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 563.0145898033334,
                "monoMz_": 563.0145898033334,
                "intensity_": 62.16062913493153
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 563.08125647,
                "monoMz_": 563.08125647,
                "intensity_": 283.16345215256155
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 563.1479231366667,
                "monoMz_": 563.1479231366667,
                "intensity_": 666.1956675625403
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 563.2145898033334,
                "monoMz_": 563.2145898033334,
                "intensity_": 1076.0507727005686
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 563.28125647,
                "monoMz_": 563.28125647,
                "intensity_": 1338.9276094941797
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 563.3479231366667,
                "monoMz_": 563.3479231366667,
                "intensity_": 1365.9383157406435
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 563.4145898033333,
                "monoMz_": 563.4145898033333,
                "intensity_": 1187.7790631748903
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 563.4812564700001,
                "monoMz_": 563.4812564700001,
                "intensity_": 903.9840011451303
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 563.5479231366667,
                "monoMz_": 563.5479231366667,
                "intensity_": 613.7660239670945
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 563.6145898033334,
                "monoMz_": 563.6145898033334,
                "intensity_": 377.15017085658695
            },
            {
                "displayLevel_": -1,
                "peakId_": "10",
                "pos_": 563.68125647,
                "monoMz_": 563.68125647,
                "intensity_": 212.10969038163634
            },
            {
                "displayLevel_": -1,
                "peakId_": "11",
                "pos_": 563.7479231366667,
                "monoMz_": 563.7479231366667,
                "intensity_": 110.15821267729359
            },
            {
                "displayLevel_": -1,
                "peakId_": "12",
                "pos_": 563.8145898033333,
                "monoMz_": 563.8145898033333,
                "intensity_": 53.21168425935671
            }
        ],
        "monoMass_": 8430.1097,
        "charge_": 15,
        "intensity_": 3195.7
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 567.6579542477778,
                "monoMz_": 567.6579542477778,
                "intensity_": 206.35884928172783
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 567.7690653588888,
                "monoMz_": 567.7690653588888,
                "intensity_": 566.5006604757391
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 567.8801764699999,
                "monoMz_": 567.8801764699999,
                "intensity_": 812.4089077366981
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 567.9912875811111,
                "monoMz_": 567.9912875811111,
                "intensity_": 807.1051611752515
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 568.1023986922222,
                "monoMz_": 568.1023986922222,
                "intensity_": 622.1472016696669
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 568.2135098033333,
                "monoMz_": 568.2135098033333,
                "intensity_": 395.4518726386235
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 568.3246209144444,
                "monoMz_": 568.3246209144444,
                "intensity_": 215.23832428559749
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 568.4357320255555,
                "monoMz_": 568.4357320255555,
                "intensity_": 102.91203649706054
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 568.5468431366667,
                "monoMz_": 568.5468431366667,
                "intensity_": 44.024739301547534
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 568.6579542477778,
                "monoMz_": 568.6579542477778,
                "intensity_": 17.081092263301933
            }
        ],
        "monoMass_": 5099.8561,
        "charge_": 9,
        "intensity_": 2617.21
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 3,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 568.3292098033332,
                "monoMz_": 568.3292098033332,
                "intensity_": 2759.4753254935777
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 568.49587647,
                "monoMz_": 568.49587647,
                "intensity_": 5064.342734094482
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 568.6625431366666,
                "monoMz_": 568.6625431366666,
                "intensity_": 5000.599573976964
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 568.8292098033332,
                "monoMz_": 568.8292098033332,
                "intensity_": 3490.9320203042257
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 568.99587647,
                "monoMz_": 568.99587647,
                "intensity_": 1918.057265905519
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 569.1625431366666,
                "monoMz_": 569.1625431366666,
                "intensity_": 877.7888017317872
            }
        ],
        "monoMass_": 3403.9316,
        "charge_": 6,
        "intensity_": 17087.81
    },
    {
        "displayColor_": "red",
        "displayLevel_": 1,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 570.3144431366666,
                "monoMz_": 570.3144431366666,
                "intensity_": 2673.562410174301
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 570.6477764699999,
                "monoMz_": 570.6477764699999,
                "intensity_": 2423.082764514855
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 570.9811098033333,
                "monoMz_": 570.9811098033333,
                "intensity_": 1206.8152992493383
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 571.3144431366666,
                "monoMz_": 571.3144431366666,
                "intensity_": 430.3492514918556
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 571.6477764699999,
                "monoMz_": 571.6477764699999,
                "intensity_": 121.95596632817814
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 571.9811098033333,
                "monoMz_": 571.9811098033333,
                "intensity_": 29.01758982569908
            }
        ],
        "monoMass_": 1707.9215,
        "charge_": 3,
        "intensity_": 4264.48
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 570.4851639699999,
                "monoMz_": 570.4851639699999,
                "intensity_": 142.83109229119495
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 570.6101639699999,
                "monoMz_": 570.6101639699999,
                "intensity_": 350.36948441919907
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 570.7351639699999,
                "monoMz_": 570.7351639699999,
                "intensity_": 451.9901903125694
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 570.8601639699999,
                "monoMz_": 570.8601639699999,
                "intensity_": 405.94922203440575
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 570.9851639699999,
                "monoMz_": 570.9851639699999,
                "intensity_": 283.95077796559417
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 571.1101639699999,
                "monoMz_": 571.1101639699999,
                "intensity_": 164.2387398341844
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 571.2351639699999,
                "monoMz_": 571.2351639699999,
                "intensity_": 81.51937327901082
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 571.3601639699999,
                "monoMz_": 571.3601639699999,
                "intensity_": 35.601260454476105
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 571.4851639699999,
                "monoMz_": 571.4851639699999,
                "intensity_": 13.926929659398432
            }
        ],
        "monoMass_": 4555.8231,
        "charge_": 8,
        "intensity_": 2516.21
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 571.57869647,
                "monoMz_": 571.57869647,
                "intensity_": 262.84181812200023
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 571.67869647,
                "monoMz_": 571.67869647,
                "intensity_": 809.1828149929452
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 571.77869647,
                "monoMz_": 571.77869647,
                "intensity_": 1305.5731897026158
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 571.87869647,
                "monoMz_": 571.87869647,
                "intensity_": 1462.8045834931524
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 571.97869647,
                "monoMz_": 571.97869647,
                "intensity_": 1274.0237422628438
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 572.07869647,
                "monoMz_": 572.07869647,
                "intensity_": 916.2686679114573
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 572.17869647,
                "monoMz_": 572.17869647,
                "intensity_": 564.902114645606
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 572.27869647,
                "monoMz_": 572.27869647,
                "intensity_": 306.210448663526
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 572.37869647,
                "monoMz_": 572.37869647,
                "intensity_": 148.6123151358617
            }
        ],
        "monoMass_": 5705.7142,
        "charge_": 10,
        "intensity_": 3830
    },
    {
        "displayColor_": "red",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 572.6400875811111,
                "monoMz_": 572.6400875811111,
                "intensity_": 1188.4533180386481
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 572.7511986922223,
                "monoMz_": 572.7511986922223,
                "intensity_": 3289.138665790833
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 572.8623098033333,
                "monoMz_": 572.8623098033333,
                "intensity_": 4754.335449125524
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 572.9734209144444,
                "monoMz_": 572.9734209144444,
                "intensity_": 4760.045313881337
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 573.0845320255555,
                "monoMz_": 573.0845320255555,
                "intensity_": 3697.3526786381535
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 573.1956431366667,
                "monoMz_": 573.1956431366667,
                "intensity_": 2367.9603689121213
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 573.3067542477778,
                "monoMz_": 573.3067542477778,
                "intensity_": 1298.5608298216785
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 573.4178653588889,
                "monoMz_": 573.4178653588889,
                "intensity_": 625.5470669780926
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 573.52897647,
                "monoMz_": 573.52897647,
                "intensity_": 269.61248901177123
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 573.6400875811111,
                "monoMz_": 573.6400875811111,
                "intensity_": 105.39468611866305
            }
        ],
        "monoMass_": 5144.6953,
        "charge_": 9,
        "intensity_": 16816.8
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 573.6223401063636,
                "monoMz_": 573.6223401063636,
                "intensity_": 184.45467394754868
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 573.7132491972727,
                "monoMz_": 573.7132491972727,
                "intensity_": 628.56988162887
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 573.8041582881817,
                "monoMz_": 573.8041582881817,
                "intensity_": 1115.8338824495813
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 573.8950673790908,
                "monoMz_": 573.8950673790908,
                "intensity_": 1369.2661175504193
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 573.98597647,
                "monoMz_": 573.98597647,
                "intensity_": 1301.5313844813857
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 574.0768855609091,
                "monoMz_": 574.0768855609091,
                "intensity_": 1018.8046455891775
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 574.1677946518181,
                "monoMz_": 574.1677946518181,
                "intensity_": 682.1942725859018
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 574.2587037427272,
                "monoMz_": 574.2587037427272,
                "intensity_": 400.9579387846637
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 574.3496128336363,
                "monoMz_": 574.3496128336363,
                "intensity_": 210.72503028435813
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 574.4405219245455,
                "monoMz_": 574.4405219245455,
                "intensity_": 100.41555520305644
            }
        ],
        "monoMass_": 6298.7657,
        "charge_": 11,
        "intensity_": 3943.89
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 573.62935147,
                "monoMz_": 573.62935147,
                "intensity_": 1333.447275713785
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 573.87935147,
                "monoMz_": 573.87935147,
                "intensity_": 1632.8207270043495
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 574.12935147,
                "monoMz_": 574.12935147,
                "intensity_": 1073.7792729956504
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 574.37935147,
                "monoMz_": 574.37935147,
                "intensity_": 498.78217242218875
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 574.62935147,
                "monoMz_": 574.62935147,
                "intensity_": 182.45316781005843
            }
        ],
        "monoMass_": 2290.4883,
        "charge_": 4,
        "intensity_": 2809.82
    },
    {
        "displayColor_": "red",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 574.6418542477777,
                "monoMz_": 574.6418542477777,
                "intensity_": 5837.70458692292
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 574.7529653588888,
                "monoMz_": 574.7529653588888,
                "intensity_": 16241.444745247467
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 574.86407647,
                "monoMz_": 574.86407647,
                "intensity_": 23589.246235615905
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 574.9751875811111,
                "monoMz_": 574.9751875811111,
                "intensity_": 23722.68085923249
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 575.0862986922222,
                "monoMz_": 575.0862986922222,
                "intensity_": 18503.420631639547
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 575.1974098033332,
                "monoMz_": 575.1974098033332,
                "intensity_": 11897.294524726498
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 575.3085209144443,
                "monoMz_": 575.3085209144443,
                "intensity_": 6548.936891258463
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 575.4196320255555,
                "monoMz_": 575.4196320255555,
                "intensity_": 3166.2033491774837
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 575.5307431366666,
                "monoMz_": 575.5307431366666,
                "intensity_": 1369.4171358677147
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 575.6418542477777,
                "monoMz_": 575.6418542477777,
                "intensity_": 537.1391407675114
            }
        ],
        "monoMass_": 5162.7112,
        "charge_": 9,
        "intensity_": 87181.44
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 576.64272647,
                "monoMz_": 576.64272647,
                "intensity_": 4731.471333786537
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 577.14272647,
                "monoMz_": 577.14272647,
                "intensity_": 2922.455630105434
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 577.64272647,
                "monoMz_": 577.64272647,
                "intensity_": 1033.8304608682788
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 578.14272647,
                "monoMz_": 578.14272647,
                "intensity_": 266.9597853156041
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 578.64272647,
                "monoMz_": 578.64272647,
                "intensity_": 55.42866621346257
            }
        ],
        "monoMass_": 1151.2709,
        "charge_": 2,
        "intensity_": 23960.12
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 7,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 576.9757320255554,
                "monoMz_": 576.9757320255554,
                "intensity_": 52167.593582559435
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 576.8646209144443,
                "monoMz_": 576.8646209144443,
                "intensity_": 51702.56747995701
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 577.0868431366666,
                "monoMz_": 577.0868431366666,
                "intensity_": 40816.740598482196
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 576.7535098033333,
                "monoMz_": 576.7535098033333,
                "intensity_": 35470.661427466635
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 577.1979542477777,
                "monoMz_": 577.1979542477777,
                "intensity_": 26321.594941486772
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 577.3090653588888,
                "monoMz_": 577.3090653588888,
                "intensity_": 14529.597135797341
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 576.6423986922222,
                "monoMz_": 576.6423986922222,
                "intensity_": 12699.421795501905
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 577.4201764699999,
                "monoMz_": 577.4201764699999,
                "intensity_": 7043.562585302201
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 577.5312875811111,
                "monoMz_": 577.5312875811111,
                "intensity_": 3054.354176554042
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 577.6423986922222,
                "monoMz_": 577.6423986922222,
                "intensity_": 1201.0764174405704
            }
        ],
        "monoMass_": 5180.7161,
        "charge_": 9,
        "intensity_": 209861.94
    },
    {
        "displayColor_": "red",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 578.08677647,
                "monoMz_": 578.08677647,
                "intensity_": 2546.6914428319324
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 578.33677647,
                "monoMz_": 578.33677647,
                "intensity_": 3146.2885571680677
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 578.58677647,
                "monoMz_": 578.58677647,
                "intensity_": 2084.852936901351
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 578.83677647,
                "monoMz_": 578.83677647,
                "intensity_": 975.0218926204142
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 579.08677647,
                "monoMz_": 579.08677647,
                "intensity_": 358.87432512412204
            }
        ],
        "monoMass_": 2308.318,
        "charge_": 4,
        "intensity_": 11376.32
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 1,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 579.92557647,
                "monoMz_": 579.92557647,
                "intensity_": 12786.613957436228
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 580.1255764699999,
                "monoMz_": 580.1255764699999,
                "intensity_": 10904.52411849048
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 579.72557647,
                "monoMz_": 579.72557647,
                "intensity_": 8207.690545963978
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 580.32557647,
                "monoMz_": 580.32557647,
                "intensity_": 6641.777228335853
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 580.5255764699999,
                "monoMz_": 580.5255764699999,
                "intensity_": 3203.1860425637715
            }
        ],
        "monoMass_": 2893.5915,
        "charge_": 5,
        "intensity_": 37288.42
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 580.52837647,
                "monoMz_": 580.52837647,
                "intensity_": 504.199839098643
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 581.52837647,
                "monoMz_": 581.52837647,
                "intensity_": 152.89016090135698
            }
        ],
        "monoMass_": 579.5211,
        "charge_": 1,
        "intensity_": 3716.05
    },
    {
        "displayColor_": "red",
        "displayLevel_": 4,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 581.5861139699999,
                "monoMz_": 581.5861139699999,
                "intensity_": 15692.507627664025
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 581.7111139699999,
                "monoMz_": 581.7111139699999,
                "intensity_": 14336.928559557116
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 581.4611139699999,
                "monoMz_": 581.4611139699999,
                "intensity_": 11946.740109682174
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 581.8361139699999,
                "monoMz_": 581.8361139699999,
                "intensity_": 10193.82941616911
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 581.9611139699999,
                "monoMz_": 581.9611139699999,
                "intensity_": 5990.28022259812
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 581.3361139699999,
                "monoMz_": 581.3361139699999,
                "intensity_": 4776.744555009309
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 582.0861139699999,
                "monoMz_": 582.0861139699999,
                "intensity_": 3019.492372335974
            }
        ],
        "monoMass_": 4642.6307,
        "charge_": 8,
        "intensity_": 63778.48
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 1,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 587.2057542477777,
                "monoMz_": 587.2057542477777,
                "intensity_": 258.0736870439532
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 587.3168653588888,
                "monoMz_": 587.3168653588888,
                "intensity_": 733.3361840765377
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 587.42797647,
                "monoMz_": 587.42797647,
                "intensity_": 1086.9334474203408
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 587.5390875811111,
                "monoMz_": 587.5390875811111,
                "intensity_": 1114.7638159234623
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 587.6501986922221,
                "monoMz_": 587.6501986922221,
                "intensity_": 886.3151890927999
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 587.7613098033332,
                "monoMz_": 587.7613098033332,
                "intensity_": 580.6851871654728
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 587.8724209144444,
                "monoMz_": 587.8724209144444,
                "intensity_": 325.6115525029328
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 587.9835320255555,
                "monoMz_": 587.9835320255555,
                "intensity_": 160.3309281204683
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 588.0946431366666,
                "monoMz_": 588.0946431366666,
                "intensity_": 70.61625229373769
            }
        ],
        "monoMass_": 5275.7863,
        "charge_": 9,
        "intensity_": 4367.13
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 588.69427647,
                "monoMz_": 588.69427647,
                "intensity_": 657.25
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 589.19427647,
                "monoMz_": 589.19427647,
                "intensity_": 413.2189442425
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 589.69427647,
                "monoMz_": 589.69427647,
                "intensity_": 148.09543463
            }
        ],
        "monoMass_": 1175.374,
        "charge_": 2,
        "intensity_": 1777.5
    },
    {
        "displayColor_": "red",
        "displayLevel_": 3,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 589.5406764699999,
                "monoMz_": 589.5406764699999,
                "intensity_": 5230.60609721796
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 589.4295653588888,
                "monoMz_": 589.4295653588888,
                "intensity_": 5083.452985130379
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 589.6517875811111,
                "monoMz_": 589.6517875811111,
                "intensity_": 4171.416731501078
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 589.3184542477777,
                "monoMz_": 589.3184542477777,
                "intensity_": 3417.7094075587984
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 589.7628986922222,
                "monoMz_": 589.7628986922222,
                "intensity_": 2740.8946085486714
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 589.8740098033333,
                "monoMz_": 589.8740098033333,
                "intensity_": 1541.1693182117374
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 589.2073431366666,
                "monoMz_": 589.2073431366666,
                "intensity_": 1198.1422997076716
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 589.9851209144443,
                "monoMz_": 589.9851209144443,
                "intensity_": 760.8854939137365
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 590.0962320255554,
                "monoMz_": 590.0962320255554,
                "intensity_": 335.9839027820399
            }
        ],
        "monoMass_": 5293.8006,
        "charge_": 9,
        "intensity_": 22062.14
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 4,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 593.7181889699999,
                "monoMz_": 593.7181889699999,
                "intensity_": 3628.645741896392
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 593.8431889699999,
                "monoMz_": 593.8431889699999,
                "intensity_": 3375.506552556406
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 593.5931889699999,
                "monoMz_": 593.5931889699999,
                "intensity_": 2711.7907784738495
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 593.9681889699999,
                "monoMz_": 593.9681889699999,
                "intensity_": 2442.860158302619
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 594.0931889699999,
                "monoMz_": 594.0931889699999,
                "intensity_": 1460.7760058675133
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 593.4681889699999,
                "monoMz_": 593.4681889699999,
                "intensity_": 1063.654258103608
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 594.2181889699999,
                "monoMz_": 594.2181889699999,
                "intensity_": 749.1704538771311
            }
        ],
        "monoMass_": 4739.6873,
        "charge_": 8,
        "intensity_": 14025.42
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 1,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 595.8478514699999,
                "monoMz_": 595.8478514699999,
                "intensity_": 777.756442319627
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 596.0978514699999,
                "monoMz_": 596.0978514699999,
                "intensity_": 989.6935576803731
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 596.3478514699999,
                "monoMz_": 596.3478514699999,
                "intensity_": 674.3086143367947
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 596.5978514699999,
                "monoMz_": 596.5978514699999,
                "intensity_": 323.8849759483943
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 596.8478514699999,
                "monoMz_": 596.8478514699999,
                "intensity_": 122.33979139732568
            }
        ],
        "monoMass_": 2379.3623,
        "charge_": 4,
        "intensity_": 2668.52
    },
    {
        "displayColor_": "red",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 598.1499431366667,
                "monoMz_": 598.1499431366667,
                "intensity_": 1145.6128689124487
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 598.3166098033333,
                "monoMz_": 598.3166098033333,
                "intensity_": 2212.5733841678225
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 598.48327647,
                "monoMz_": 598.48327647,
                "intensity_": 2287.5015023789047
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 598.6499431366667,
                "monoMz_": 598.6499431366667,
                "intensity_": 1666.6717875070829
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 598.8166098033333,
                "monoMz_": 598.8166098033333,
                "intensity_": 953.744654397774
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 598.98327647,
                "monoMz_": 598.98327647,
                "intensity_": 453.9677347806722
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 599.1499431366667,
                "monoMz_": 599.1499431366667,
                "intensity_": 186.18281253063225
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 599.3166098033333,
                "monoMz_": 599.3166098033333,
                "intensity_": 67.36849762109539
            }
        ],
        "monoMass_": 3582.856,
        "charge_": 6,
        "intensity_": 4875.34
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 598.44738897,
                "monoMz_": 598.44738897,
                "intensity_": 303.6628152586588
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 598.57238897,
                "monoMz_": 598.57238897,
                "intensity_": 781.93748812377
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 598.69738897,
                "monoMz_": 598.69738897,
                "intensity_": 1055.757296347167
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 598.82238897,
                "monoMz_": 598.82238897,
                "intensity_": 990.2622944806824
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 598.94738897,
                "monoMz_": 598.94738897,
                "intensity_": 722.2099933490731
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 599.07238897,
                "monoMz_": 599.07238897,
                "intensity_": 435.026578189681
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 599.19738897,
                "monoMz_": 599.19738897,
                "intensity_": 224.66255549972814
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 599.32238897,
                "monoMz_": 599.32238897,
                "intensity_": 102.018049255059
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 599.44738897,
                "monoMz_": 599.44738897,
                "intensity_": 41.47722017440224
            }
        ],
        "monoMass_": 4779.5209,
        "charge_": 8,
        "intensity_": 2597.79
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 600.1277653588888,
                "monoMz_": 600.1277653588888,
                "intensity_": 230.5132652736793
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 600.2388764699999,
                "monoMz_": 600.2388764699999,
                "intensity_": 670.4488887027658
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 600.349987581111,
                "monoMz_": 600.349987581111,
                "intensity_": 1025.9234734970332
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 600.4610986922222,
                "monoMz_": 600.4610986922222,
                "intensity_": 1093.2678084200363
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 600.5722098033333,
                "monoMz_": 600.5722098033333,
                "intensity_": 907.5378099983928
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 600.6833209144444,
                "monoMz_": 600.6833209144444,
                "intensity_": 623.093348553045
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 600.7944320255555,
                "monoMz_": 600.7944320255555,
                "intensity_": 367.1783665291029
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 600.9055431366666,
                "monoMz_": 600.9055431366666,
                "intensity_": 190.41566771978503
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 601.0166542477778,
                "monoMz_": 601.0166542477778,
                "intensity_": 88.47526657574122
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 601.1277653588888,
                "monoMz_": 601.1277653588888,
                "intensity_": 37.328590714084136
            }
        ],
        "monoMass_": 5392.0844,
        "charge_": 9,
        "intensity_": 2825.44
    },
    {
        "displayColor_": "red",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 600.81926397,
                "monoMz_": 600.81926397,
                "intensity_": 1724.4345948990363
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 600.94426397,
                "monoMz_": 600.94426397,
                "intensity_": 1623.2273905696445
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 600.69426397,
                "monoMz_": 600.69426397,
                "intensity_": 1272.2569250042225
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 601.06926397,
                "monoMz_": 601.06926397,
                "intensity_": 1187.7830649203142
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 601.19426397,
                "monoMz_": 601.19426397,
                "intensity_": 717.7177832395749
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 600.56926397,
                "monoMz_": 600.56926397,
                "intensity_": 491.96653223063845
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 601.31926397,
                "monoMz_": 601.31926397,
                "intensity_": 371.7654051009634
            }
        ],
        "monoMass_": 4796.4959,
        "charge_": 8,
        "intensity_": 7223.97
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 601.6519764699999,
                "monoMz_": 601.6519764699999,
                "intensity_": 3440.6905648744564
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 601.4853098033333,
                "monoMz_": 601.4853098033333,
                "intensity_": 3306.4098052126487
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 601.8186431366667,
                "monoMz_": 601.8186431366667,
                "intensity_": 2521.5782749552413
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 601.3186431366667,
                "monoMz_": 601.3186431366667,
                "intensity_": 1699.1421988657175
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 601.9853098033333,
                "monoMz_": 601.9853098033333,
                "intensity_": 1450.7594351255436
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 602.1519764699999,
                "monoMz_": 602.1519764699999,
                "intensity_": 694.0445744330831
            }
        ],
        "monoMass_": 3601.8682,
        "charge_": 6,
        "intensity_": 12944.01
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 602.10110147,
                "monoMz_": 602.10110147,
                "intensity_": 240.78875360144914
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 602.1844348033333,
                "monoMz_": 602.1844348033333,
                "intensity_": 937.2347356476928
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 602.2677681366666,
                "monoMz_": 602.2677681366666,
                "intensity_": 1887.913573540969
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 602.35110147,
                "monoMz_": 602.35110147,
                "intensity_": 2615.1562094261
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 602.4344348033333,
                "monoMz_": 602.4344348033333,
                "intensity_": 2794.404218093117
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 602.5177681366666,
                "monoMz_": 602.5177681366666,
                "intensity_": 2450.7775888761043
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 602.60110147,
                "monoMz_": 602.60110147,
                "intensity_": 1833.7200717850428
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 602.6844348033333,
                "monoMz_": 602.6844348033333,
                "intensity_": 1201.6911149349141
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 602.7677681366666,
                "monoMz_": 602.7677681366666,
                "intensity_": 702.9397024003749
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 602.85110147,
                "monoMz_": 602.85110147,
                "intensity_": 372.31113262765393
            },
            {
                "displayLevel_": -1,
                "peakId_": "10",
                "pos_": 602.9344348033333,
                "monoMz_": 602.9344348033333,
                "intensity_": 180.5357819068832
            }
        ],
        "monoMass_": 7213.1259,
        "charge_": 12,
        "intensity_": 8102.25
    },
    {
        "displayColor_": "red",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 602.3288964699999,
                "monoMz_": 602.3288964699999,
                "intensity_": 1342.4076651135542
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 602.52889647,
                "monoMz_": 602.52889647,
                "intensity_": 2170.5592309692315
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 602.72889647,
                "monoMz_": 602.72889647,
                "intensity_": 1911.6399707461198
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 602.9288964699999,
                "monoMz_": 602.9288964699999,
                "intensity_": 1199.0821879034788
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 603.12889647,
                "monoMz_": 603.12889647,
                "intensity_": 594.5580868612225
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 603.3288964699999,
                "monoMz_": 603.3288964699999,
                "intensity_": 246.18076903076832
            }
        ],
        "monoMass_": 3006.6081,
        "charge_": 5,
        "intensity_": 5000.77
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 602.54335147,
                "monoMz_": 602.54335147,
                "intensity_": 2185.6161378151082
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 602.79335147,
                "monoMz_": 602.79335147,
                "intensity_": 2806.166852549164
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 603.04335147,
                "monoMz_": 603.04335147,
                "intensity_": 1931.1979104805962
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 603.29335147,
                "monoMz_": 603.29335147,
                "intensity_": 937.6031474508358
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 603.54335147,
                "monoMz_": 603.54335147,
                "intensity_": 358.1637873266919
            }
        ],
        "monoMass_": 2406.1443,
        "charge_": 4,
        "intensity_": 7942.47
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 3,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 603.5108431366666,
                "monoMz_": 603.5108431366666,
                "intensity_": 9323.038755840582
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 603.6775098033334,
                "monoMz_": 603.6775098033334,
                "intensity_": 18146.594839408903
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 603.84417647,
                "monoMz_": 603.84417647,
                "intensity_": 18906.914275925545
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 604.0108431366666,
                "monoMz_": 604.0108431366666,
                "intensity_": 13882.31937295209
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 604.1775098033334,
                "monoMz_": 604.1775098033334,
                "intensity_": 8005.847922942535
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 604.34417647,
                "monoMz_": 604.34417647,
                "intensity_": 3840.5448587841934
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 604.5108431366666,
                "monoMz_": 604.5108431366666,
                "intensity_": 1587.6128354728971
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 604.6775098033334,
                "monoMz_": 604.6775098033334,
                "intensity_": 579.1030915327506
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 604.84417647,
                "monoMz_": 604.84417647,
                "intensity_": 189.51572407445283
            }
        ],
        "monoMass_": 3615.0214,
        "charge_": 6,
        "intensity_": 61517.4
    },
    {
        "displayColor_": "red",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 603.9174050414285,
                "monoMz_": 603.9174050414285,
                "intensity_": 0.4363776712444914
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 604.0602621842856,
                "monoMz_": 604.0602621842856,
                "intensity_": 0.9914065907636356
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 604.2031193271428,
                "monoMz_": 604.2031193271428,
                "intensity_": 1.1909932375046677
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 604.34597647,
                "monoMz_": 604.34597647,
                "intensity_": 1
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 604.4888336128571,
                "monoMz_": 604.4888336128571,
                "intensity_": 0.6557639493326227
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 604.6316907557142,
                "monoMz_": 604.6316907557142,
                "intensity_": 0.35633249258341715
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 604.7745478985714,
                "monoMz_": 604.7745478985714,
                "intensity_": 0.16640958879136256
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 604.9174050414285,
                "monoMz_": 604.9174050414285,
                "intensity_": 0.06845489696104141
            }
        ],
        "monoMass_": 4220.3709,
        "charge_": 7,
        "intensity_": 7823.08
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 604.0128681366666,
                "monoMz_": 604.0128681366666,
                "intensity_": 127.95622338939818
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 604.09620147,
                "monoMz_": 604.09620147,
                "intensity_": 499.5127852396208
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 604.1795348033334,
                "monoMz_": 604.1795348033334,
                "intensity_": 1009.1999293359638
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 604.2628681366666,
                "monoMz_": 604.2628681366666,
                "intensity_": 1402.1981904225488
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 604.34620147,
                "monoMz_": 604.34620147,
                "intensity_": 1502.9172252881176
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 604.4295348033334,
                "monoMz_": 604.4295348033334,
                "intensity_": 1322.2092234871689
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 604.5128681366666,
                "monoMz_": 604.5128681366666,
                "intensity_": 992.4209556451964
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 604.59620147,
                "monoMz_": 604.59620147,
                "intensity_": 652.436927405214
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 604.6795348033334,
                "monoMz_": 604.6795348033334,
                "intensity_": 382.8802916841561
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 604.7628681366666,
                "monoMz_": 604.7628681366666,
                "intensity_": 203.45472915154568
            },
            {
                "displayLevel_": -1,
                "peakId_": "10",
                "pos_": 604.84620147,
                "monoMz_": 604.84620147,
                "intensity_": 98.98277471188231
            }
        ],
        "monoMass_": 7236.0671,
        "charge_": 12,
        "intensity_": 10415.49
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 605.23217647,
                "monoMz_": 605.23217647,
                "intensity_": 98.78401180195858
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 605.3090995469231,
                "monoMz_": 605.3090995469231,
                "intensity_": 419.22336943757534
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 605.3860226238461,
                "monoMz_": 605.3860226238461,
                "intensity_": 917.4230574249548
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 605.4629457007692,
                "monoMz_": 605.4629457007692,
                "intensity_": 1376.591691764126
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 605.5398687776923,
                "monoMz_": 605.5398687776923,
                "intensity_": 1589.545594789157
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 605.6167918546154,
                "monoMz_": 605.6167918546154,
                "intensity_": 1503.497912438182
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 605.6937149315385,
                "monoMz_": 605.6937149315385,
                "intensity_": 1211.250449353539
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 605.7706380084616,
                "monoMz_": 605.7706380084616,
                "intensity_": 853.5007847578966
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 605.8475610853847,
                "monoMz_": 605.8475610853847,
                "intensity_": 536.2261965060491
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 605.9244841623076,
                "monoMz_": 605.9244841623076,
                "intensity_": 304.75440521084306
            },
            {
                "displayLevel_": -1,
                "peakId_": "10",
                "pos_": 606.0014072392307,
                "monoMz_": 606.0014072392307,
                "intensity_": 158.4518179981958
            }
        ],
        "monoMass_": 7854.9237,
        "charge_": 13,
        "intensity_": 7017.25
    },
    {
        "displayColor_": "red",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 605.5910514699999,
                "monoMz_": 605.5910514699999,
                "intensity_": 559.1374605159533
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 605.7160514699999,
                "monoMz_": 605.7160514699999,
                "intensity_": 1460.4445250449148
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 605.8410514699999,
                "monoMz_": 605.8410514699999,
                "intensity_": 1998.5928236843201
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 605.9660514699999,
                "monoMz_": 605.9660514699999,
                "intensity_": 1898.9150851396012
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 606.0910514699999,
                "monoMz_": 606.0910514699999,
                "intensity_": 1402.2505584590715
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 606.2160514699999,
                "monoMz_": 606.2160514699999,
                "intensity_": 854.9558396635724
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 606.3410514699999,
                "monoMz_": 606.3410514699999,
                "intensity_": 446.8061911140678
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 606.4660514699999,
                "monoMz_": 606.4660514699999,
                "intensity_": 205.28092507902517
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 606.5910514699999,
                "monoMz_": 606.5910514699999,
                "intensity_": 84.43251245751132
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 606.7160514699999,
                "monoMz_": 606.7160514699999,
                "intensity_": 31.507176315679725
            }
        ],
        "monoMass_": 4836.6702,
        "charge_": 8,
        "intensity_": 6999.3
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 1,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 606.13289647,
                "monoMz_": 606.13289647,
                "intensity_": 8012.029994703643
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 606.3328964699999,
                "monoMz_": 606.3328964699999,
                "intensity_": 7107.417727099139
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 605.93289647,
                "monoMz_": 605.93289647,
                "intensity_": 4920.709202900162
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 606.53289647,
                "monoMz_": 606.53289647,
                "intensity_": 4491.1051961480525
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 606.7328964699999,
                "monoMz_": 606.7328964699999,
                "intensity_": 2243.7380735809757
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 606.93289647,
                "monoMz_": 606.93289647,
                "intensity_": 936.2632312564826
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 607.13289647,
                "monoMz_": 607.13289647,
                "intensity_": 337.48128527150783
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 607.3328964699999,
                "monoMz_": 607.3328964699999,
                "intensity_": 107.47000529635689
            }
        ],
        "monoMass_": 3024.6281,
        "charge_": 5,
        "intensity_": 18950.96
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 606.4355582881818,
                "monoMz_": 606.4355582881818,
                "intensity_": 439.1419114406341
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 606.526467379091,
                "monoMz_": 606.526467379091,
                "intensity_": 1580.9840810845635
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 606.61737647,
                "monoMz_": 606.61737647,
                "intensity_": 2956.741041560272
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 606.7082855609091,
                "monoMz_": 606.7082855609091,
                "intensity_": 3814.0839037538526
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 606.7991946518182,
                "monoMz_": 606.7991946518182,
                "intensity_": 3804.490567355775
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 606.8901037427273,
                "monoMz_": 606.8901037427273,
                "intensity_": 3120.890325962342
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 606.9810128336363,
                "monoMz_": 606.9810128336363,
                "intensity_": 2187.6048577527276
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 607.0719219245454,
                "monoMz_": 607.0719219245454,
                "intensity_": 1344.7960962461475
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 607.1628310154546,
                "monoMz_": 607.1628310154546,
                "intensity_": 738.7057823673174
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 607.2537401063637,
                "monoMz_": 607.2537401063637,
                "intensity_": 367.7247159764047
            }
        ],
        "monoMass_": 6659.7111,
        "charge_": 11,
        "intensity_": 15474.4
    },
    {
        "displayColor_": "red",
        "displayLevel_": 6,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 607.9725764699999,
                "monoMz_": 607.9725764699999,
                "intensity_": 108812.14852522618
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 608.0975764699999,
                "monoMz_": 608.0975764699999,
                "intensity_": 103392.90418808484
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 607.8475764699999,
                "monoMz_": 607.8475764699999,
                "intensity_": 79506.6365313863
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 608.2225764699999,
                "monoMz_": 608.2225764699999,
                "intensity_": 76355.61189765488
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 608.3475764699999,
                "monoMz_": 608.3475764699999,
                "intensity_": 46557.283154415345
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 607.7225764699999,
                "monoMz_": 607.7225764699999,
                "intensity_": 30436.76552664605
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 608.4725764699999,
                "monoMz_": 608.4725764699999,
                "intensity_": 24332.716285247134
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 608.5975764699999,
                "monoMz_": 608.5975764699999,
                "intensity_": 11180.107678942106
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 608.7225764699999,
                "monoMz_": 608.7225764699999,
                "intensity_": 4598.676691411827
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 608.8475764699999,
                "monoMz_": 608.8475764699999,
                "intensity_": 1716.1514747738247
            }
        ],
        "monoMass_": 4853.7224,
        "charge_": 8,
        "intensity_": 443089.42
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 1,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 608.0142098033333,
                "monoMz_": 608.0142098033333,
                "intensity_": 979.5911686307311
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 608.0975431366667,
                "monoMz_": 608.0975431366667,
                "intensity_": 3849.2138757829753
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 608.1808764699999,
                "monoMz_": 608.1808764699999,
                "intensity_": 7824.321014676941
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 608.2642098033333,
                "monoMz_": 608.2642098033333,
                "intensity_": 10933.579127021105
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 608.3475431366667,
                "monoMz_": 608.3475431366667,
                "intensity_": 11782.562475817398
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 608.4308764699999,
                "monoMz_": 608.4308764699999,
                "intensity_": 10419.501448827452
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 608.5142098033333,
                "monoMz_": 608.5142098033333,
                "intensity_": 7859.451313826716
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 608.5975431366667,
                "monoMz_": 608.5975431366667,
                "intensity_": 5191.652237148371
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 608.6808764699999,
                "monoMz_": 608.6808764699999,
                "intensity_": 3060.7865355286485
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 608.7642098033333,
                "monoMz_": 608.7642098033333,
                "intensity_": 1633.7362969425346
            },
            {
                "displayLevel_": -1,
                "peakId_": "10",
                "pos_": 608.8475431366667,
                "monoMz_": 608.8475431366667,
                "intensity_": 798.3028949934334
            }
        ],
        "monoMass_": 7284.0832,
        "charge_": 12,
        "intensity_": 73984.46
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 3,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 610.97417647,
                "monoMz_": 610.97417647,
                "intensity_": 1466.0323600480424
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 611.09917647,
                "monoMz_": 611.09917647,
                "intensity_": 3845.580930436293
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 611.22417647,
                "monoMz_": 611.22417647,
                "intensity_": 5282.961991397908
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 611.34917647,
                "monoMz_": 611.34917647,
                "intensity_": 5037.368605274961
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 611.47417647,
                "monoMz_": 611.47417647,
                "intensity_": 3732.2249961126877
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 611.59917647,
                "monoMz_": 611.59917647,
                "intensity_": 2282.708390331013
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 611.72417647,
                "monoMz_": 611.72417647,
                "intensity_": 1196.538008602092
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 611.84917647,
                "monoMz_": 611.84917647,
                "intensity_": 551.3205625795608
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 611.97417647,
                "monoMz_": 611.97417647,
                "intensity_": 227.3898839891877
            }
        ],
        "monoMass_": 4879.7352,
        "charge_": 8,
        "intensity_": 19255.85
    },
    {
        "displayColor_": "red",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 613.1019139699999,
                "monoMz_": 613.1019139699999,
                "intensity_": 450.6771122072969
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 613.2269139699999,
                "monoMz_": 613.2269139699999,
                "intensity_": 1187.3300267597417
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 613.3519139699999,
                "monoMz_": 613.3519139699999,
                "intensity_": 1638.486904317492
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 613.4769139699999,
                "monoMz_": 613.4769139699999,
                "intensity_": 1569.5559077166763
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 613.6019139699999,
                "monoMz_": 613.6019139699999,
                "intensity_": 1168.4031763634705
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 613.7269139699999,
                "monoMz_": 613.7269139699999,
                "intensity_": 718.0730956825081
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 613.8519139699999,
                "monoMz_": 613.8519139699999,
                "intensity_": 378.250043329001
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 613.9769139699999,
                "monoMz_": 613.9769139699999,
                "intensity_": 175.15872314078865
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 614.1019139699999,
                "monoMz_": 614.1019139699999,
                "intensity_": 72.61336483931672
            }
        ],
        "monoMass_": 4896.7571,
        "charge_": 8,
        "intensity_": 5734.85
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 614.47849647,
                "monoMz_": 614.47849647,
                "intensity_": 616.0612492715911
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 614.6784964699999,
                "monoMz_": 614.6784964699999,
                "intensity_": 1018.8768850413645
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 614.87849647,
                "monoMz_": 614.87849647,
                "intensity_": 915.6669184933592
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 615.07849647,
                "monoMz_": 615.07849647,
                "intensity_": 585.2877439170082
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 615.2784964699999,
                "monoMz_": 615.2784964699999,
                "intensity_": 295.5044452290241
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 615.47849647,
                "monoMz_": 615.47849647,
                "intensity_": 124.53316464091503
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 615.6784964699999,
                "monoMz_": 615.6784964699999,
                "intensity_": 45.31397889115907
            }
        ],
        "monoMass_": 3067.3561,
        "charge_": 5,
        "intensity_": 4109.29
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 614.95383897,
                "monoMz_": 614.95383897,
                "intensity_": 2232.6447472622044
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 615.07883897,
                "monoMz_": 615.07883897,
                "intensity_": 2148.535371041627
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 614.82883897,
                "monoMz_": 614.82883897,
                "intensity_": 1609.8750411594576
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 615.20383897,
                "monoMz_": 615.20383897,
                "intensity_": 1606.2706594794774
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 615.32883897,
                "monoMz_": 615.32883897,
                "intensity_": 991.182456935476
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 614.70383897,
                "monoMz_": 614.70383897,
                "intensity_": 607.7252527377954
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 615.45383897,
                "monoMz_": 615.45383897,
                "intensity_": 524.1326667968664
            }
        ],
        "monoMass_": 4909.5725,
        "charge_": 8,
        "intensity_": 5312.34
    },
    {
        "displayColor_": "red",
        "displayLevel_": 4,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 620.35877647,
                "monoMz_": 620.35877647,
                "intensity_": 20524.48082887186
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 620.85877647,
                "monoMz_": 620.85877647,
                "intensity_": 13662.128668286075
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 621.35877647,
                "monoMz_": 621.35877647,
                "intensity_": 5153.829313786262
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 621.85877647,
                "monoMz_": 621.85877647,
                "intensity_": 1413.0191711281389
            }
        ],
        "monoMass_": 1238.703,
        "charge_": 2,
        "intensity_": 40636.09
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 620.67527647,
                "monoMz_": 620.67527647,
                "intensity_": 245.01121319620188
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 620.8181336128572,
                "monoMz_": 620.8181336128572,
                "intensity_": 572.0281317630628
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 620.9609907557143,
                "monoMz_": 620.9609907557143,
                "intensity_": 704.5716937795547
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 621.1038478985714,
                "monoMz_": 621.1038478985714,
                "intensity_": 605.57
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 621.2467050414286,
                "monoMz_": 621.2467050414286,
                "intensity_": 406.0239738545111
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 621.3895621842858,
                "monoMz_": 621.3895621842858,
                "intensity_": 225.38453022565676
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 621.5324193271429,
                "monoMz_": 621.5324193271429,
                "intensity_": 107.45594817325272
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 621.67527647,
                "monoMz_": 621.67527647,
                "intensity_": 45.10544683530298
            }
        ],
        "monoMass_": 4337.676,
        "charge_": 7,
        "intensity_": 3833.76
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 621.9741014699999,
                "monoMz_": 621.9741014699999,
                "intensity_": 472.2800228256128
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 622.0991014699999,
                "monoMz_": 622.0991014699999,
                "intensity_": 1261.7451677305132
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 622.2241014699999,
                "monoMz_": 622.2241014699999,
                "intensity_": 1764.3431336620345
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 622.3491014699999,
                "monoMz_": 622.3491014699999,
                "intensity_": 1711.6401752649456
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 622.4741014699999,
                "monoMz_": 622.4741014699999,
                "intensity_": 1289.852653541413
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 622.5991014699999,
                "monoMz_": 622.5991014699999,
                "intensity_": 802.2113419357091
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 622.7241014699999,
                "monoMz_": 622.7241014699999,
                "intensity_": 427.5305115538965
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 622.8491014699999,
                "monoMz_": 622.8491014699999,
                "intensity_": 200.2668663379655
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 622.9741014699999,
                "monoMz_": 622.9741014699999,
                "intensity_": 83.9708238461606
            }
        ],
        "monoMass_": 4967.7346,
        "charge_": 8,
        "intensity_": 4820.27
    },
    {
        "displayColor_": "red",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 622.67972647,
                "monoMz_": 622.67972647,
                "intensity_": 469.3786018604847
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 622.8463931366666,
                "monoMz_": 622.8463931366666,
                "intensity_": 941.319125577469
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 623.0130598033334,
                "monoMz_": 623.0130598033334,
                "intensity_": 1007.4645318015499
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 623.17972647,
                "monoMz_": 623.17972647,
                "intensity_": 758.3748804289996
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 623.3463931366666,
                "monoMz_": 623.3463931366666,
                "intensity_": 447.7811817064423
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 623.5130598033334,
                "monoMz_": 623.5130598033334,
                "intensity_": 219.73088565983366
            }
        ],
        "monoMass_": 3730.0347,
        "charge_": 6,
        "intensity_": 2656.4
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 623.9774139699999,
                "monoMz_": 623.9774139699999,
                "intensity_": 2905.518802640306
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 624.1024139699999,
                "monoMz_": 624.1024139699999,
                "intensity_": 7794.155112572602
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 624.2274139699999,
                "monoMz_": 624.2274139699999,
                "intensity_": 10939.28651554939
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 624.3524139699999,
                "monoMz_": 624.3524139699999,
                "intensity_": 10648.85644424071
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 624.4774139699999,
                "monoMz_": 624.4774139699999,
                "intensity_": 8050.432512594007
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 624.6024139699999,
                "monoMz_": 624.6024139699999,
                "intensity_": 5022.045394096892
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 624.7274139699999,
                "monoMz_": 624.7274139699999,
                "intensity_": 2684.1709820345113
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 624.8524139699999,
                "monoMz_": 624.8524139699999,
                "intensity_": 1260.8203757077333
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 624.9774139699999,
                "monoMz_": 624.9774139699999,
                "intensity_": 530.0680507898777
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 625.1024139699999,
                "monoMz_": 625.1024139699999,
                "intensity_": 202.1634844506105
            }
        ],
        "monoMass_": 4983.7611,
        "charge_": 8,
        "intensity_": 44365.45
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 625.18329647,
                "monoMz_": 625.18329647,
                "intensity_": 28.927692602213472
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 625.2499631366666,
                "monoMz_": 625.2499631366666,
                "intensity_": 146.43094370012358
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 625.3166298033333,
                "monoMz_": 625.3166298033333,
                "intensity_": 381.1389516125536
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 625.3832964699999,
                "monoMz_": 625.3832964699999,
                "intensity_": 678.6027936373328
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 625.4499631366666,
                "monoMz_": 625.4499631366666,
                "intensity_": 927.9351121350364
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 625.5166298033332,
                "monoMz_": 625.5166298033332,
                "intensity_": 1037.6635034075052
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 625.5832964699999,
                "monoMz_": 625.5832964699999,
                "intensity_": 986.928430267205
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 625.6499631366665,
                "monoMz_": 625.6499631366665,
                "intensity_": 820.0453879541604
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 625.7166298033333,
                "monoMz_": 625.7166298033333,
                "intensity_": 606.9161010502062
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 625.7832964699999,
                "monoMz_": 625.7832964699999,
                "intensity_": 405.9861777605075
            },
            {
                "displayLevel_": -1,
                "peakId_": "10",
                "pos_": 625.8499631366666,
                "monoMz_": 625.8499631366666,
                "intensity_": 248.27875798942395
            },
            {
                "displayLevel_": -1,
                "peakId_": "11",
                "pos_": 625.9166298033332,
                "monoMz_": 625.9166298033332,
                "intensity_": 140.07837810889785
            },
            {
                "displayLevel_": -1,
                "peakId_": "12",
                "pos_": 625.9832964699999,
                "monoMz_": 625.9832964699999,
                "intensity_": 73.4553485221445
            },
            {
                "displayLevel_": -1,
                "peakId_": "13",
                "pos_": 626.0499631366665,
                "monoMz_": 626.0499631366665,
                "intensity_": 36.018462386398326
            }
        ],
        "monoMass_": 9362.6403,
        "charge_": 15,
        "intensity_": 5404.92
    },
    {
        "displayColor_": "red",
        "displayLevel_": 1,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 625.2076764699999,
                "monoMz_": 625.2076764699999,
                "intensity_": 2625.5716325832504
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 625.4576764699999,
                "monoMz_": 625.4576764699999,
                "intensity_": 3525.427903613859
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 625.7076764699999,
                "monoMz_": 625.7076764699999,
                "intensity_": 2527.3159121613476
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 625.9576764699999,
                "monoMz_": 625.9576764699999,
                "intensity_": 1274.8685171527943
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 626.2076764699999,
                "monoMz_": 626.2076764699999,
                "intensity_": 505.0520963861407
            }
        ],
        "monoMass_": 2496.8016,
        "charge_": 4,
        "intensity_": 5750.4
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 625.5156598033333,
                "monoMz_": 625.5156598033333,
                "intensity_": 1035.3180342060964
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 625.6823264699999,
                "monoMz_": 625.6823264699999,
                "intensity_": 2087.6037169159904
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 625.8489931366666,
                "monoMz_": 625.8489931366666,
                "intensity_": 2244.8819657939034
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 626.0156598033333,
                "monoMz_": 626.0156598033333,
                "intensity_": 1697.0556972872548
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 626.1823264699999,
                "monoMz_": 626.1823264699999,
                "intensity_": 1005.9675494325205
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 626.3489931366666,
                "monoMz_": 626.3489931366666,
                "intensity_": 495.46307217414596
            }
        ],
        "monoMass_": 3747.0503,
        "charge_": 6,
        "intensity_": 6327.23
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 3,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 627.58427647,
                "monoMz_": 627.58427647,
                "intensity_": 13169.091311549388
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 627.70927647,
                "monoMz_": 627.70927647,
                "intensity_": 12878.415885120648
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 627.83427647,
                "monoMz_": 627.83427647,
                "intensity_": 9777.87317454728
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 627.45927647,
                "monoMz_": 627.45927647,
                "intensity_": 9336.40190747373
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 627.95927647,
                "monoMz_": 627.95927647,
                "intensity_": 6124.504389660901
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 627.33427647,
                "monoMz_": 627.33427647,
                "intensity_": 3461.4724736400476
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 628.08427647,
                "monoMz_": 628.08427647,
                "intensity_": 3286.1269446631595
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 628.20927647,
                "monoMz_": 628.20927647,
                "intensity_": 1549.3323990761705
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 628.33427647,
                "monoMz_": 628.33427647,
                "intensity_": 653.7086884506132
            }
        ],
        "monoMass_": 5010.616,
        "charge_": 8,
        "intensity_": 52177.22
    },
    {
        "displayColor_": "red",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 628.13719647,
                "monoMz_": 628.13719647,
                "intensity_": 1254.7413908059088
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 628.33719647,
                "monoMz_": 628.33719647,
                "intensity_": 2121.658609194091
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 628.53719647,
                "monoMz_": 628.53719647,
                "intensity_": 1945.0627579729564
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 628.73719647,
                "monoMz_": 628.73719647,
                "intensity_": 1266.5825796696818
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 628.93719647,
                "monoMz_": 628.93719647,
                "intensity_": 650.9694797812984
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 629.13719647,
                "monoMz_": 629.13719647,
                "intensity_": 279.1406825675608
            }
        ],
        "monoMass_": 3135.6496,
        "charge_": 5,
        "intensity_": 8882.84
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 631.0019931366667,
                "monoMz_": 631.0019931366667,
                "intensity_": 2207.765455850464
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 631.1686598033333,
                "monoMz_": 631.1686598033333,
                "intensity_": 4509.139324204373
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 631.3353264699999,
                "monoMz_": 631.3353264699999,
                "intensity_": 4907.9044504005
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 631.5019931366667,
                "monoMz_": 631.5019931366667,
                "intensity_": 3753.625985365309
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 631.6686598033333,
                "monoMz_": 631.6686598033333,
                "intensity_": 2250.409059722926
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 631.8353264699999,
                "monoMz_": 631.8353264699999,
                "intensity_": 1120.8181884170067
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 632.0019931366667,
                "monoMz_": 632.0019931366667,
                "intensity_": 480.5955495995
            }
        ],
        "monoMass_": 3779.9683,
        "charge_": 6,
        "intensity_": 16647.54
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 631.94113647,
                "monoMz_": 631.94113647,
                "intensity_": 3422.1508397544326
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 632.14113647,
                "monoMz_": 632.14113647,
                "intensity_": 3154.3269840561975
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 632.3411364699999,
                "monoMz_": 632.3411364699999,
                "intensity_": 2063.9027647142802
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 631.74113647,
                "monoMz_": 631.74113647,
                "intensity_": 2010.8491602455672
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 632.54113647,
                "monoMz_": 632.54113647,
                "intensity_": 1065.4347499293035
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 632.74113647,
                "monoMz_": 632.74113647,
                "intensity_": 458.75680725987286
            }
        ],
        "monoMass_": 3153.6693,
        "charge_": 5,
        "intensity_": 12432.03
    },
    {
        "displayColor_": "red",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 634.00417647,
                "monoMz_": 634.00417647,
                "intensity_": 100.99661809468738
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 635.00417647,
                "monoMz_": 635.00417647,
                "intensity_": 33.98752533642189
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 636.00417647,
                "monoMz_": 636.00417647,
                "intensity_": 7.208935596396414
            }
        ],
        "monoMass_": 632.9969,
        "charge_": 1,
        "intensity_": 2845.67
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 634.320133612857,
                "monoMz_": 634.320133612857,
                "intensity_": 651.2692344914682
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 634.4629907557141,
                "monoMz_": 634.4629907557141,
                "intensity_": 1558.8162122792244
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 634.6058478985714,
                "monoMz_": 634.6058478985714,
                "intensity_": 1964.516834805644
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 634.7487050414285,
                "monoMz_": 634.7487050414285,
                "intensity_": 1725.194836310564
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 634.8915621842856,
                "monoMz_": 634.8915621842856,
                "intensity_": 1180.6531651943562
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 635.0344193271428,
                "monoMz_": 635.0344193271428,
                "intensity_": 668.4361059169124
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 635.1772764699999,
                "monoMz_": 635.1772764699999,
                "intensity_": 324.8492230603472
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 635.320133612857,
                "monoMz_": 635.320133612857,
                "intensity_": 138.93377378439084
            }
        ],
        "monoMass_": 4433.19,
        "charge_": 7,
        "intensity_": 8841.98
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 5,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 634.3377598033333,
                "monoMz_": 634.3377598033333,
                "intensity_": 13768.98414343111
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 634.1710931366666,
                "monoMz_": 634.1710931366666,
                "intensity_": 12649.636229768976
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 634.50442647,
                "monoMz_": 634.50442647,
                "intensity_": 10531.16015679783
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 634.6710931366666,
                "monoMz_": 634.6710931366666,
                "intensity_": 6314.001614010932
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 634.00442647,
                "monoMz_": 634.00442647,
                "intensity_": 6193.167762965009
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 634.8377598033333,
                "monoMz_": 634.8377598033333,
                "intensity_": 3144.8175279209136
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 635.00442647,
                "monoMz_": 635.00442647,
                "intensity_": 1348.5158565688907
            }
        ],
        "monoMass_": 3797.9829,
        "charge_": 6,
        "intensity_": 46902.57
    },
    {
        "displayColor_": "red",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 636.7509193271428,
                "monoMz_": 636.7509193271428,
                "intensity_": 1547.0592367485795
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 636.8937764699999,
                "monoMz_": 636.8937764699999,
                "intensity_": 3703.661493068996
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 637.0366336128571,
                "monoMz_": 637.0366336128571,
                "intensity_": 4671.63283213971
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 637.1794907557143,
                "monoMz_": 637.1794907557143,
                "intensity_": 4108.043403396767
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 637.3223478985714,
                "monoMz_": 637.3223478985714,
                "intensity_": 2816.209389733819
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 637.4652050414286,
                "monoMz_": 637.4652050414286,
                "intensity_": 1597.6521794267426
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 637.6080621842857,
                "monoMz_": 637.6080621842857,
                "intensity_": 778.2160602825073
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 637.7509193271428,
                "monoMz_": 637.7509193271428,
                "intensity_": 333.67716786029064
            }
        ],
        "monoMass_": 4450.2055,
        "charge_": 7,
        "intensity_": 15880.08
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 1,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 638.27985647,
                "monoMz_": 638.27985647,
                "intensity_": 696.9718929287155
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 638.47985647,
                "monoMz_": 638.47985647,
                "intensity_": 1194.1013075777546
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 638.67985647,
                "monoMz_": 638.67985647,
                "intensity_": 1108.298692422245
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 638.8798564699999,
                "monoMz_": 638.8798564699999,
                "intensity_": 730.3148076222352
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 639.07985647,
                "monoMz_": 639.07985647,
                "intensity_": 379.74104069714974
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 639.27985647,
                "monoMz_": 639.27985647,
                "intensity_": 164.72718289734502
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 639.47985647,
                "monoMz_": 639.47985647,
                "intensity_": 61.684574880511676
            }
        ],
        "monoMass_": 3186.3629,
        "charge_": 5,
        "intensity_": 3439.79
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 3,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 642.84362647,
                "monoMz_": 642.84362647,
                "intensity_": 7441.077975528122
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 643.34362647,
                "monoMz_": 643.34362647,
                "intensity_": 5119.5197619823375
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 643.84362647,
                "monoMz_": 643.84362647,
                "intensity_": 1995.5220244718782
            }
        ],
        "monoMass_": 1283.6727,
        "charge_": 2,
        "intensity_": 15058.63
    },
    {
        "displayColor_": "red",
        "displayLevel_": 1,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 645.48751397,
                "monoMz_": 645.48751397,
                "intensity_": 331.1708126740215
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 645.61251397,
                "monoMz_": 645.61251397,
                "intensity_": 921.3711271606077
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 645.73751397,
                "monoMz_": 645.73751397,
                "intensity_": 1338.2091762088116
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 645.86251397,
                "monoMz_": 645.86251397,
                "intensity_": 1345.7788728393923
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 645.98751397,
                "monoMz_": 645.98751397,
                "intensity_": 1049.6921789355756
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 646.11251397,
                "monoMz_": 646.11251397,
                "intensity_": 674.9290988793716
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 646.23751397,
                "monoMz_": 646.23751397,
                "intensity_": 371.5187570962953
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 646.36251397,
                "monoMz_": 646.36251397,
                "intensity_": 179.6175398438607
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 646.48751397,
                "monoMz_": 646.48751397,
                "intensity_": 77.68652541804786
            }
        ],
        "monoMass_": 5155.8419,
        "charge_": 8,
        "intensity_": 4413.03
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 647.16475647,
                "monoMz_": 647.16475647,
                "intensity_": 1822.0928104414427
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 647.36475647,
                "monoMz_": 647.36475647,
                "intensity_": 3168.4375054083143
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 647.56475647,
                "monoMz_": 647.56475647,
                "intensity_": 2977.805864776669
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 647.76475647,
                "monoMz_": 647.76475647,
                "intensity_": 1984.1614622274583
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 647.96475647,
                "monoMz_": 647.96475647,
                "intensity_": 1042.2853362853627
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 648.16475647,
                "monoMz_": 648.16475647,
                "intensity_": 456.4824945916856
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 648.36475647,
                "monoMz_": 648.36475647,
                "intensity_": 172.503425444452
            }
        ],
        "monoMass_": 3230.7874,
        "charge_": 5,
        "intensity_": 8122.1
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 4,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 648.47233897,
                "monoMz_": 648.47233897,
                "intensity_": 2401.0075456647783
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 648.59733897,
                "monoMz_": 648.59733897,
                "intensity_": 6706.236481351715
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 648.72233897,
                "monoMz_": 648.72233897,
                "intensity_": 9775.110760837044
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 648.84733897,
                "monoMz_": 648.84733897,
                "intensity_": 9863.030604689695
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 648.97233897,
                "monoMz_": 648.97233897,
                "intensity_": 7716.989304277569
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 649.09733897,
                "monoMz_": 649.09733897,
                "intensity_": 4976.4744479017345
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 649.22233897,
                "monoMz_": 649.22233897,
                "intensity_": 2747.028401787171
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 649.34733897,
                "monoMz_": 649.34733897,
                "intensity_": 1331.6863702930752
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 649.47233897,
                "monoMz_": 649.47233897,
                "intensity_": 577.4693953103043
            }
        ],
        "monoMass_": 5179.7205,
        "charge_": 8,
        "intensity_": 34087.68
    },
    {
        "displayColor_": "red",
        "displayLevel_": 1,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 649.26567647,
                "monoMz_": 649.26567647,
                "intensity_": 1415.2998654446244
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 649.76567647,
                "monoMz_": 649.76567647,
                "intensity_": 973.9001345553755
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 650.26567647,
                "monoMz_": 650.26567647,
                "intensity_": 379.66214289076436
            }
        ],
        "monoMass_": 1296.5168,
        "charge_": 2,
        "intensity_": 2395.92
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 652.3293764699999,
                "monoMz_": 652.3293764699999,
                "intensity_": 3101.6
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 653.3293764699999,
                "monoMz_": 653.3293764699999,
                "intensity_": 1044.111551288
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 654.3293764699999,
                "monoMz_": 654.3293764699999,
                "intensity_": 221.50603672
            }
        ],
        "monoMass_": 651.3221,
        "charge_": 1,
        "intensity_": 4147.42
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 654.94643647,
                "monoMz_": 654.94643647,
                "intensity_": 4194.085770551533
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 655.1464364699999,
                "monoMz_": 655.1464364699999,
                "intensity_": 3989.497889196309
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 655.34643647,
                "monoMz_": 655.34643647,
                "intensity_": 2689.2861719624434
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 654.7464364699999,
                "monoMz_": 654.7464364699999,
                "intensity_": 2381.301661869244
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 655.54643647,
                "monoMz_": 655.54643647,
                "intensity_": 1428.8455021606974
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 655.7464364699999,
                "monoMz_": 655.7464364699999,
                "intensity_": 632.8842294484675
            }
        ],
        "monoMass_": 3268.6958,
        "charge_": 5,
        "intensity_": 13676.48
    },
    {
        "displayColor_": "red",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 657.20657647,
                "monoMz_": 657.20657647,
                "intensity_": 1204.4529962407964
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 657.3732431366666,
                "monoMz_": 657.3732431366666,
                "intensity_": 2562.405366502833
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 657.5399098033333,
                "monoMz_": 657.5399098033333,
                "intensity_": 2895.523502809997
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 657.70657647,
                "monoMz_": 657.70657647,
                "intensity_": 2293.976497190003
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 657.8732431366666,
                "monoMz_": 657.8732431366666,
                "intensity_": 1422.473291519917
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 658.0399098033333,
                "monoMz_": 658.0399098033333,
                "intensity_": 731.9940456916678
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 658.20657647,
                "monoMz_": 658.20657647,
                "intensity_": 324.06311043300104
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 658.3732431366666,
                "monoMz_": 658.3732431366666,
                "intensity_": 126.45741405809217
            }
        ],
        "monoMass_": 3937.1958,
        "charge_": 6,
        "intensity_": 10421.07
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 4,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 657.5603764699999,
                "monoMz_": 657.5603764699999,
                "intensity_": 5200.759990124634
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 658.0603764699999,
                "monoMz_": 658.0603764699999,
                "intensity_": 3654.6097222741128
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 658.5603764699999,
                "monoMz_": 658.5603764699999,
                "intensity_": 1447.5400098753657
            }
        ],
        "monoMass_": 1313.1062,
        "charge_": 2,
        "intensity_": 8990.46
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 657.89380147,
                "monoMz_": 657.89380147,
                "intensity_": 4207.563356422314
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 658.14380147,
                "monoMz_": 658.14380147,
                "intensity_": 3145.1249133435895
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 657.64380147,
                "monoMz_": 657.64380147,
                "intensity_": 2993.579967893221
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 658.39380147,
                "monoMz_": 658.39380147,
                "intensity_": 1649.9366435776851
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 658.64380147,
                "monoMz_": 658.64380147,
                "intensity_": 678.4727889712491
            }
        ],
        "monoMass_": 2626.5461,
        "charge_": 4,
        "intensity_": 11396.92
    },
    {
        "displayColor_": "red",
        "displayLevel_": 1,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 666.8761664699999,
                "monoMz_": 666.8761664699999,
                "intensity_": 392.0230270322743
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 666.97616647,
                "monoMz_": 666.97616647,
                "intensity_": 1411.3482430391869
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 667.07616647,
                "monoMz_": 667.07616647,
                "intensity_": 2639.4897482239357
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 667.17616647,
                "monoMz_": 667.17616647,
                "intensity_": 3404.841756960814
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 667.27616647,
                "monoMz_": 667.27616647,
                "intensity_": 3396.2777627800356
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 667.3761664699999,
                "monoMz_": 667.3761664699999,
                "intensity_": 2786.0262041622354
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 667.47616647,
                "monoMz_": 667.47616647,
                "intensity_": 1952.8800507183348
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 667.57616647,
                "monoMz_": 667.57616647,
                "intensity_": 1200.502668174202
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 667.67616647,
                "monoMz_": 667.67616647,
                "intensity_": 659.444405886612
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 667.77616647,
                "monoMz_": 667.77616647,
                "intensity_": 328.26872706988576
            }
        ],
        "monoMass_": 6658.6889,
        "charge_": 10,
        "intensity_": 15735.78
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 3,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 668.55197647,
                "monoMz_": 668.55197647,
                "intensity_": 2212.864619207753
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 668.7519764699999,
                "monoMz_": 668.7519764699999,
                "intensity_": 3979.42552355711
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 668.95197647,
                "monoMz_": 668.95197647,
                "intensity_": 3857.395723580637
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 669.15197647,
                "monoMz_": 669.15197647,
                "intensity_": 2646.5921953783063
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 669.35197647,
                "monoMz_": 669.35197647,
                "intensity_": 1430.1810199052
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 669.55197647,
                "monoMz_": 669.55197647,
                "intensity_": 644.011682609543
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 669.7519764699999,
                "monoMz_": 669.7519764699999,
                "intensity_": 250.16447644289025
            }
        ],
        "monoMass_": 3337.7235,
        "charge_": 5,
        "intensity_": 12223.63
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 669.3299320255555,
                "monoMz_": 669.3299320255555,
                "intensity_": 15.569214793606307
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 669.4410431366666,
                "monoMz_": 669.4410431366666,
                "intensity_": 50.58532480329507
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 669.5521542477778,
                "monoMz_": 669.5521542477778,
                "intensity_": 85.87831739045703
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 669.6632653588889,
                "monoMz_": 669.6632653588889,
                "intensity_": 101.01737997664678
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 669.77437647,
                "monoMz_": 669.77437647,
                "intensity_": 92.20922811966345
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 669.8854875811111,
                "monoMz_": 669.8854875811111,
                "intensity_": 69.41332491886502
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 669.9965986922222,
                "monoMz_": 669.9965986922222,
                "intensity_": 44.74929922876805
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 670.1077098033334,
                "monoMz_": 670.1077098033334,
                "intensity_": 25.345416202905838
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 670.2188209144444,
                "monoMz_": 670.2188209144444,
                "intensity_": 12.845676120491731
            }
        ],
        "monoMass_": 6014.9039,
        "charge_": 9,
        "intensity_": 3580.85
    },
    {
        "displayColor_": "red",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 670.6888431366666,
                "monoMz_": 670.6888431366666,
                "intensity_": 2677.5544328166484
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 671.02217647,
                "monoMz_": 671.02217647,
                "intensity_": 2881.8105043852065
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 671.3555098033332,
                "monoMz_": 671.3555098033332,
                "intensity_": 1679.5116692483869
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 671.6888431366666,
                "monoMz_": 671.6888431366666,
                "intensity_": 694.9271336442349
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 672.02217647,
                "monoMz_": 672.02217647,
                "intensity_": 227.23012427246258
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 672.3555098033332,
                "monoMz_": 672.3555098033332,
                "intensity_": 62.149874608302504
            }
        ],
        "monoMass_": 2009.0447,
        "charge_": 3,
        "intensity_": 7131.53
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 4,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 672.35273647,
                "monoMz_": 672.35273647,
                "intensity_": 10508.339621006491
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 672.55273647,
                "monoMz_": 672.55273647,
                "intensity_": 10256.417520768708
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 672.75273647,
                "monoMz_": 672.75273647,
                "intensity_": 7080.336811858355
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 672.1527364699999,
                "monoMz_": 672.1527364699999,
                "intensity_": 5796.431239632459
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 672.95273647,
                "monoMz_": 672.95273647,
                "intensity_": 3847.772621567731
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 673.1527364699999,
                "monoMz_": 673.1527364699999,
                "intensity_": 1741.860378993508
            }
        ],
        "monoMass_": 3355.7273,
        "charge_": 5,
        "intensity_": 36198.58
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 676.6911098033333,
                "monoMz_": 676.6911098033333,
                "intensity_": 3896.6839821493686
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 677.0244431366666,
                "monoMz_": 677.0244431366666,
                "intensity_": 4236.53440879224
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 677.35777647,
                "monoMz_": 677.35777647,
                "intensity_": 2490.0655912077596
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 677.6911098033333,
                "monoMz_": 677.6911098033333,
                "intensity_": 1038.0599361844374
            }
        ],
        "monoMass_": 2027.0515,
        "charge_": 3,
        "intensity_": 11757.51
    },
    {
        "displayColor_": "red",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 683.3733764699999,
                "monoMz_": 683.3733764699999,
                "intensity_": 3436.7435115295375
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 683.8733764699999,
                "monoMz_": 683.8733764699999,
                "intensity_": 2492.2552196030965
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 684.3733764699999,
                "monoMz_": 684.3733764699999,
                "intensity_": 1018.3564884704625
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 684.8733764699999,
                "monoMz_": 684.8733764699999,
                "intensity_": 301.6506419449782
            }
        ],
        "monoMass_": 1364.7322,
        "charge_": 2,
        "intensity_": 7692.52
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 683.6428907557142,
                "monoMz_": 683.6428907557142,
                "intensity_": 456.5580162676018
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 683.7857478985715,
                "monoMz_": 683.7857478985715,
                "intensity_": 1175.645520242473
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 683.9286050414286,
                "monoMz_": 683.9286050414286,
                "intensity_": 1587.334479757527
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 684.0714621842857,
                "monoMz_": 684.0714621842857,
                "intensity_": 1488.86253447791
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 684.2143193271429,
                "monoMz_": 684.2143193271429,
                "intensity_": 1085.845040365668
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 684.35717647,
                "monoMz_": 684.35717647,
                "intensity_": 654.0638549793598
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 684.5000336128571,
                "monoMz_": 684.5000336128571,
                "intensity_": 337.7808724495816
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 684.6428907557142,
                "monoMz_": 684.6428907557142,
                "intensity_": 153.3844641192106
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 684.7857478985715,
                "monoMz_": 684.7857478985715,
                "intensity_": 62.36113350589003
            }
        ],
        "monoMass_": 4778.4493,
        "charge_": 7,
        "intensity_": 4984.11
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 3,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 686.5059336128571,
                "monoMz_": 686.5059336128571,
                "intensity_": 3812.4432153063754
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 686.6487907557142,
                "monoMz_": 686.6487907557142,
                "intensity_": 3588.6906180045867
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 686.36307647,
                "monoMz_": 686.36307647,
                "intensity_": 2812.752247146193
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 686.7916478985715,
                "monoMz_": 686.7916478985715,
                "intensity_": 2625.994340699476
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 686.9345050414286,
                "monoMz_": 686.9345050414286,
                "intensity_": 1586.7567846936251
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 686.2202193271428,
                "monoMz_": 686.2202193271428,
                "intensity_": 1087.657643559579
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 687.0773621842857,
                "monoMz_": 687.0773621842857,
                "intensity_": 821.912585467341
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 687.2202193271428,
                "monoMz_": 687.2202193271428,
                "intensity_": 374.2996130940675
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 687.36307647,
                "monoMz_": 687.36307647,
                "intensity_": 152.60077895236475
            }
        ],
        "monoMass_": 4796.4906,
        "charge_": 7,
        "intensity_": 14460.74
    },
    {
        "displayColor_": "red",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 689.3429431366667,
                "monoMz_": 689.3429431366667,
                "intensity_": 2728.4337700076094
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 689.67627647,
                "monoMz_": 689.67627647,
                "intensity_": 3007.8534275066677
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 690.0096098033333,
                "monoMz_": 690.0096098033333,
                "intensity_": 1794.3465724933328
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 690.3429431366667,
                "monoMz_": 690.3429431366667,
                "intensity_": 759.6639729861668
            }
        ],
        "monoMass_": 2065.007,
        "charge_": 3,
        "intensity_": 7993.43
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 691.6623864699999,
                "monoMz_": 691.6623864699999,
                "intensity_": 211.65240765521867
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 691.7623864699999,
                "monoMz_": 691.7623864699999,
                "intensity_": 790.1405471797673
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 691.8623864699999,
                "monoMz_": 691.8623864699999,
                "intensity_": 1529.4599732179506
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 691.96238647,
                "monoMz_": 691.96238647,
                "intensity_": 2039.016636400072
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 692.06238647,
                "monoMz_": 692.06238647,
                "intensity_": 2099.5127742984337
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 692.1623864699999,
                "monoMz_": 692.1623864699999,
                "intensity_": 1776.1404840385421
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 692.2623864699999,
                "monoMz_": 692.2623864699999,
                "intensity_": 1282.9445890770571
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 692.3623864699999,
                "monoMz_": 692.3623864699999,
                "intensity_": 812.1968736545207
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 692.46238647,
                "monoMz_": 692.46238647,
                "intensity_": 459.21817803709337
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 692.56238647,
                "monoMz_": 692.56238647,
                "intensity_": 235.1972257015665
            }
        ],
        "monoMass_": 6906.5511,
        "charge_": 10,
        "intensity_": 5699.91
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 4,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 692.37872647,
                "monoMz_": 692.37872647,
                "intensity_": 12073.193986016873
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 692.87872647,
                "monoMz_": 692.87872647,
                "intensity_": 8887.199187242699
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 693.37872647,
                "monoMz_": 693.37872647,
                "intensity_": 3673.1770310435804
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 693.87872647,
                "monoMz_": 693.87872647,
                "intensity_": 1098.806013983127
            }
        ],
        "monoMass_": 1382.7429,
        "charge_": 2,
        "intensity_": 27332
    },
    {
        "displayColor_": "red",
        "displayLevel_": 6,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 694.68187647,
                "monoMz_": 694.68187647,
                "intensity_": 14912.216436397595
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 694.8247336128571,
                "monoMz_": 694.8247336128571,
                "intensity_": 14169.533330030672
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 694.5390193271428,
                "monoMz_": 694.5390193271428,
                "intensity_": 10896.027586581136
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 694.9675907557142,
                "monoMz_": 694.9675907557142,
                "intensity_": 10464.193807251524
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 695.1104478985715,
                "monoMz_": 695.1104478985715,
                "intensity_": 6380.4666344616335
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 694.3961621842857,
                "monoMz_": 694.3961621842857,
                "intensity_": 4171.222067653658
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 695.2533050414286,
                "monoMz_": 695.2533050414286,
                "intensity_": 3334.689523632926
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 695.3961621842857,
                "monoMz_": 695.3961621842857,
                "intensity_": 1532.1835636024068
            }
        ],
        "monoMass_": 4853.7222,
        "charge_": 7,
        "intensity_": 60366.41
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 697.9696478985715,
                "monoMz_": 697.9696478985715,
                "intensity_": 712.6055342643276
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 698.1125050414286,
                "monoMz_": 698.1125050414286,
                "intensity_": 1869.250862511972
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 698.2553621842857,
                "monoMz_": 698.2553621842857,
                "intensity_": 2567.9296412358
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 698.3982193271429,
                "monoMz_": 698.3982193271429,
                "intensity_": 2448.5521903013655
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 698.54107647,
                "monoMz_": 698.54107647,
                "intensity_": 1814.1510786722363
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 698.6839336128571,
                "monoMz_": 698.6839336128571,
                "intensity_": 1109.5734830902288
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 698.8267907557142,
                "monoMz_": 698.8267907557142,
                "intensity_": 581.6103587642
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 698.9696478985715,
                "monoMz_": 698.9696478985715,
                "intensity_": 267.98459212390304
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 699.1125050414286,
                "monoMz_": 699.1125050414286,
                "intensity_": 110.52913576962825
            }
        ],
        "monoMass_": 4878.7366,
        "charge_": 7,
        "intensity_": 7051.96
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 3,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 705.38667647,
                "monoMz_": 705.38667647,
                "intensity_": 2395.235464508786
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 705.88667647,
                "monoMz_": 705.88667647,
                "intensity_": 1798.3645354912142
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 706.38667647,
                "monoMz_": 706.38667647,
                "intensity_": 754.751485424099
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 706.88667647,
                "monoMz_": 706.88667647,
                "intensity_": 228.78254415267256
            }
        ],
        "monoMass_": 1408.7588,
        "charge_": 2,
        "intensity_": 5002.41
    },
    {
        "displayColor_": "red",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 710.3965320255555,
                "monoMz_": 710.3965320255555,
                "intensity_": 1472.0426075175792
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 710.5076431366666,
                "monoMz_": 710.5076431366666,
                "intensity_": 1412.5647995629618
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 710.2854209144444,
                "monoMz_": 710.2854209144444,
                "intensity_": 1187.8573924824207
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 710.6187542477778,
                "monoMz_": 710.6187542477778,
                "intensity_": 1115.9654223262514
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 710.7298653588889,
                "monoMz_": 710.7298653588889,
                "intensity_": 754.0204103336893
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 710.1743098033334,
                "monoMz_": 710.1743098033334,
                "intensity_": 662.312906627074
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 710.84097647,
                "monoMz_": 710.84097647,
                "intensity_": 447.1133344259325
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 710.9520875811111,
                "monoMz_": 710.9520875811111,
                "intensity_": 237.04037141186225
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 710.0631986922223,
                "monoMz_": 710.0631986922223,
                "intensity_": 192.2657519134866
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 711.0631986922223,
                "monoMz_": 711.0631986922223,
                "intensity_": 113.93374255368859
            }
        ],
        "monoMass_": 6381.5033,
        "charge_": 9,
        "intensity_": 4434.18
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 5,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 713.4017193271428,
                "monoMz_": 713.4017193271428,
                "intensity_": 10190.474928091407
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 713.5445764699999,
                "monoMz_": 713.5445764699999,
                "intensity_": 9919.925257797278
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 713.687433612857,
                "monoMz_": 713.687433612857,
                "intensity_": 7499.367583367576
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 713.2588621842856,
                "monoMz_": 713.2588621842856,
                "intensity_": 7260.632779608447
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 713.8302907557143,
                "monoMz_": 713.8302907557143,
                "intensity_": 4678.278387126706
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 713.1160050414285,
                "monoMz_": 713.1160050414285,
                "intensity_": 2706.6314123245365
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 713.9731478985714,
                "monoMz_": 713.9731478985714,
                "intensity_": 2500.4352026298
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 714.1160050414285,
                "monoMz_": 714.1160050414285,
                "intensity_": 1174.5152125975903
            },
            {
                "displayLevel_": -1,
                "peakId_": "8",
                "pos_": 714.2588621842856,
                "monoMz_": 714.2588621842856,
                "intensity_": 493.7840483543868
            },
            {
                "displayLevel_": -1,
                "peakId_": "9",
                "pos_": 714.4017193271428,
                "monoMz_": 714.4017193271428,
                "intensity_": 188.3250719085933
            }
        ],
        "monoMass_": 4984.7611,
        "charge_": 7,
        "intensity_": 34729.58
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 714.03567647,
                "monoMz_": 714.03567647,
                "intensity_": 2840.4
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 715.03567647,
                "monoMz_": 715.03567647,
                "intensity_": 1061.110638756
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 716.03567647,
                "monoMz_": 716.03567647,
                "intensity_": 245.43220384800003
            }
        ],
        "monoMass_": 713.0284,
        "charge_": 1,
        "intensity_": 3823.42
    },
    {
        "displayColor_": "red",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 714.1866598033332,
                "monoMz_": 714.1866598033332,
                "intensity_": 922.888947857368
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 714.35332647,
                "monoMz_": 714.35332647,
                "intensity_": 2130.4531789052185
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 714.5199931366666,
                "monoMz_": 714.5199931366666,
                "intensity_": 2595.920459913414
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 714.6866598033332,
                "monoMz_": 714.6866598033332,
                "intensity_": 2207.9977101484096
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 714.85332647,
                "monoMz_": 714.85332647,
                "intensity_": 1465.4178150663924
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 715.0199931366666,
                "monoMz_": 715.0199931366666,
                "intensity_": 805.3422947358907
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 715.1866598033332,
                "monoMz_": 715.1866598033332,
                "intensity_": 380.169540086586
            },
            {
                "displayLevel_": -1,
                "peakId_": "7",
                "pos_": 715.35332647,
                "monoMz_": 715.35332647,
                "intensity_": 158.01139398492478
            }
        ],
        "monoMass_": 4279.0763,
        "charge_": 6,
        "intensity_": 5001.97
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 4,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 719.3707098033333,
                "monoMz_": 719.3707098033333,
                "intensity_": 6845.427252127065
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 719.7040431366665,
                "monoMz_": 719.7040431366665,
                "intensity_": 7874.96506198708
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 720.0373764699999,
                "monoMz_": 720.0373764699999,
                "intensity_": 4884.310449053661
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 720.3707098033333,
                "monoMz_": 720.3707098033333,
                "intensity_": 2145.078585659893
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 720.7040431366665,
                "monoMz_": 720.7040431366665,
                "intensity_": 743.1349380129199
            }
        ],
        "monoMass_": 2155.0903,
        "charge_": 3,
        "intensity_": 18055.96
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 2,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 722.6626014699999,
                "monoMz_": 722.6626014699999,
                "intensity_": 3802.1622732108667
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 722.9126014699999,
                "monoMz_": 722.9126014699999,
                "intensity_": 3242.3156268487687
            },
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 722.4126014699999,
                "monoMz_": 722.4126014699999,
                "intensity_": 2440.7771213735905
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 723.1626014699999,
                "monoMz_": 723.1626014699999,
                "intensity_": 1974.7377267891334
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 723.4126014699999,
                "monoMz_": 723.4126014699999,
                "intensity_": 952.3261777710846
            }
        ],
        "monoMass_": 2885.6213,
        "charge_": 4,
        "intensity_": 10903.77
    },
    {
        "displayColor_": "red",
        "displayLevel_": 3,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 760.80783647,
                "monoMz_": 760.80783647,
                "intensity_": 591.0642108595653
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 761.00783647,
                "monoMz_": 761.00783647,
                "intensity_": 1207.2573426025579
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 761.20783647,
                "monoMz_": 761.20783647,
                "intensity_": 1314.0857891404348
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 761.40783647,
                "monoMz_": 761.40783647,
                "intensity_": 1005.073995368947
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 761.60783647,
                "monoMz_": 761.60783647,
                "intensity_": 602.5963649279039
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 761.80783647,
                "monoMz_": 761.80783647,
                "intensity_": 300.13543336471787
            }
        ],
        "monoMass_": 3799.0028,
        "charge_": 5,
        "intensity_": 2719.53
    },
    {
        "displayColor_": "darkorange",
        "displayLevel_": 5,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 768.1243264699999,
                "monoMz_": 768.1243264699999,
                "intensity_": 1800.5907889525758
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 768.3743264699999,
                "monoMz_": 768.3743264699999,
                "intensity_": 2977.9187320275632
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 768.6243264699999,
                "monoMz_": 768.6243264699999,
                "intensity_": 2676.2621754527563
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 768.8743264699999,
                "monoMz_": 768.8743264699999,
                "intensity_": 1710.6476374384035
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 769.1243264699999,
                "monoMz_": 769.1243264699999,
                "intensity_": 863.6845489032742
            },
            {
                "displayLevel_": -1,
                "peakId_": "5",
                "pos_": 769.3743264699999,
                "monoMz_": 769.3743264699999,
                "intensity_": 363.9788567073025
            },
            {
                "displayLevel_": -1,
                "peakId_": "6",
                "pos_": 769.6243264699999,
                "monoMz_": 769.6243264699999,
                "intensity_": 132.44126797243595
            }
        ],
        "monoMass_": 3068.4682,
        "charge_": 4,
        "intensity_": 7721.96
    },
    {
        "displayColor_": "blue",
        "displayLevel_": 4,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 863.9684264699999,
                "monoMz_": 863.9684264699999,
                "intensity_": 1956.033976669988
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 864.4684264699999,
                "monoMz_": 864.4684264699999,
                "intensity_": 1801.3047256344353
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 864.9684264699999,
                "monoMz_": 864.9684264699999,
                "intensity_": 908.8660233300116
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 865.4684264699999,
                "monoMz_": 865.4684264699999,
                "intensity_": 327.8021495836377
            }
        ],
        "monoMass_": 1725.9223,
        "charge_": 2,
        "intensity_": 4603.61
    },
    {
        "displayColor_": "red",
        "displayLevel_": 0,
        "peaks_": [
            {
                "displayLevel_": -1,
                "peakId_": "0",
                "pos_": 876.4128431366667,
                "monoMz_": 876.4128431366667,
                "intensity_": 1083.6951194613075
            },
            {
                "displayLevel_": -1,
                "peakId_": "1",
                "pos_": 876.7461764699999,
                "monoMz_": 876.7461764699999,
                "intensity_": 1523.1648805386922
            },
            {
                "displayLevel_": -1,
                "peakId_": "2",
                "pos_": 877.0795098033333,
                "monoMz_": 877.0795098033333,
                "intensity_": 1138.5553602181872
            },
            {
                "displayLevel_": -1,
                "peakId_": "3",
                "pos_": 877.4128431366667,
                "monoMz_": 877.4128431366667,
                "intensity_": 597.2876312784326
            },
            {
                "displayLevel_": -1,
                "peakId_": "4",
                "pos_": 877.7461764699999,
                "monoMz_": 877.7461764699999,
                "intensity_": 245.61149459217333
            }
        ],
        "monoMass_": 2626.2167,
        "charge_": 3,
        "intensity_": 4490.93
    }
];
const sampleMatchedPeakList = [
    {
        "ion": "B42",
        "ionPos": "42",
        "position": 42,
        "massError": -0.0437,
        "charge": 8,
        "mass": 4853.7224,
        "PPMerror": -9.0033,
        "thMass": 4853.7661,
        "peakId": "1",
        "matchedInd": "Y",
        "intensity": 443089.42
    },
    {
        "ion": "B43",
        "ionPos": "43",
        "position": 43,
        "massError": -0.042,
        "charge": 9,
        "mass": 4984.7646,
        "PPMerror": -8.4256,
        "thMass": 4984.8066,
        "peakId": "3",
        "matchedInd": "Y",
        "intensity": 135674.42
    },
    {
        "ion": "B42",
        "ionPos": "42",
        "position": 42,
        "massError": -0.0439,
        "charge": 7,
        "mass": 4853.7222,
        "PPMerror": -9.0445,
        "thMass": 4853.7661,
        "peakId": "6",
        "matchedInd": "Y",
        "intensity": 60366.41
    },
    {
        "ion": "B40",
        "ionPos": "40",
        "position": 40,
        "massError": -0.0397,
        "charge": 8,
        "mass": 4642.6307,
        "PPMerror": -8.5511,
        "thMass": 4642.6704,
        "peakId": "10",
        "matchedInd": "Y",
        "intensity": 63778.48
    },
    {
        "ion": "B11",
        "ionPos": "11",
        "position": 11,
        "massError": -0.0105,
        "charge": 2,
        "mass": 1238.703,
        "PPMerror": -8.4765,
        "thMass": 1238.7135,
        "peakId": "13",
        "matchedInd": "Y",
        "intensity": 40636.09
    },
    {
        "ion": "B43",
        "ionPos": "43",
        "position": 43,
        "massError": -0.0455,
        "charge": 7,
        "mass": 4984.7611,
        "PPMerror": -9.1277,
        "thMass": 4984.8066,
        "peakId": "15",
        "matchedInd": "Y",
        "intensity": 34729.58
    },
    {
        "ion": "B41",
        "ionPos": "41",
        "position": 41,
        "massError": -0.0359,
        "charge": 8,
        "mass": 4739.6873,
        "PPMerror": -7.5743,
        "thMass": 4739.7232,
        "peakId": "27",
        "matchedInd": "Y",
        "intensity": 14025.42
    },
    {
        "ion": "B23",
        "ionPos": "23",
        "position": 23,
        "massError": -0.0226,
        "charge": 4,
        "mass": 2626.5461,
        "PPMerror": -8.6044,
        "thMass": 2626.5687,
        "peakId": "29",
        "matchedInd": "Y",
        "intensity": 11396.92
    },
    {
        "ion": "B43",
        "ionPos": "43",
        "position": 43,
        "massError": -0.0422,
        "charge": 10,
        "mass": 4984.7644,
        "PPMerror": -8.4657,
        "thMass": 4984.8066,
        "peakId": "42",
        "matchedInd": "Y",
        "intensity": 10559.51
    },
    {
        "ion": "B26",
        "ionPos": "26",
        "position": 26,
        "massError": -0.0278,
        "charge": 4,
        "mass": 2885.6213,
        "PPMerror": -9.6339,
        "thMass": 2885.6491,
        "peakId": "44",
        "matchedInd": "Y",
        "intensity": 10903.77
    },
    {
        "ion": "B11",
        "ionPos": "11",
        "position": 11,
        "massError": -0.0103,
        "charge": 3,
        "mass": 1238.7032,
        "PPMerror": -8.3151,
        "thMass": 1238.7135,
        "peakId": "63",
        "matchedInd": "Y",
        "intensity": 11233.24
    },
    {
        "ion": "B56",
        "ionPos": "56",
        "position": 56,
        "massError": -0.072,
        "charge": 9,
        "mass": 6381.5033,
        "PPMerror": -11.2825,
        "thMass": 6381.5753,
        "peakId": "89",
        "matchedInd": "Y",
        "intensity": 4434.18
    },
    {
        "ion": "B13",
        "ionPos": "13",
        "position": 13,
        "massError": -0.0153,
        "charge": 3,
        "mass": 1488.8412,
        "PPMerror": -10.2763,
        "thMass": 1488.8565,
        "peakId": "93",
        "matchedInd": "Y",
        "intensity": 3088.8
    },
    {
        "ion": "Y47",
        "ionPos": "47",
        "position": 44,
        "massError": -0.0522,
        "charge": 9,
        "mass": 5180.7161,
        "PPMerror": -10.0757,
        "thMass": 5180.7683,
        "peakId": "2",
        "matchedInd": "Y",
        "intensity": 209861.94
    },
    {
        "ion": "Y30",
        "ionPos": "30",
        "position": 61,
        "massError": -0.0285,
        "charge": 6,
        "mass": 3355.7287,
        "PPMerror": -8.4929,
        "thMass": 3355.7572,
        "peakId": "4",
        "matchedInd": "Y",
        "intensity": 83840.49
    },
    {
        "ion": "Y34",
        "ionPos": "34",
        "position": 57,
        "massError": -0.0323,
        "charge": 6,
        "mass": 3797.9829,
        "PPMerror": -8.5044,
        "thMass": 3798.0152,
        "peakId": "9",
        "matchedInd": "Y",
        "intensity": 46902.57
    },
    {
        "ion": "Y45",
        "ionPos": "45",
        "position": 46,
        "massError": -0.0468,
        "charge": 8,
        "mass": 5010.616,
        "PPMerror": -9.3401,
        "thMass": 5010.6628,
        "peakId": "11",
        "matchedInd": "Y",
        "intensity": 52177.22
    },
    {
        "ion": "Y30",
        "ionPos": "30",
        "position": 61,
        "massError": -0.0299,
        "charge": 5,
        "mass": 3355.7273,
        "PPMerror": -8.9101,
        "thMass": 3355.7572,
        "peakId": "12",
        "matchedInd": "Y",
        "intensity": 36198.58
    },
    {
        "ion": "Y26",
        "ionPos": "26",
        "position": 65,
        "massError": -0.0236,
        "charge": 5,
        "mass": 2893.5915,
        "PPMerror": -8.1559,
        "thMass": 2893.6151,
        "peakId": "14",
        "matchedInd": "Y",
        "intensity": 37288.42
    },
    {
        "ion": "Y25",
        "ionPos": "25",
        "position": 66,
        "massError": -0.0238,
        "charge": 5,
        "mass": 2794.5229,
        "PPMerror": -8.5166,
        "thMass": 2794.5467,
        "peakId": "19",
        "matchedInd": "Y",
        "intensity": 35088.79
    },
    {
        "ion": "Y48",
        "ionPos": "48",
        "position": 43,
        "massError": -0.0518,
        "charge": 9,
        "mass": 5293.8006,
        "PPMerror": -9.7849,
        "thMass": 5293.8524,
        "peakId": "20",
        "matchedInd": "Y",
        "intensity": 22062.14
    },
    {
        "ion": "Y27",
        "ionPos": "27",
        "position": 64,
        "massError": -0.0275,
        "charge": 5,
        "mass": 3024.6281,
        "PPMerror": -9.0919,
        "thMass": 3024.6556,
        "peakId": "30",
        "matchedInd": "Y",
        "intensity": 18950.96
    },
    {
        "ion": "Y34",
        "ionPos": "34",
        "position": 57,
        "massError": -0.0304,
        "charge": 7,
        "mass": 3797.9848,
        "PPMerror": -8.0042,
        "thMass": 3798.0152,
        "peakId": "32",
        "matchedInd": "Y",
        "intensity": 10257.47
    },
    {
        "ion": "Y29",
        "ionPos": "29",
        "position": 62,
        "massError": -0.0294,
        "charge": 5,
        "mass": 3268.6958,
        "PPMerror": -8.9943,
        "thMass": 3268.7252,
        "peakId": "33",
        "matchedInd": "Y",
        "intensity": 13676.48
    },
    {
        "ion": "Y32",
        "ionPos": "32",
        "position": 59,
        "massError": -0.0258,
        "charge": 6,
        "mass": 3601.8682,
        "PPMerror": -7.1629,
        "thMass": 3601.894,
        "peakId": "35",
        "matchedInd": "Y",
        "intensity": 12944.01
    },
    {
        "ion": "Y26",
        "ionPos": "26",
        "position": 65,
        "massError": -0.0246,
        "charge": 6,
        "mass": 2893.5905,
        "PPMerror": -8.5015,
        "thMass": 2893.6151,
        "peakId": "37",
        "matchedInd": "Y",
        "intensity": 10203.04
    },
    {
        "ion": "Y6",
        "ionPos": "6",
        "position": 85,
        "massError": -0.0053,
        "charge": 2,
        "mass": 729.5172,
        "PPMerror": -7.265,
        "thMass": 729.5225,
        "peakId": "38",
        "matchedInd": "Y",
        "intensity": 23177.39
    },
    {
        "ion": "Y27",
        "ionPos": "27",
        "position": 64,
        "massError": -0.0258,
        "charge": 6,
        "mass": 3024.6298,
        "PPMerror": -8.5299,
        "thMass": 3024.6556,
        "peakId": "43",
        "matchedInd": "Y",
        "intensity": 12599.81
    },
    {
        "ion": "Y28",
        "ionPos": "28",
        "position": 63,
        "massError": -0.0289,
        "charge": 5,
        "mass": 3153.6693,
        "PPMerror": -9.1638,
        "thMass": 3153.6982,
        "peakId": "47",
        "matchedInd": "Y",
        "intensity": 12432.03
    },
    {
        "ion": "Y43",
        "ionPos": "43",
        "position": 48,
        "massError": -0.0404,
        "charge": 7,
        "mass": 4796.4906,
        "PPMerror": -8.4228,
        "thMass": 4796.531,
        "peakId": "49",
        "matchedInd": "Y",
        "intensity": 14460.74
    },
    {
        "ion": "Y43",
        "ionPos": "43",
        "position": 48,
        "massError": -0.0351,
        "charge": 8,
        "mass": 4796.4959,
        "PPMerror": -7.3178,
        "thMass": 4796.531,
        "peakId": "51",
        "matchedInd": "Y",
        "intensity": 7223.97
    },
    {
        "ion": "Y25",
        "ionPos": "25",
        "position": 66,
        "massError": -0.0239,
        "charge": 6,
        "mass": 2794.5228,
        "PPMerror": -8.5524,
        "thMass": 2794.5467,
        "peakId": "77",
        "matchedInd": "Y",
        "intensity": 6131.37
    },
    {
        "ion": "Y44",
        "ionPos": "44",
        "position": 47,
        "massError": -0.0426,
        "charge": 8,
        "mass": 4909.5725,
        "PPMerror": -8.6769,
        "thMass": 4909.6151,
        "peakId": "79",
        "matchedInd": "Y",
        "intensity": 5312.34
    }
];
