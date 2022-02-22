"use strict";
let sampleMonoPeaks = [
    {
        "displayLevel_": 7,
        "peakId_": "0",
        "pos_": 4853.7224,
        "monoMass_": 4853.7224,
        "monoMz_": 607.7226,
        "charge_": 8,
        "intensity_": 443089.42
    },
    {
        "displayLevel_": 7,
        "peakId_": "1",
        "pos_": 5180.7161,
        "monoMass_": 5180.7161,
        "monoMz_": 576.6424,
        "charge_": 9,
        "intensity_": 209861.94
    },
    {
        "displayLevel_": 6,
        "peakId_": "2",
        "pos_": 4984.7646,
        "monoMass_": 4984.7646,
        "monoMz_": 554.87,
        "charge_": 9,
        "intensity_": 135674.42
    },
    {
        "displayLevel_": 5,
        "peakId_": "4",
        "pos_": 5162.7112,
        "monoMass_": 5162.7112,
        "monoMz_": 574.6419,
        "charge_": 9,
        "intensity_": 87181.44
    },
    {
        "displayLevel_": 7,
        "peakId_": "3",
        "pos_": 3355.7287,
        "monoMass_": 3355.7287,
        "monoMz_": 560.2954,
        "charge_": 6,
        "intensity_": 83840.49
    },
    {
        "displayLevel_": 7,
        "peakId_": "24",
        "pos_": 7284.0832,
        "monoMass_": 7284.0832,
        "monoMz_": 608.0142,
        "charge_": 12,
        "intensity_": 73984.46
    },
    {
        "displayLevel_": 6,
        "peakId_": "9",
        "pos_": 4642.6307,
        "monoMass_": 4642.6307,
        "monoMz_": 581.3361,
        "charge_": 8,
        "intensity_": 63778.48
    },
    {
        "displayLevel_": 7,
        "peakId_": "6",
        "pos_": 3615.0214,
        "monoMass_": 3615.0214,
        "monoMz_": 603.5108,
        "charge_": 6,
        "intensity_": 61517.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "5",
        "pos_": 4853.7222,
        "monoMass_": 4853.7222,
        "monoMz_": 694.3962,
        "charge_": 7,
        "intensity_": 60366.41
    },
    {
        "displayLevel_": 5,
        "peakId_": "10",
        "pos_": 5010.616,
        "monoMass_": 5010.616,
        "monoMz_": 627.3343,
        "charge_": 8,
        "intensity_": 52177.22
    },
    {
        "displayLevel_": 6,
        "peakId_": "8",
        "pos_": 3797.9829,
        "monoMass_": 3797.9829,
        "monoMz_": 634.0044,
        "charge_": 6,
        "intensity_": 46902.57
    },
    {
        "displayLevel_": 0,
        "peakId_": "7",
        "pos_": 4983.7611,
        "monoMass_": 4983.7611,
        "monoMz_": 623.9774,
        "charge_": 8,
        "intensity_": 44365.45
    },
    {
        "displayLevel_": 7,
        "peakId_": "12",
        "pos_": 1238.703,
        "monoMass_": 1238.703,
        "monoMz_": 620.3588,
        "charge_": 2,
        "intensity_": 40636.09
    },
    {
        "displayLevel_": 7,
        "peakId_": "13",
        "pos_": 2893.5915,
        "monoMass_": 2893.5915,
        "monoMz_": 579.7256,
        "charge_": 5,
        "intensity_": 37288.42
    },
    {
        "displayLevel_": 0,
        "peakId_": "11",
        "pos_": 3355.7273,
        "monoMass_": 3355.7273,
        "monoMz_": 672.1527,
        "charge_": 5,
        "intensity_": 36198.58
    },
    {
        "displayLevel_": 6,
        "peakId_": "18",
        "pos_": 2794.5229,
        "monoMass_": 2794.5229,
        "monoMz_": 559.9118,
        "charge_": 5,
        "intensity_": 35088.79
    },
    {
        "displayLevel_": 0,
        "peakId_": "14",
        "pos_": 4984.7611,
        "monoMass_": 4984.7611,
        "monoMz_": 713.116,
        "charge_": 7,
        "intensity_": 34729.58
    },
    {
        "displayLevel_": 0,
        "peakId_": "15",
        "pos_": 5179.7205,
        "monoMass_": 5179.7205,
        "monoMz_": 648.4723,
        "charge_": 8,
        "intensity_": 34087.68
    },
    {
        "displayLevel_": 6,
        "peakId_": "20",
        "pos_": 1382.7429,
        "monoMass_": 1382.7429,
        "monoMz_": 692.3787,
        "charge_": 2,
        "intensity_": 27332
    },
    {
        "displayLevel_": 6,
        "peakId_": "72",
        "pos_": 1151.2709,
        "monoMass_": 1151.2709,
        "monoMz_": 576.6427,
        "charge_": 2,
        "intensity_": 23960.12
    },
    {
        "displayLevel_": 7,
        "peakId_": "37",
        "pos_": 729.5172,
        "monoMass_": 729.5172,
        "monoMz_": 365.7659,
        "charge_": 2,
        "intensity_": 23177.39
    },
    {
        "displayLevel_": 6,
        "peakId_": "19",
        "pos_": 5293.8006,
        "monoMass_": 5293.8006,
        "monoMz_": 589.2073,
        "charge_": 9,
        "intensity_": 22062.14
    },
    {
        "displayLevel_": 4,
        "peakId_": "21",
        "pos_": 3337.7214,
        "monoMass_": 3337.7214,
        "monoMz_": 557.2942,
        "charge_": 6,
        "intensity_": 20303.89
    },
    {
        "displayLevel_": 5,
        "peakId_": "17",
        "pos_": 4879.7352,
        "monoMass_": 4879.7352,
        "monoMz_": 610.9742,
        "charge_": 8,
        "intensity_": 19255.85
    },
    {
        "displayLevel_": 6,
        "peakId_": "29",
        "pos_": 3024.6281,
        "monoMass_": 3024.6281,
        "monoMz_": 605.9329,
        "charge_": 5,
        "intensity_": 18950.96
    },
    {
        "displayLevel_": 7,
        "peakId_": "33",
        "pos_": 2155.0903,
        "monoMass_": 2155.0903,
        "monoMz_": 719.3707,
        "charge_": 3,
        "intensity_": 18055.96
    },
    {
        "displayLevel_": 6,
        "peakId_": "22",
        "pos_": 3403.9316,
        "monoMass_": 3403.9316,
        "monoMz_": 568.3292,
        "charge_": 6,
        "intensity_": 17087.81
    },
    {
        "displayLevel_": 4,
        "peakId_": "23",
        "pos_": 5144.6953,
        "monoMass_": 5144.6953,
        "monoMz_": 572.6401,
        "charge_": 9,
        "intensity_": 16816.8
    },
    {
        "displayLevel_": 4,
        "peakId_": "16",
        "pos_": 3779.9683,
        "monoMass_": 3779.9683,
        "monoMz_": 631.002,
        "charge_": 6,
        "intensity_": 16647.54
    },
    {
        "displayLevel_": 7,
        "peakId_": "27",
        "pos_": 4450.2055,
        "monoMass_": 4450.2055,
        "monoMz_": 636.7509,
        "charge_": 7,
        "intensity_": 15880.08
    },
    {
        "displayLevel_": 7,
        "peakId_": "35",
        "pos_": 6658.6889,
        "monoMass_": 6658.6889,
        "monoMz_": 666.8762,
        "charge_": 10,
        "intensity_": 15735.78
    },
    {
        "displayLevel_": 0,
        "peakId_": "38",
        "pos_": 6659.7111,
        "monoMass_": 6659.7111,
        "monoMz_": 606.4356,
        "charge_": 11,
        "intensity_": 15474.4
    },
    {
        "displayLevel_": 5,
        "peakId_": "40",
        "pos_": 1283.6727,
        "monoMass_": 1283.6727,
        "monoMz_": 642.8436,
        "charge_": 2,
        "intensity_": 15058.63
    },
    {
        "displayLevel_": 2,
        "peakId_": "49",
        "pos_": 4880.7448,
        "monoMass_": 4880.7448,
        "monoMz_": 543.3123,
        "charge_": 9,
        "intensity_": 14599.78
    },
    {
        "displayLevel_": 6,
        "peakId_": "48",
        "pos_": 4796.4906,
        "monoMass_": 4796.4906,
        "monoMz_": 686.2202,
        "charge_": 7,
        "intensity_": 14460.74
    },
    {
        "displayLevel_": 5,
        "peakId_": "26",
        "pos_": 4739.6873,
        "monoMass_": 4739.6873,
        "monoMz_": 593.4682,
        "charge_": 8,
        "intensity_": 14025.42
    },
    {
        "displayLevel_": 2,
        "peakId_": "57",
        "pos_": 4455.3645,
        "monoMass_": 4455.3645,
        "monoMz_": 557.9278,
        "charge_": 8,
        "intensity_": 13759.93
    },
    {
        "displayLevel_": 5,
        "peakId_": "32",
        "pos_": 3268.6958,
        "monoMass_": 3268.6958,
        "monoMz_": 654.7464,
        "charge_": 5,
        "intensity_": 13676.48
    },
    {
        "displayLevel_": 5,
        "peakId_": "34",
        "pos_": 3601.8682,
        "monoMass_": 3601.8682,
        "monoMz_": 601.3186,
        "charge_": 6,
        "intensity_": 12944.01
    },
    {
        "displayLevel_": 0,
        "peakId_": "42",
        "pos_": 3024.6298,
        "monoMass_": 3024.6298,
        "monoMz_": 505.1122,
        "charge_": 6,
        "intensity_": 12599.81
    },
    {
        "displayLevel_": 6,
        "peakId_": "46",
        "pos_": 3153.6693,
        "monoMass_": 3153.6693,
        "monoMz_": 631.7411,
        "charge_": 5,
        "intensity_": 12432.03
    },
    {
        "displayLevel_": 0,
        "peakId_": "25",
        "pos_": 3337.7235,
        "monoMass_": 3337.7235,
        "monoMz_": 668.552,
        "charge_": 5,
        "intensity_": 12223.63
    },
    {
        "displayLevel_": 7,
        "peakId_": "45",
        "pos_": 1967.0958,
        "monoMass_": 1967.0958,
        "monoMz_": 492.7812,
        "charge_": 4,
        "intensity_": 12218.46
    },
    {
        "displayLevel_": 5,
        "peakId_": "30",
        "pos_": 2027.0515,
        "monoMass_": 2027.0515,
        "monoMz_": 676.6911,
        "charge_": 3,
        "intensity_": 11757.51
    },
    {
        "displayLevel_": 6,
        "peakId_": "28",
        "pos_": 2626.5461,
        "monoMass_": 2626.5461,
        "monoMz_": 657.6438,
        "charge_": 4,
        "intensity_": 11396.92
    },
    {
        "displayLevel_": 6,
        "peakId_": "125",
        "pos_": 2308.318,
        "monoMass_": 2308.318,
        "monoMz_": 578.0868,
        "charge_": 4,
        "intensity_": 11376.32
    },
    {
        "displayLevel_": 0,
        "peakId_": "62",
        "pos_": 1238.7032,
        "monoMass_": 1238.7032,
        "monoMz_": 413.9083,
        "charge_": 3,
        "intensity_": 11233.24
    },
    {
        "displayLevel_": 4,
        "peakId_": "43",
        "pos_": 2885.6213,
        "monoMass_": 2885.6213,
        "monoMz_": 722.4126,
        "charge_": 4,
        "intensity_": 10903.77
    },
    {
        "displayLevel_": 0,
        "peakId_": "41",
        "pos_": 4984.7644,
        "monoMass_": 4984.7644,
        "monoMz_": 499.4837,
        "charge_": 10,
        "intensity_": 10559.51
    },
    {
        "displayLevel_": 6,
        "peakId_": "51",
        "pos_": 3937.1958,
        "monoMass_": 3937.1958,
        "monoMz_": 657.2066,
        "charge_": 6,
        "intensity_": 10421.07
    },
    {
        "displayLevel_": 5,
        "peakId_": "71",
        "pos_": 7236.0671,
        "monoMass_": 7236.0671,
        "monoMz_": 604.0129,
        "charge_": 12,
        "intensity_": 10415.49
    },
    {
        "displayLevel_": 0,
        "peakId_": "31",
        "pos_": 3797.9848,
        "monoMass_": 3797.9848,
        "monoMz_": 543.5765,
        "charge_": 7,
        "intensity_": 10257.47
    },
    {
        "displayLevel_": 0,
        "peakId_": "36",
        "pos_": 2893.5905,
        "monoMass_": 2893.5905,
        "monoMz_": 483.2724,
        "charge_": 6,
        "intensity_": 10203.04
    },
    {
        "displayLevel_": 4,
        "peakId_": "52",
        "pos_": 2776.5077,
        "monoMass_": 2776.5077,
        "monoMz_": 556.3088,
        "charge_": 5,
        "intensity_": 9759.65
    },
    {
        "displayLevel_": 5,
        "peakId_": "39",
        "pos_": 1313.1062,
        "monoMass_": 1313.1062,
        "monoMz_": 657.5604,
        "charge_": 2,
        "intensity_": 8990.46
    },
    {
        "displayLevel_": 4,
        "peakId_": "58",
        "pos_": 3135.6496,
        "monoMass_": 3135.6496,
        "monoMz_": 628.1372,
        "charge_": 5,
        "intensity_": 8882.84
    },
    {
        "displayLevel_": 6,
        "peakId_": "56",
        "pos_": 4433.19,
        "monoMass_": 4433.19,
        "monoMz_": 634.3201,
        "charge_": 7,
        "intensity_": 8841.98
    },
    {
        "displayLevel_": 0,
        "peakId_": "54",
        "pos_": 1382.7414,
        "monoMass_": 1382.7414,
        "monoMz_": 461.9211,
        "charge_": 3,
        "intensity_": 8455.83
    },
    {
        "displayLevel_": 5,
        "peakId_": "47",
        "pos_": 3230.7874,
        "monoMass_": 3230.7874,
        "monoMz_": 647.1647,
        "charge_": 5,
        "intensity_": 8122.1
    },
    {
        "displayLevel_": 4,
        "peakId_": "107",
        "pos_": 7213.1259,
        "monoMass_": 7213.1259,
        "monoMz_": 602.1011,
        "charge_": 12,
        "intensity_": 8102.25
    },
    {
        "displayLevel_": 5,
        "peakId_": "55",
        "pos_": 2065.007,
        "monoMass_": 2065.007,
        "monoMz_": 689.3429,
        "charge_": 3,
        "intensity_": 7993.43
    },
    {
        "displayLevel_": 6,
        "peakId_": "64",
        "pos_": 2406.1443,
        "monoMass_": 2406.1443,
        "monoMz_": 602.5433,
        "charge_": 4,
        "intensity_": 7942.47
    },
    {
        "displayLevel_": 6,
        "peakId_": "102",
        "pos_": 4220.3709,
        "monoMass_": 4220.3709,
        "monoMz_": 603.9174,
        "charge_": 7,
        "intensity_": 7823.08
    },
    {
        "displayLevel_": 5,
        "peakId_": "44",
        "pos_": 3068.4682,
        "monoMass_": 3068.4682,
        "monoMz_": 768.1243,
        "charge_": 4,
        "intensity_": 7721.96
    },
    {
        "displayLevel_": 4,
        "peakId_": "61",
        "pos_": 1364.7322,
        "monoMass_": 1364.7322,
        "monoMz_": 683.3734,
        "charge_": 2,
        "intensity_": 7692.52
    },
    {
        "displayLevel_": 4,
        "peakId_": "79",
        "pos_": 711.5067,
        "monoMass_": 711.5067,
        "monoMz_": 356.7606,
        "charge_": 2,
        "intensity_": 7359.43
    },
    {
        "displayLevel_": 0,
        "peakId_": "50",
        "pos_": 4796.4959,
        "monoMass_": 4796.4959,
        "monoMz_": 600.5693,
        "charge_": 8,
        "intensity_": 7223.97
    },
    {
        "displayLevel_": 4,
        "peakId_": "59",
        "pos_": 2009.0447,
        "monoMass_": 2009.0447,
        "monoMz_": 670.6888,
        "charge_": 3,
        "intensity_": 7131.53
    },
    {
        "displayLevel_": 0,
        "peakId_": "73",
        "pos_": 4878.7366,
        "monoMass_": 4878.7366,
        "monoMz_": 697.9697,
        "charge_": 7,
        "intensity_": 7051.96
    },
    {
        "displayLevel_": 7,
        "peakId_": "80",
        "pos_": 7854.9237,
        "monoMass_": 7854.9237,
        "monoMz_": 605.2322,
        "charge_": 13,
        "intensity_": 7017.25
    },
    {
        "displayLevel_": 4,
        "peakId_": "67",
        "pos_": 4836.6702,
        "monoMass_": 4836.6702,
        "monoMz_": 605.5911,
        "charge_": 8,
        "intensity_": 6999.3
    },
    {
        "displayLevel_": 5,
        "peakId_": "68",
        "pos_": 3747.0503,
        "monoMass_": 3747.0503,
        "monoMz_": 625.5157,
        "charge_": 6,
        "intensity_": 6327.23
    },
    {
        "displayLevel_": 3,
        "peakId_": "94",
        "pos_": 2017.0879,
        "monoMass_": 2017.0879,
        "monoMz_": 505.2793,
        "charge_": 4,
        "intensity_": 6256.65
    },
    {
        "displayLevel_": 0,
        "peakId_": "76",
        "pos_": 2794.5228,
        "monoMass_": 2794.5228,
        "monoMz_": 466.7611,
        "charge_": 6,
        "intensity_": 6131.37
    },
    {
        "displayLevel_": 0,
        "peakId_": "53",
        "pos_": 3354.7251,
        "monoMass_": 3354.7251,
        "monoMz_": 480.2537,
        "charge_": 7,
        "intensity_": 6025.54
    },
    {
        "displayLevel_": 5,
        "peakId_": "85",
        "pos_": 2496.8016,
        "monoMass_": 2496.8016,
        "monoMz_": 625.2077,
        "charge_": 4,
        "intensity_": 5750.4
    },
    {
        "displayLevel_": 4,
        "peakId_": "69",
        "pos_": 4896.7571,
        "monoMass_": 4896.7571,
        "monoMz_": 613.1019,
        "charge_": 8,
        "intensity_": 5734.85
    },
    {
        "displayLevel_": 6,
        "peakId_": "96",
        "pos_": 6906.5511,
        "monoMass_": 6906.5511,
        "monoMz_": 691.6624,
        "charge_": 10,
        "intensity_": 5699.91
    },
    {
        "displayLevel_": 5,
        "peakId_": "74",
        "pos_": 5052.836,
        "monoMass_": 5052.836,
        "monoMz_": 506.2909,
        "charge_": 10,
        "intensity_": 5470.44
    },
    {
        "displayLevel_": 7,
        "peakId_": "91",
        "pos_": 9362.6403,
        "monoMass_": 9362.6403,
        "monoMz_": 625.1833,
        "charge_": 15,
        "intensity_": 5404.92
    },
    {
        "displayLevel_": 1,
        "peakId_": "60",
        "pos_": 4854.7301,
        "monoMass_": 4854.7301,
        "monoMz_": 486.4803,
        "charge_": 10,
        "intensity_": 5359.37
    },
    {
        "displayLevel_": 5,
        "peakId_": "78",
        "pos_": 4909.5725,
        "monoMass_": 4909.5725,
        "monoMz_": 614.7038,
        "charge_": 8,
        "intensity_": 5312.34
    },
    {
        "displayLevel_": 0,
        "peakId_": "82",
        "pos_": 4852.7073,
        "monoMass_": 4852.7073,
        "monoMz_": 540.197,
        "charge_": 9,
        "intensity_": 5142.61
    },
    {
        "displayLevel_": 5,
        "peakId_": "77",
        "pos_": 1408.7588,
        "monoMass_": 1408.7588,
        "monoMz_": 705.3867,
        "charge_": 2,
        "intensity_": 5002.41
    },
    {
        "displayLevel_": 5,
        "peakId_": "83",
        "pos_": 4279.0763,
        "monoMass_": 4279.0763,
        "monoMz_": 714.1867,
        "charge_": 6,
        "intensity_": 5001.97
    },
    {
        "displayLevel_": 4,
        "peakId_": "63",
        "pos_": 3006.6081,
        "monoMass_": 3006.6081,
        "monoMz_": 602.3289,
        "charge_": 5,
        "intensity_": 5000.77
    },
    {
        "displayLevel_": 4,
        "peakId_": "90",
        "pos_": 4778.4493,
        "monoMass_": 4778.4493,
        "monoMz_": 683.6429,
        "charge_": 7,
        "intensity_": 4984.11
    },
    {
        "displayLevel_": 4,
        "peakId_": "84",
        "pos_": 3582.856,
        "monoMass_": 3582.856,
        "monoMz_": 598.1499,
        "charge_": 6,
        "intensity_": 4875.34
    },
    {
        "displayLevel_": 4,
        "peakId_": "65",
        "pos_": 4967.7346,
        "monoMass_": 4967.7346,
        "monoMz_": 621.9741,
        "charge_": 8,
        "intensity_": 4820.27
    },
    {
        "displayLevel_": 0,
        "peakId_": "114",
        "pos_": 1364.7334,
        "monoMass_": 1364.7334,
        "monoMz_": 455.9184,
        "charge_": 3,
        "intensity_": 4776.39
    },
    {
        "displayLevel_": 0,
        "peakId_": "98",
        "pos_": 10022.255,
        "monoMass_": 10022.255,
        "monoMz_": 557.7992,
        "charge_": 18,
        "intensity_": 4746.9
    },
    {
        "displayLevel_": 3,
        "peakId_": "75",
        "pos_": 4957.8003,
        "monoMass_": 4957.8003,
        "monoMz_": 551.874,
        "charge_": 9,
        "intensity_": 4615.72
    },
    {
        "displayLevel_": 6,
        "peakId_": "70",
        "pos_": 1725.9223,
        "monoMass_": 1725.9223,
        "monoMz_": 863.9684,
        "charge_": 2,
        "intensity_": 4603.61
    },
    {
        "displayLevel_": 0,
        "peakId_": "112",
        "pos_": 2626.2167,
        "monoMass_": 2626.2167,
        "monoMz_": 876.4128,
        "charge_": 3,
        "intensity_": 4490.93
    },
    {
        "displayLevel_": 6,
        "peakId_": "88",
        "pos_": 6381.5033,
        "monoMass_": 6381.5033,
        "monoMz_": 710.0632,
        "charge_": 9,
        "intensity_": 4434.18
    },
    {
        "displayLevel_": 3,
        "peakId_": "81",
        "pos_": 5155.8419,
        "monoMass_": 5155.8419,
        "monoMz_": 645.4875,
        "charge_": 8,
        "intensity_": 4413.03
    },
    {
        "displayLevel_": 4,
        "peakId_": "89",
        "pos_": 5275.7863,
        "monoMass_": 5275.7863,
        "monoMz_": 587.2058,
        "charge_": 9,
        "intensity_": 4367.13
    },
    {
        "displayLevel_": 4,
        "peakId_": "87",
        "pos_": 1707.9215,
        "monoMass_": 1707.9215,
        "monoMz_": 570.3145,
        "charge_": 3,
        "intensity_": 4264.48
    },
    {
        "displayLevel_": 6,
        "peakId_": "95",
        "pos_": 651.3221,
        "monoMass_": 651.3221,
        "monoMz_": 652.3293,
        "charge_": 1,
        "intensity_": 4147.42
    },
    {
        "displayLevel_": 0,
        "peakId_": "120",
        "pos_": 3067.3561,
        "monoMass_": 3067.3561,
        "monoMz_": 614.4785,
        "charge_": 5,
        "intensity_": 4109.29
    },
    {
        "displayLevel_": 5,
        "peakId_": "110",
        "pos_": 6298.7657,
        "monoMass_": 6298.7657,
        "monoMz_": 573.6223,
        "charge_": 11,
        "intensity_": 3943.89
    },
    {
        "displayLevel_": 5,
        "peakId_": "129",
        "pos_": 4337.676,
        "monoMass_": 4337.676,
        "monoMz_": 620.6753,
        "charge_": 7,
        "intensity_": 3833.76
    },
    {
        "displayLevel_": 7,
        "peakId_": "118",
        "pos_": 5705.7142,
        "monoMass_": 5705.7142,
        "monoMz_": 571.5787,
        "charge_": 10,
        "intensity_": 3830
    },
    {
        "displayLevel_": 1,
        "peakId_": "100",
        "pos_": 713.0284,
        "monoMass_": 713.0284,
        "monoMz_": 714.0356,
        "charge_": 1,
        "intensity_": 3823.42
    },
    {
        "displayLevel_": 5,
        "peakId_": "104",
        "pos_": 1784.0734,
        "monoMass_": 1784.0734,
        "monoMz_": 199.2377,
        "charge_": 9,
        "intensity_": 3737.81
    },
    {
        "displayLevel_": 5,
        "peakId_": "99",
        "pos_": 579.5211,
        "monoMass_": 579.5211,
        "monoMz_": 580.5284,
        "charge_": 1,
        "intensity_": 3716.05
    },
    {
        "displayLevel_": 5,
        "peakId_": "117",
        "pos_": 4563.6582,
        "monoMass_": 4563.6582,
        "monoMz_": 508.0804,
        "charge_": 9,
        "intensity_": 3702.48
    },
    {
        "displayLevel_": 6,
        "peakId_": "124",
        "pos_": 6014.9039,
        "monoMass_": 6014.9039,
        "monoMz_": 669.3299,
        "charge_": 9,
        "intensity_": 3580.85
    },
    {
        "displayLevel_": 5,
        "peakId_": "113",
        "pos_": 3186.3629,
        "monoMass_": 3186.3629,
        "monoMz_": 638.2799,
        "charge_": 5,
        "intensity_": 3439.79
    },
    {
        "displayLevel_": 3,
        "peakId_": "103",
        "pos_": 2875.5882,
        "monoMass_": 2875.5882,
        "monoMz_": 480.272,
        "charge_": 6,
        "intensity_": 3408.01
    },
    {
        "displayLevel_": 4,
        "peakId_": "86",
        "pos_": 1951.0758,
        "monoMass_": 1951.0758,
        "monoMz_": 488.7762,
        "charge_": 4,
        "intensity_": 3373.24
    },
    {
        "displayLevel_": 7,
        "peakId_": "126",
        "pos_": 8430.1097,
        "monoMass_": 8430.1097,
        "monoMz_": 563.0146,
        "charge_": 15,
        "intensity_": 3195.7
    },
    {
        "displayLevel_": 6,
        "peakId_": "92",
        "pos_": 1488.8412,
        "monoMass_": 1488.8412,
        "monoMz_": 497.2877,
        "charge_": 3,
        "intensity_": 3088.8
    },
    {
        "displayLevel_": 4,
        "peakId_": "66",
        "pos_": 632.9969,
        "monoMass_": 632.9969,
        "monoMz_": 634.0041,
        "charge_": 1,
        "intensity_": 2845.67
    },
    {
        "displayLevel_": 6,
        "peakId_": "109",
        "pos_": 5392.0844,
        "monoMass_": 5392.0844,
        "monoMz_": 600.1278,
        "charge_": 9,
        "intensity_": 2825.44
    },
    {
        "displayLevel_": 4,
        "peakId_": "105",
        "pos_": 2290.4883,
        "monoMass_": 2290.4883,
        "monoMz_": 573.6293,
        "charge_": 4,
        "intensity_": 2809.82
    },
    {
        "displayLevel_": 4,
        "peakId_": "123",
        "pos_": 593.4544,
        "monoMass_": 593.4544,
        "monoMz_": 198.8254,
        "charge_": 3,
        "intensity_": 2793.16
    },
    {
        "displayLevel_": 0,
        "peakId_": "116",
        "pos_": 3799.0028,
        "monoMass_": 3799.0028,
        "monoMz_": 760.8078,
        "charge_": 5,
        "intensity_": 2719.53
    },
    {
        "displayLevel_": 5,
        "peakId_": "111",
        "pos_": 2379.3623,
        "monoMass_": 2379.3623,
        "monoMz_": 595.8478,
        "charge_": 4,
        "intensity_": 2668.52
    },
    {
        "displayLevel_": 4,
        "peakId_": "115",
        "pos_": 3730.0347,
        "monoMass_": 3730.0347,
        "monoMz_": 622.6797,
        "charge_": 6,
        "intensity_": 2656.4
    },
    {
        "displayLevel_": 5,
        "peakId_": "97",
        "pos_": 5099.8561,
        "monoMass_": 5099.8561,
        "monoMz_": 567.658,
        "charge_": 9,
        "intensity_": 2617.21
    },
    {
        "displayLevel_": 0,
        "peakId_": "108",
        "pos_": 4779.5209,
        "monoMass_": 4779.5209,
        "monoMz_": 598.4474,
        "charge_": 8,
        "intensity_": 2597.79
    },
    {
        "displayLevel_": 3,
        "peakId_": "122",
        "pos_": 4555.8231,
        "monoMass_": 4555.8231,
        "monoMz_": 570.4852,
        "charge_": 8,
        "intensity_": 2516.21
    },
    {
        "displayLevel_": 6,
        "peakId_": "93",
        "pos_": 970.9439,
        "monoMass_": 970.9439,
        "monoMz_": 486.4792,
        "charge_": 2,
        "intensity_": 2411.68
    },
    {
        "displayLevel_": 4,
        "peakId_": "101",
        "pos_": 1296.5168,
        "monoMass_": 1296.5168,
        "monoMz_": 649.2657,
        "charge_": 2,
        "intensity_": 2395.92
    },
    {
        "displayLevel_": 5,
        "peakId_": "121",
        "pos_": 1095.6816,
        "monoMass_": 1095.6816,
        "monoMz_": 366.2345,
        "charge_": 3,
        "intensity_": 2302.72
    },
    {
        "displayLevel_": 0,
        "peakId_": "119",
        "pos_": 1408.752,
        "monoMass_": 1408.752,
        "monoMz_": 470.5913,
        "charge_": 3,
        "intensity_": 2186.59
    },
    {
        "displayLevel_": 5,
        "peakId_": "106",
        "pos_": 1210.7113,
        "monoMass_": 1210.7113,
        "monoMz_": 404.5777,
        "charge_": 3,
        "intensity_": 1786.02
    },
    {
        "displayLevel_": 4,
        "peakId_": "127",
        "pos_": 1175.374,
        "monoMass_": 1175.374,
        "monoMz_": 588.6943,
        "charge_": 2,
        "intensity_": 1777.5
    },
    {
        "displayLevel_": 7,
        "peakId_": "128",
        "pos_": 555.7954,
        "monoMass_": 555.7954,
        "monoMz_": 556.8027,
        "charge_": 1,
        "intensity_": 1581.36
    }
];
