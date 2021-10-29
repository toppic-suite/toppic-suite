"use strict";
let sampleMatchedIons = [
    {
        "mz": 348.20245361328125,
        "intensity": 178241.99734228774,
        "text": "B6",
        "error": -0.0045,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 7,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 348.20245361328125,
                    "monoMz_": 348.20245361328125,
                    "intensity_": 178241.99734228774
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 348.70389461328125,
                    "monoMz_": 348.70389461328125,
                    "intensity_": 69661.40187500001
                }
            ],
            "monoMass_": 694.3903542928045,
            "charge_": 2,
            "intensity_": -1
        }
    },
    {
        "mz": 390.21282958984375,
        "intensity": 33501.95410156249,
        "text": "B10",
        "error": -0.0056,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 4,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 390.21282958984375,
                    "monoMz_": 390.21282958984375,
                    "intensity_": 33501.95410156249
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 390.5471359231771,
                    "monoMz_": 390.5471359231771,
                    "intensity_": 21830.865620458604
                }
            ],
            "monoMass_": 1167.6166593688943,
            "charge_": 3,
            "intensity_": -1
        }
    },
    {
        "mz": 453.74981689453125,
        "intensity": 116528.892578125,
        "text": "B8",
        "error": -0.0055,
        "env": {
            "displayColor_": "darkorange",
            "displayLevel_": 5,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 453.74981689453125,
                    "monoMz_": 453.74981689453125,
                    "intensity_": 116528.892578125
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 454.25127789453126,
                    "monoMz_": 454.25127789453126,
                    "intensity_": 59871.618601944625
                }
            ],
            "monoMass_": 905.4850808553045,
            "charge_": 2,
            "intensity_": -1
        }
    },
    {
        "mz": 510.0311279296875,
        "intensity": 1142514.9105677975,
        "text": "Y18",
        "error": -0.015,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 7,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 510.0311279296875,
                    "monoMz_": 510.0311279296875,
                    "intensity_": 1142514.9105677975
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 509.7804079296875,
                    "monoMz_": 509.7804079296875,
                    "intensity_": 1003233.1800000002
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 510.28167017968747,
                    "monoMz_": 510.28167017968747,
                    "intensity_": 744438.2659068907
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 510.5321449296875,
                    "monoMz_": 510.5321449296875,
                    "intensity_": 353537.44892023434
                }
            ],
            "monoMass_": 2035.092525851234,
            "charge_": 4,
            "intensity_": -1
        }
    },
    {
        "mz": 535.2796020507812,
        "intensity": 216742.64105173617,
        "text": "B9",
        "error": -0.0092,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 5,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 535.2796020507812,
                    "monoMz_": 535.2796020507812,
                    "intensity_": 216742.64105173617
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 535.7810600507812,
                    "monoMz_": 535.7810600507812,
                    "intensity_": 130425.4109375
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 536.2824050507812,
                    "monoMz_": 536.2824050507812,
                    "intensity_": 44763.61273723068
                }
            ],
            "monoMass_": 1068.5446511678044,
            "charge_": 2,
            "intensity_": -1
        }
    },
    {
        "mz": 579.4706420898438,
        "intensity": 90980.54405307952,
        "text": "Y56",
        "error": -0.0007,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 0,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 579.4706420898438,
                    "monoMz_": 579.4706420898438,
                    "intensity_": 90980.54405307952
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 579.5617759989346,
                    "monoMz_": 579.5617759989346,
                    "intensity_": 89989.5393867868
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 579.6529020898438,
                    "monoMz_": 579.6529020898438,
                    "intensity_": 73156.52019308838
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 579.3794982716619,
                    "monoMz_": 579.3794982716619,
                    "intensity_": 71071.99216385302
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "5",
                    "pos_": 579.7440219080255,
                    "monoMz_": 579.7440219080255,
                    "intensity_": 50796.23279096846
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 579.2883419080255,
                    "monoMz_": 579.2883419080255,
                    "intensity_": 38257.4325
                }
            ],
            "monoMass_": 6360.088825852614,
            "charge_": 11,
            "intensity_": -1
        }
    },
    {
        "mz": 584.8133544921875,
        "intensity": 51735.1578125,
        "text": "B10",
        "error": -0.0101,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 0,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 584.8133544921875,
                    "monoMz_": 584.8133544921875,
                    "intensity_": 51735.1578125
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 585.3148139921875,
                    "monoMz_": 585.3148139921875,
                    "intensity_": 33712.1612259994
                }
            ],
            "monoMass_": 1167.612156050617,
            "charge_": 2,
            "intensity_": -1
        }
    },
    {
        "mz": 584.9208374023438,
        "intensity": 43490.39234375,
        "text": "Y26",
        "error": -0.0302,
        "env": {
            "displayColor_": "darkorange",
            "displayLevel_": 0,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 584.9208374023438,
                    "monoMz_": 584.9208374023438,
                    "intensity_": 43490.39234375
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 585.1213360023437,
                    "monoMz_": 585.1213360023437,
                    "intensity_": 38458.30996847707
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 584.7202602023438,
                    "monoMz_": 584.7202602023438,
                    "intensity_": 26693.306862706686
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 585.3217894023437,
                    "monoMz_": 585.3217894023437,
                    "intensity_": 24168.509971831616
                }
            ],
            "monoMass_": 2918.5649186773235,
            "charge_": 5,
            "intensity_": -1
        }
    },
    {
        "mz": 585.9281616210938,
        "intensity": 79228.11902553569,
        "text": "Y57",
        "error": -0.0051,
        "env": {
            "displayColor_": "blue",
            "displayLevel_": 3,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 585.9281616210938,
                    "monoMz_": 585.9281616210938,
                    "intensity_": 79228.11902553569
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 586.0192959847301,
                    "monoMz_": 586.0192959847301,
                    "intensity_": 79112.70425542643
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 586.1104226210938,
                    "monoMz_": 586.1104226210938,
                    "intensity_": 64911.14473278057
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 585.8370173483664,
                    "monoMz_": 585.8370173483664,
                    "intensity_": 61286.985154353824
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "5",
                    "pos_": 586.2015430756393,
                    "monoMz_": 586.2015430756393,
                    "intensity_": 45479.857500000006
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 585.7458607120029,
                    "monoMz_": 585.7458607120029,
                    "intensity_": 32655.24465652082
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "6",
                    "pos_": 586.2926587120029,
                    "monoMz_": 586.2926587120029,
                    "intensity_": 27932.507979686
                }
            ],
            "monoMass_": 6431.121533696363,
            "charge_": 11,
            "intensity_": -1
        }
    },
    {
        "mz": 605.4883422851562,
        "intensity": 47023.259765625,
        "text": "Y59",
        "error": 0.0227,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 1,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 605.4883422851562,
                    "monoMz_": 605.4883422851562,
                    "intensity_": 47023.259765625
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 605.3972063760654,
                    "monoMz_": 605.3972063760654,
                    "intensity_": 45798.90104412949
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 605.5794706487926,
                    "monoMz_": 605.5794706487926,
                    "intensity_": 39642.042191844725
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 605.3060610124289,
                    "monoMz_": 605.3060610124289,
                    "intensity_": 34423.81368727641
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 605.6705928306109,
                    "monoMz_": 605.6705928306109,
                    "intensity_": 28520.94015726592
                }
            ],
            "monoMass_": 6645.281010001049,
            "charge_": 11,
            "intensity_": -1
        }
    },
    {
        "mz": 614.3359375,
        "intensity": 311983.38,
        "text": "Y27",
        "error": -0.0231,
        "env": {
            "displayColor_": "darkorange",
            "displayLevel_": 1,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 614.3359375,
                    "monoMz_": 614.3359375,
                    "intensity_": 311983.38
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 614.5364447,
                    "monoMz_": 614.5364447,
                    "intensity_": 288625.77710782195
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 614.7369080999999,
                    "monoMz_": 614.7369080999999,
                    "intensity_": 189163.3353578352
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 614.1353579,
                    "monoMz_": 614.1353579,
                    "intensity_": 182068.867614807
                }
            ],
            "monoMass_": 3065.640407165605,
            "charge_": 5,
            "intensity_": -1
        }
    },
    {
        "mz": 627.2124633789062,
        "intensity": 98483.19359375,
        "text": "Y55",
        "error": 0.0019,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 0,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 627.2124633789062,
                    "monoMz_": 627.2124633789062,
                    "intensity_": 98483.19359375
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 627.3127093789062,
                    "monoMz_": 627.3127093789062,
                    "intensity_": 95982.31236046225
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 627.1122062789062,
                    "monoMz_": 627.1122062789062,
                    "intensity_": 78127.38716363831
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 627.4129466789062,
                    "monoMz_": 627.4129466789062,
                    "intensity_": 76924.23472964638
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "5",
                    "pos_": 627.5131770789062,
                    "monoMz_": 627.5131770789062,
                    "intensity_": 52679.07880686964
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 627.0119349789063,
                    "monoMz_": 627.0119349789063,
                    "intensity_": 42742.14623956286
                }
            ],
            "monoMass_": 6259.043693120271,
            "charge_": 10,
            "intensity_": -1
        }
    },
    {
        "mz": 632.5548095703125,
        "intensity": 39890.74685346494,
        "text": "Y51",
        "error": -0.0175,
        "env": {
            "displayColor_": "darkorange",
            "displayLevel_": 3,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 632.5548095703125,
                    "monoMz_": 632.5548095703125,
                    "intensity_": 39890.74685346494
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 632.6661866814236,
                    "monoMz_": 632.6661866814236,
                    "intensity_": 35718.885585701435
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 632.4434187925348,
                    "monoMz_": 632.4434187925348,
                    "intensity_": 34553.793750000004
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 632.7775533480904,
                    "monoMz_": 632.7775533480904,
                    "intensity_": 26367.1019372766
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 632.3320102369793,
                    "monoMz_": 632.3320102369793,
                    "intensity_": 20724.601270305393
                }
            ],
            "monoMass_": 5680.919714930902,
            "charge_": 9,
            "intensity_": -1
        }
    },
    {
        "mz": 637.3155517578125,
        "intensity": 80980.103515625,
        "text": "Y56",
        "error": -0.0149,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 2,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 637.3155517578125,
                    "monoMz_": 637.3155517578125,
                    "intensity_": 80980.103515625
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 637.4157990578125,
                    "monoMz_": 637.4157990578125,
                    "intensity_": 80098.0285478821
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 637.5160377578125,
                    "monoMz_": 637.5160377578125,
                    "intensity_": 65115.26876144992
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 637.2152935578125,
                    "monoMz_": 637.2152935578125,
                    "intensity_": 63259.86882572062
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "5",
                    "pos_": 637.6162695578125,
                    "monoMz_": 637.6162695578125,
                    "intensity_": 45212.7895302158
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 637.1150215578125,
                    "monoMz_": 637.1150215578125,
                    "intensity_": 34052.23475344971
                }
            ],
            "monoMass_": 6360.074556909335,
            "charge_": 10,
            "intensity_": -1
        }
    },
    {
        "mz": 643.7481689453125,
        "intensity": 57383.539813100186,
        "text": "Y28",
        "error": -0.0304,
        "env": {
            "displayColor_": "darkorange",
            "displayLevel_": 2,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 643.7481689453125,
                    "monoMz_": 643.7481689453125,
                    "intensity_": 57383.539813100186
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 643.9486817453125,
                    "monoMz_": 643.9486817453125,
                    "intensity_": 55618.071515020965
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 644.1491531453124,
                    "monoMz_": 644.1491531453124,
                    "intensity_": 38099.3008744407
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 643.5475893453124,
                    "monoMz_": 643.5475893453124,
                    "intensity_": 31841.157608137244
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 644.3495993453125,
                    "monoMz_": 644.3495993453125,
                    "intensity_": 20539.704140625003
                }
            ],
            "monoMass_": 3212.701564392167,
            "charge_": 5,
            "intensity_": -1
        }
    },
    {
        "mz": 649.3348999023438,
        "intensity": 170417.0790625,
        "text": "B11",
        "error": -0.0096,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 4,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 649.3348999023438,
                    "monoMz_": 649.3348999023438,
                    "intensity_": 170417.0790625
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 649.8363574023438,
                    "monoMz_": 649.8363574023438,
                    "intensity_": 126038.76920800518
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 650.3377219023438,
                    "monoMz_": 650.3377219023438,
                    "intensity_": 51872.482706466224
                }
            ],
            "monoMass_": 1296.6552468709294,
            "charge_": 2,
            "intensity_": -1
        }
    },
    {
        "mz": 663.1588745117188,
        "intensity": 100302.04121510527,
        "text": "Y29",
        "error": -0.0296,
        "env": {
            "displayColor_": "blue",
            "displayLevel_": 4,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 663.1588745117188,
                    "monoMz_": 663.1588745117188,
                    "intensity_": 100302.04121510527
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 663.3593887117188,
                    "monoMz_": 663.3593887117188,
                    "intensity_": 99793.71147324753
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 663.5598629117187,
                    "monoMz_": 663.5598629117187,
                    "intensity_": 70038.9003125
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 662.9582969117188,
                    "monoMz_": 662.9582969117188,
                    "intensity_": 54054.135121601175
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 663.7603123117187,
                    "monoMz_": 663.7603123117187,
                    "intensity_": 38636.968148714084
                }
            ],
            "monoMass_": 3309.7551022241987,
            "charge_": 5,
            "intensity_": -1
        }
    },
    {
        "mz": 679.7057495117188,
        "intensity": 107270.38210067827,
        "text": "Y18",
        "error": -0.015,
        "env": {
            "displayColor_": "blue",
            "displayLevel_": 5,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 679.7057495117188,
                    "monoMz_": 679.7057495117188,
                    "intensity_": 107270.38210067827
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 679.3714561783854,
                    "monoMz_": 679.3714561783854,
                    "intensity_": 94193.26221414114
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 680.039805845052,
                    "monoMz_": 680.039805845052,
                    "intensity_": 69895.08539062501
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 680.3737721783854,
                    "monoMz_": 680.3737721783854,
                    "intensity_": 33193.52498754293
                }
            ],
            "monoMass_": 2035.0925391345193,
            "charge_": 3,
            "intensity_": -1
        }
    },
    {
        "mz": 695.3983154296875,
        "intensity": 68693.47237043973,
        "text": "B6",
        "error": -0.0038,
        "env": {
            "displayColor_": "darkorange",
            "displayLevel_": 0,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 695.3983154296875,
                    "monoMz_": 695.3983154296875,
                    "intensity_": 68693.47237043973
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 696.4011974296875,
                    "monoMz_": 696.4011974296875,
                    "intensity_": 26847.116035156254
                }
            ],
            "monoMass_": 694.3910389628085,
            "charge_": 1,
            "intensity_": -1
        }
    },
    {
        "mz": 705.8740234375,
        "intensity": 71870.9345703125,
        "text": "B12",
        "error": -0.0155,
        "env": {
            "displayColor_": "blue",
            "displayLevel_": 6,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 705.8740234375,
                    "monoMz_": 705.8740234375,
                    "intensity_": 71870.9345703125
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 706.3754824375001,
                    "monoMz_": 706.3754824375001,
                    "intensity_": 56766.689889978225
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 706.8768499375001,
                    "monoMz_": 706.8768499375001,
                    "intensity_": 24913.2192975737
                }
            ],
            "monoMass_": 1409.733493941242,
            "charge_": 2,
            "intensity_": -1
        }
    }
];
