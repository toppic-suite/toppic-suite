"use strict";
let sampleMatchedIon = [
    {
        "mz": 365.7658996582031,
        "intensity": 15429.819365234374,
        "text": "Y6",
        "error": -0.0053,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 6,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 365.7658996582031,
                    "monoMz_": 365.7658996582031,
                    "intensity_": 15429.819365234374
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 366.2673491582031,
                    "monoMz_": 366.2673491582031,
                    "intensity_": 6213.491359114269
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 366.7686341582031,
                    "monoMz_": 366.7686341582031,
                    "intensity_": 1534.0841768882815
                }
            ],
            "monoMass_": 729.5172463826483,
            "charge_": 2,
            "intensity_": -1
        }
    },
    {
        "mz": 413.9083251953125,
        "intensity": 5740.253556604234,
        "text": "B11",
        "error": -0.0104,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 5,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 413.9083251953125,
                    "monoMz_": 413.9083251953125,
                    "intensity_": 5740.253556604234
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 414.24262952864575,
                    "monoMz_": 414.24262952864575,
                    "intensity_": 3959.9875634765626
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 414.57686719531245,
                    "monoMz_": 414.57686719531245,
                    "intensity_": 1532.9963052118233
                }
            ],
            "monoMass_": 1238.7031461853005,
            "charge_": 3,
            "intensity_": -1
        }
    },
    {
        "mz": 466.92822265625,
        "intensity": 2035.4222119140622,
        "text": "Y25",
        "error": -0.024,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 4,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 466.92822265625,
                    "monoMz_": 466.92822265625,
                    "intensity_": 2035.4222119140622
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 467.09530132291667,
                    "monoMz_": 467.09530132291667,
                    "intensity_": 1737.517674498765
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 466.76107382291667,
                    "monoMz_": 466.76107382291667,
                    "intensity_": 1301.1176362701958
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 467.2623408229167,
                    "monoMz_": 467.2623408229167,
                    "intensity_": 1057.3091256076027
                }
            ],
            "monoMass_": 2794.522784136226,
            "charge_": 6,
            "intensity_": -1
        }
    },
    {
        "mz": 483.43951416015625,
        "intensity": 3341.1076171875,
        "text": "Y26",
        "error": -0.0246,
        "env": {
            "displayColor_": "darkorange",
            "displayLevel_": 4,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 483.43951416015625,
                    "monoMz_": 483.43951416015625,
                    "intensity_": 3341.1076171875
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 483.6065963268229,
                    "monoMz_": 483.6065963268229,
                    "intensity_": 2954.522722265174
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 483.27236649348964,
                    "monoMz_": 483.27236649348964,
                    "intensity_": 2050.687659517734
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 483.7736408268229,
                    "monoMz_": 483.7736408268229,
                    "intensity_": 1856.722563565541
                }
            ],
            "monoMass_": 2893.5905401596638,
            "charge_": 6,
            "intensity_": -1
        }
    },
    {
        "mz": 497.2876892089844,
        "intensity": 1377.0401953125,
        "text": "B13",
        "error": -0.0153,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 2,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 497.2876892089844,
                    "monoMz_": 497.2876892089844,
                    "intensity_": 1377.0401953125
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 497.62198020898444,
                    "monoMz_": 497.62198020898444,
                    "intensity_": 1135.231344888341
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 497.9558722089844,
                    "monoMz_": 497.9558722089844,
                    "intensity_": 576.5303208345797
                }
            ],
            "monoMass_": 1488.8412382263161,
            "charge_": 3,
            "intensity_": -1
        }
    },
    {
        "mz": 499.7845153808594,
        "intensity": 2474.93408203125,
        "text": "B43",
        "error": -0.0421,
        "env": {
            "displayColor_": "darkorange",
            "displayLevel_": 3,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 499.7845153808594,
                    "monoMz_": 499.7845153808594,
                    "intensity_": 2474.93408203125
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 499.6842723808594,
                    "monoMz_": 499.6842723808594,
                    "intensity_": 2414.5169564449952
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 499.88474448085935,
                    "monoMz_": 499.88474448085935,
                    "intensity_": 1977.5939498037599
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 499.5840103808594,
                    "monoMz_": 499.5840103808594,
                    "intensity_": 1641.943736089453
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 499.9849635808594,
                    "monoMz_": 499.9849635808594,
                    "intensity_": 1307.785115647705
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "5",
                    "pos_": 500.0851754808594,
                    "monoMz_": 500.0851754808594,
                    "intensity_": 742.733657859375
                }
            ],
            "monoMass_": 4984.764451139804,
            "charge_": 10,
            "intensity_": -1
        }
    },
    {
        "mz": 505.2793884277344,
        "intensity": 4054.833323129341,
        "text": "Y27",
        "error": -0.0259,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 4,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 505.2793884277344,
                    "monoMz_": 505.2793884277344,
                    "intensity_": 4054.833323129341
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 505.44647642773435,
                    "monoMz_": 505.44647642773435,
                    "intensity_": 3730.2532417294806
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 505.61352809440103,
                    "monoMz_": 505.61352809440103,
                    "intensity_": 2432.6088162966703
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 505.1122392610677,
                    "monoMz_": 505.1122392610677,
                    "intensity_": 2382.1104492187496
                }
            ],
            "monoMass_": 3024.629776765132,
            "charge_": 6,
            "intensity_": -1
        }
    },
    {
        "mz": 543.863037109375,
        "intensity": 2865.5688318815555,
        "text": "Y34",
        "error": -0.0304,
        "env": {
            "displayColor_": "blue",
            "displayLevel_": 0,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 543.863037109375,
                    "monoMz_": 543.863037109375,
                    "intensity_": 2865.5688318815555
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 543.7198029665178,
                    "monoMz_": 543.7198029665178,
                    "intensity_": 2537.6189501953127
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 544.0062466808035,
                    "monoMz_": 544.0062466808035,
                    "intensity_": 2264.6115940161376
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 544.1494395379465,
                    "monoMz_": 544.1494395379465,
                    "intensity_": 1399.0901119778
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 543.5765328236606,
                    "monoMz_": 543.5765328236606,
                    "intensity_": 1190.5821826055244
                }
            ],
            "monoMass_": 3797.9847944974713,
            "charge_": 7,
            "intensity_": -1
        }
    },
    {
        "mz": 555.2042236328125,
        "intensity": 31799.327257465357,
        "text": "B43",
        "error": -0.042,
        "env": {
            "displayColor_": "darkorange",
            "displayLevel_": 6,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 555.2042236328125,
                    "monoMz_": 555.2042236328125,
                    "intensity_": 31799.327257465357
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 555.0928425217014,
                    "monoMz_": 555.0928425217014,
                    "intensity_": 31023.054482193744
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 555.3155892994791,
                    "monoMz_": 555.3155892994791,
                    "intensity_": 25409.225097656254
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 554.9814402994792,
                    "monoMz_": 554.9814402994792,
                    "intensity_": 21096.604786904616
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 555.4269438550348,
                    "monoMz_": 555.4269438550348,
                    "intensity_": 16803.149294704508
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "5",
                    "pos_": 555.5382904105903,
                    "monoMz_": 555.5382904105903,
                    "intensity_": 9543.054428350772
                }
            ],
            "monoMass_": 4984.764586493402,
            "charge_": 9,
            "intensity_": -1
        }
    },
    {
        "mz": 560.1124267578125,
        "intensity": 11648.381240234376,
        "text": "Y25",
        "error": -0.0239,
        "env": {
            "displayColor_": "darkorange",
            "displayLevel_": 0,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 560.1124267578125,
                    "monoMz_": 560.1124267578125,
                    "intensity_": 11648.381240234376
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 560.3129211578125,
                    "monoMz_": 560.3129211578125,
                    "intensity_": 9943.523346526985
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 559.9118481578125,
                    "monoMz_": 559.9118481578125,
                    "intensity_": 7446.078841507574
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 560.5133685578126,
                    "monoMz_": 560.5133685578126,
                    "intensity_": 6050.803470536265
                }
            ],
            "monoMass_": 2794.5228584546676,
            "charge_": 5,
            "intensity_": -1
        }
    },
    {
        "mz": 560.629638671875,
        "intensity": 23117.07370864125,
        "text": "Y30",
        "error": -0.0285,
        "env": {
            "displayColor_": "blue",
            "displayLevel_": 3,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 560.629638671875,
                    "monoMz_": 560.629638671875,
                    "intensity_": 23117.07370864125
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 560.4625413385417,
                    "monoMz_": 560.4625413385417,
                    "intensity_": 22991.35797605744
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 560.7967028385416,
                    "monoMz_": 560.7967028385416,
                    "intensity_": 16378.352404911595
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 560.2953923385417,
                    "monoMz_": 560.2953923385417,
                    "intensity_": 12239.713563878193
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 560.9637463385417,
                    "monoMz_": 560.9637463385417,
                    "intensity_": 9113.990224609373
                }
            ],
            "monoMass_": 3355.728695229976,
            "charge_": 6,
            "intensity_": -1
        }
    },
    {
        "mz": 576.9766235351562,
        "intensity": 49042.98085248435,
        "text": "Y47",
        "error": -0.0522,
        "env": {
            "displayColor_": "blue",
            "displayLevel_": 7,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 576.9766235351562,
                    "monoMz_": 576.9766235351562,
                    "intensity_": 49042.98085248435
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 576.8652390907117,
                    "monoMz_": 576.8652390907117,
                    "intensity_": 46258.03300781249
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 577.0879929796006,
                    "monoMz_": 577.0879929796006,
                    "intensity_": 40466.291394582564
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 576.7538344240452,
                    "monoMz_": 576.7538344240452,
                    "intensity_": 30346.21099194635
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 577.1993515351562,
                    "monoMz_": 577.1993515351562,
                    "intensity_": 27598.828599318076
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "5",
                    "pos_": 577.3107020907117,
                    "monoMz_": 577.3107020907117,
                    "intensity_": 16149.592195635152
                }
            ],
            "monoMass_": 5180.716129614495,
            "charge_": 9,
            "intensity_": -1
        }
    },
    {
        "mz": 579.9261474609375,
        "intensity": 12210.53853515625,
        "text": "Y26",
        "error": -0.0237,
        "env": {
            "displayColor_": "darkorange",
            "displayLevel_": 1,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 579.9261474609375,
                    "monoMz_": 579.9261474609375,
                    "intensity_": 12210.53853515625
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 580.1266460609374,
                    "monoMz_": 580.1266460609374,
                    "intensity_": 10797.71072551749
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 579.7255702609375,
                    "monoMz_": 579.7255702609375,
                    "intensity_": 7494.52084730782
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 580.3270994609375,
                    "monoMz_": 580.3270994609375,
                    "intensity_": 6785.648655817849
                }
            ],
            "monoMass_": 2893.591468970292,
            "charge_": 5,
            "intensity_": -1
        }
    },
    {
        "mz": 581.5867919921875,
        "intensity": 15372.449792195068,
        "text": "B40",
        "error": -0.0397,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 4,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 581.5867919921875,
                    "monoMz_": 581.5867919921875,
                    "intensity_": 15372.449792195068
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 581.7120889921874,
                    "monoMz_": 581.7120889921874,
                    "intensity_": 14821.6712109375
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 581.8373678671875,
                    "monoMz_": 581.8373678671875,
                    "intensity_": 11168.757011259131
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 581.4614689921874,
                    "monoMz_": 581.4614689921874,
                    "intensity_": 11152.87296633785
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "5",
                    "pos_": 581.9626338671875,
                    "monoMz_": 581.9626338671875,
                    "intensity_": 6978.463626184558
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 581.3361083671874,
                    "monoMz_": 581.3361083671874,
                    "intensity_": 4284.263787133778
                }
            ],
            "monoMass_": 4642.630655202467,
            "charge_": 8,
            "intensity_": -1
        }
    },
    {
        "mz": 589.5415649414062,
        "intensity": 5139.471598570317,
        "text": "Y48",
        "error": -0.0518,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 3,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 589.5415649414062,
                    "monoMz_": 589.5415649414062,
                    "intensity_": 5139.471598570317
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 589.4301793858508,
                    "monoMz_": 589.4301793858508,
                    "intensity_": 4751.717225529521
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 589.6529360525174,
                    "monoMz_": 589.6529360525174,
                    "intensity_": 4323.252172851563
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 589.3187741636285,
                    "monoMz_": 589.3187741636285,
                    "intensity_": 3052.698729505321
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 589.7642962747396,
                    "monoMz_": 589.7642962747396,
                    "intensity_": 3004.3839446644934
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "5",
                    "pos_": 589.8756484969618,
                    "monoMz_": 589.8756484969618,
                    "intensity_": 1790.6117946969848
                }
            ],
            "monoMass_": 5293.800591270746,
            "charge_": 9,
            "intensity_": -1
        }
    },
    {
        "mz": 593.7188720703125,
        "intensity": 3346.368110230592,
        "text": "B41",
        "error": -0.0359,
        "env": {
            "displayColor_": "darkorange",
            "displayLevel_": 3,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 593.7188720703125,
                    "monoMz_": 593.7188720703125,
                    "intensity_": 3346.368110230592
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 593.8441709453125,
                    "monoMz_": 593.8441709453125,
                    "intensity_": 3282.6902722431005
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 593.9694519453126,
                    "monoMz_": 593.9694519453126,
                    "intensity_": 2515.4608593549224
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 593.5935478203126,
                    "monoMz_": 593.5935478203126,
                    "intensity_": 2384.500341796875
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "5",
                    "pos_": 594.0947203203125,
                    "monoMz_": 594.0947203203125,
                    "intensity_": 1597.7099349024306
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 593.4681869453126,
                    "monoMz_": 593.4681869453126,
                    "intensity_": 898.6854940996695
                }
            ],
            "monoMass_": 4739.6872838274685,
            "charge_": 8,
            "intensity_": -1
        }
    },
    {
        "mz": 600.8199462890625,
        "intensity": 1825.3150678803663,
        "text": "Y43",
        "error": -0.0352,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 2,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 600.8199462890625,
                    "monoMz_": 600.8199462890625,
                    "intensity_": 1825.3150678803663
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 600.9452464140625,
                    "monoMz_": 600.9452464140625,
                    "intensity_": 1811.3325521064305
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 601.0705287890625,
                    "monoMz_": 601.0705287890625,
                    "intensity_": 1402.9167541503905
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 600.6946211640625,
                    "monoMz_": 600.6946211640625,
                    "intensity_": 1284.3191892586988
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 601.1957985390625,
                    "monoMz_": 601.1957985390625,
                    "intensity_": 900.0855851444037
                }
            ],
            "monoMass_": 4796.495871577467,
            "charge_": 8,
            "intensity_": -1
        }
    },
    {
        "mz": 601.6528930664062,
        "intensity": 3605.2427633611314,
        "text": "Y32",
        "error": -0.0258,
        "env": {
            "displayColor_": "darkorange",
            "displayLevel_": 4,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 601.6528930664062,
                    "monoMz_": 601.6528930664062,
                    "intensity_": 3605.2427633611314
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 601.4857890664063,
                    "monoMz_": 601.4857890664063,
                    "intensity_": 3341.0501049804684
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 601.8199665664063,
                    "monoMz_": 601.8199665664063,
                    "intensity_": 2730.4640995368245
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 601.318639233073,
                    "monoMz_": 601.318639233073,
                    "intensity_": 1647.355581321004
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 601.9870200664062,
                    "monoMz_": 601.9870200664062,
                    "intensity_": 1619.8983048022387
                }
            ],
            "monoMass_": 3601.8681765971637,
            "charge_": 6,
            "intensity_": -1
        }
    },
    {
        "mz": 606.1334838867188,
        "intensity": 6098.743366517672,
        "text": "Y27",
        "error": -0.0275,
        "env": {
            "displayColor_": "darkorange",
            "displayLevel_": 2,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 606.1334838867188,
                    "monoMz_": 606.1334838867188,
                    "intensity_": 6098.743366517672
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 606.3339894867188,
                    "monoMz_": 606.3339894867188,
                    "intensity_": 5610.552987137676
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 606.5344514867187,
                    "monoMz_": 606.5344514867187,
                    "intensity_": 3658.8080691494533
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 605.9329048867187,
                    "monoMz_": 605.9329048867187,
                    "intensity_": 3582.8551120994803
                }
            ],
            "monoMass_": 3024.628142099198,
            "charge_": 5,
            "intensity_": -1
        }
    },
    {
        "mz": 608.0985717773438,
        "intensity": 103948.04120254172,
        "text": "B42",
        "error": -0.0436,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 6,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 608.0985717773438,
                    "monoMz_": 608.0985717773438,
                    "intensity_": 103948.04120254172
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 607.9732697773437,
                    "monoMz_": 607.9732697773437,
                    "intensity_": 103575.78039842636
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 608.2238561523437,
                    "monoMz_": 608.2238561523437,
                    "intensity_": 81386.97539138954
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 607.8479431523438,
                    "monoMz_": 607.8479431523438,
                    "intensity_": 72016.0584375
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 608.3491279023438,
                    "monoMz_": 608.3491279023438,
                    "intensity_": 52767.91960803362
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "5",
                    "pos_": 608.4743905273438,
                    "monoMz_": 608.4743905273438,
                    "intensity_": 29394.641224219624
                }
            ],
            "monoMass_": 4853.722441483717,
            "charge_": 8,
            "intensity_": -1
        }
    },
    {
        "mz": 615.079833984375,
        "intensity": 1336.7757145749995,
        "text": "Y44",
        "error": -0.0426,
        "env": {
            "displayColor_": "blue",
            "displayLevel_": 2,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 615.079833984375,
                    "monoMz_": 615.079833984375,
                    "intensity_": 1336.7757145749995
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 614.954531359375,
                    "monoMz_": 614.954531359375,
                    "intensity_": 1320.1367202103559
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 615.205118984375,
                    "monoMz_": 615.205118984375,
                    "intensity_": 1055.7429766391253
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 614.8292044843749,
                    "monoMz_": 614.8292044843749,
                    "intensity_": 909.3696719231067
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 615.3303916093749,
                    "monoMz_": 615.3303916093749,
                    "intensity_": 690.3156978248021
                }
            ],
            "monoMass_": 4909.572536139967,
            "charge_": 8,
            "intensity_": -1
        }
    },
    {
        "mz": 620.3587646484375,
        "intensity": 20765.290613552115,
        "text": "B11",
        "error": -0.0106,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 4,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 620.3587646484375,
                    "monoMz_": 620.3587646484375,
                    "intensity_": 20765.290613552115
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 620.8602211484374,
                    "monoMz_": 620.8602211484374,
                    "intensity_": 14325.2021484375
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 621.3615776484374,
                    "monoMz_": 621.3615776484374,
                    "intensity_": 5545.593669917376
                }
            ],
            "monoMass_": 1238.702976363117,
            "charge_": 2,
            "intensity_": -1
        }
    },
    {
        "mz": 627.7102661132812,
        "intensity": 12226.31466825993,
        "text": "Y45",
        "error": -0.0468,
        "env": {
            "displayColor_": "blue",
            "displayLevel_": 3,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 627.7102661132812,
                    "monoMz_": 627.7102661132812,
                    "intensity_": 12226.31466825993
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 627.5849612382813,
                    "monoMz_": 627.5849612382813,
                    "intensity_": 11848.959124812674
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 627.8355537382813,
                    "monoMz_": 627.8355537382813,
                    "intensity_": 9829.952958597145
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 627.4596327382812,
                    "monoMz_": 627.4596327382812,
                    "intensity_": 7999.462871093751
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 627.9608288632813,
                    "monoMz_": 627.9608288632813,
                    "intensity_": 6538.498839650353
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "5",
                    "pos_": 628.0860948632813,
                    "monoMz_": 628.0860948632813,
                    "intensity_": 3734.036195307185
                }
            ],
            "monoMass_": 5010.615958171218,
            "charge_": 8,
            "intensity_": -1
        }
    },
    {
        "mz": 631.9417114257812,
        "intensity": 3552.3355913512937,
        "text": "Y28",
        "error": -0.0289,
        "env": {
            "displayColor_": "blue",
            "displayLevel_": 0,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 631.9417114257812,
                    "monoMz_": 631.9417114257812,
                    "intensity_": 3552.3355913512937
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 632.1422204257813,
                    "monoMz_": 632.1422204257813,
                    "intensity_": 3382.4284101942826
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 632.3426872257812,
                    "monoMz_": 632.3426872257812,
                    "intensity_": 2278.596323130924
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 631.7411340257812,
                    "monoMz_": 631.7411340257812,
                    "intensity_": 2009.8369140624993
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 632.5431288257812,
                    "monoMz_": 632.5431288257812,
                    "intensity_": 1208.836889172281
                }
            ],
            "monoMass_": 3153.669287794511,
            "charge_": 5,
            "intensity_": -1
        }
    },
    {
        "mz": 634.3386840820312,
        "intensity": 13102.891319201544,
        "text": "Y34",
        "error": -0.0323,
        "env": {
            "displayColor_": "blue",
            "displayLevel_": 5,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 634.3386840820312,
                    "monoMz_": 634.3386840820312,
                    "intensity_": 13102.891319201544
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 634.1715775820312,
                    "monoMz_": 634.1715775820312,
                    "intensity_": 11603.33157732009
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 634.5057619153646,
                    "monoMz_": 634.5057619153646,
                    "intensity_": 10354.998025684736
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 634.6728202486979,
                    "monoMz_": 634.6728202486979,
                    "intensity_": 6397.377539515465
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 634.0044290820312,
                    "monoMz_": 634.0044290820312,
                    "intensity_": 5443.969368907053
                }
            ],
            "monoMass_": 3797.982915690913,
            "charge_": 6,
            "intensity_": -1
        }
    },
    {
        "mz": 645.8632202148438,
        "intensity": 1223.3082301699733,
        "text": "B45",
        "error": -0.0725,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 1,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 645.8632202148438,
                    "monoMz_": 645.8632202148438,
                    "intensity_": 1223.3082301699733
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 645.7379133398437,
                    "monoMz_": 645.7379133398437,
                    "intensity_": 1157.6332763671876
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 645.9885103398437,
                    "monoMz_": 645.9885103398437,
                    "intensity_": 1006.2911726491703
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 645.6125835898437,
                    "monoMz_": 645.6125835898437,
                    "intensity_": 762.1479891092864
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 645.4872223398437,
                    "monoMz_": 645.4872223398437,
                    "intensity_": 263.64950608947083
                }
            ],
            "monoMass_": 5155.839566983717,
            "charge_": 8,
            "intensity_": -1
        }
    },
    {
        "mz": 648.8486328125,
        "intensity": 8900.361702115466,
        "text": "Y47",
        "error": -0.0455,
        "env": {
            "displayColor_": "blue",
            "displayLevel_": 4,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 648.8486328125,
                    "monoMz_": 648.8486328125,
                    "intensity_": 8900.361702115466
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 648.7233253124999,
                    "monoMz_": 648.7233253124999,
                    "intensity_": 8394.947008549776
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 648.9739234374999,
                    "monoMz_": 648.9739234374999,
                    "intensity_": 7343.8568352588745
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 649.0992018124999,
                    "monoMz_": 649.0992018124999,
                    "intensity_": 5008.6587890625
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 649.2244711875,
                    "monoMz_": 649.2244711875,
                    "intensity_": 2930.841669578751
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "5",
                    "pos_": 649.3497339375,
                    "monoMz_": 649.3497339375,
                    "intensity_": 1509.0183220492092
                }
            ],
            "monoMass_": 5179.722856764967,
            "charge_": 8,
            "intensity_": -1
        }
    },
    {
        "mz": 654.947021484375,
        "intensity": 3806.907514069733,
        "text": "Y29",
        "error": -0.0293,
        "env": {
            "displayColor_": "blue",
            "displayLevel_": 2,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 654.947021484375,
                    "monoMz_": 654.947021484375,
                    "intensity_": 3806.907514069733
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 655.147536484375,
                    "monoMz_": 655.147536484375,
                    "intensity_": 3756.1049804687495
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 655.348010684375,
                    "monoMz_": 655.348010684375,
                    "intensity_": 2614.3262388171647
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 654.746442284375,
                    "monoMz_": 654.746442284375,
                    "intensity_": 2068.9268928405545
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 655.548459884375,
                    "monoMz_": 655.548459884375,
                    "intensity_": 1430.2160180267542
                }
            ],
            "monoMass_": 3268.69582908748,
            "charge_": 5,
            "intensity_": -1
        }
    },
    {
        "mz": 657.89453125,
        "intensity": 3857.4065596962387,
        "text": "B23",
        "error": -0.0225,
        "env": {
            "displayColor_": "blue",
            "displayLevel_": 0,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 657.89453125,
                    "monoMz_": 657.89453125,
                    "intensity_": 3857.4065596962387
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 658.145136,
                    "monoMz_": 658.145136,
                    "intensity_": 3110.079594726563
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 657.643809,
                    "monoMz_": 657.643809,
                    "intensity_": 2633.4732526516846
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 658.3956792500001,
                    "monoMz_": 658.3956792500001,
                    "intensity_": 1795.956612076341
                }
            ],
            "monoMass_": 2626.546130132484,
            "charge_": 4,
            "intensity_": -1
        }
    },
    {
        "mz": 672.5538330078125,
        "intensity": 9980.919140625,
        "text": "Y30",
        "error": -0.0299,
        "env": {
            "displayColor_": "darkorange",
            "displayLevel_": 4,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 672.5538330078125,
                    "monoMz_": 672.5538330078125,
                    "intensity_": 9980.919140625
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 672.3533162078124,
                    "monoMz_": 672.3533162078124,
                    "intensity_": 9926.640706536069
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 672.7543100078125,
                    "monoMz_": 672.7543100078125,
                    "intensity_": 7071.440488982719
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 672.1527374078125,
                    "monoMz_": 672.1527374078125,
                    "intensity_": 5284.56122626862
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 672.9547622078125,
                    "monoMz_": 672.9547622078125,
                    "intensity_": 3935.013601927886
                }
            ],
            "monoMass_": 3355.7273047046674,
            "charge_": 5,
            "intensity_": -1
        }
    },
    {
        "mz": 686.5067138671875,
        "intensity": 3418.98515625,
        "text": "Y43",
        "error": -0.0405,
        "env": {
            "displayColor_": "blue",
            "displayLevel_": 2,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 686.5067138671875,
                    "monoMz_": 686.5067138671875,
                    "intensity_": 3418.98515625
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 686.6499140100447,
                    "monoMz_": 686.6499140100447,
                    "intensity_": 3392.7946016880237
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 686.7930938671874,
                    "monoMz_": 686.7930938671874,
                    "intensity_": 2627.793766839703
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 686.3634851529017,
                    "monoMz_": 686.3634851529017,
                    "intensity_": 2405.6494800438045
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 686.9362592957589,
                    "monoMz_": 686.9362592957589,
                    "intensity_": 1685.9441469120702
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "5",
                    "pos_": 687.079414438616,
                    "monoMz_": 687.079414438616,
                    "intensity_": 929.5724203199062
                }
            ],
            "monoMass_": 4796.490574802159,
            "charge_": 7,
            "intensity_": -1
        }
    },
    {
        "mz": 694.8258666992188,
        "intensity": 14161.8598046875,
        "text": "B42",
        "error": -0.0439,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 6,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 694.8258666992188,
                    "monoMz_": 694.8258666992188,
                    "intensity_": 14161.8598046875
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 694.6826644135044,
                    "monoMz_": 694.6826644135044,
                    "intensity_": 14111.143069117757
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 694.9690488420758,
                    "monoMz_": 694.9690488420758,
                    "intensity_": 11088.144827804867
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 694.539433984933,
                    "monoMz_": 694.539433984933,
                    "intensity_": 9811.453024793693
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 695.1122165563617,
                    "monoMz_": 695.1122165563617,
                    "intensity_": 7189.090540127658
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "5",
                    "pos_": 695.2553738420759,
                    "monoMz_": 695.2553738420759,
                    "intensity_": 4004.719889000729
                }
            ],
            "monoMass_": 4853.722210626378,
            "charge_": 7,
            "intensity_": -1
        }
    },
    {
        "mz": 710.3971828884548,
        "intensity": 1389.7427501195593,
        "text": "B56",
        "error": -0.072,
        "env": {
            "displayColor_": "red",
            "displayLevel_": 2,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 710.3971828884548,
                    "monoMz_": 710.3971828884548,
                    "intensity_": 1389.7427501195593
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 710.5085687773438,
                    "monoMz_": 710.5085687773438,
                    "intensity_": 1375.8054647962904
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 710.2857848884548,
                    "monoMz_": 710.2857848884548,
                    "intensity_": 1084.900231933594
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 710.1743715551215,
                    "monoMz_": 710.1743715551215,
                    "intensity_": 583.7338696915585
                }
            ],
            "monoMass_": 6381.500961794184,
            "charge_": 9,
            "intensity_": -1
        }
    },
    {
        "mz": 713.5457153320312,
        "intensity": 8139.907200715501,
        "text": "B43",
        "error": -0.0454,
        "env": {
            "displayColor_": "darkorange",
            "displayLevel_": 4,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 713.5457153320312,
                    "monoMz_": 713.5457153320312,
                    "intensity_": 8139.907200715501
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 713.402511046317,
                    "monoMz_": 713.402511046317,
                    "intensity_": 7941.198960695443
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 713.6888997606027,
                    "monoMz_": 713.6888997606027,
                    "intensity_": 6504.185848411529
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 713.2592796177456,
                    "monoMz_": 713.2592796177456,
                    "intensity_": 5400.252773437501
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "4",
                    "pos_": 713.8320699034599,
                    "monoMz_": 713.8320699034599,
                    "intensity_": 4301.225457735199
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "5",
                    "pos_": 713.9752297606027,
                    "monoMz_": 713.9752297606027,
                    "intensity_": 2442.8056867120035
                }
            ],
            "monoMass_": 4984.761134056066,
            "charge_": 7,
            "intensity_": -1
        }
    },
    {
        "mz": 722.663330078125,
        "intensity": 3573.3419243234366,
        "text": "B26",
        "error": -0.0278,
        "env": {
            "displayColor_": "blue",
            "displayLevel_": 2,
            "peaks_": [
                {
                    "displayLevel_": -1,
                    "peakId_": "1",
                    "pos_": 722.663330078125,
                    "monoMz_": 722.663330078125,
                    "intensity_": 3573.3419243234366
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "2",
                    "pos_": 722.913955828125,
                    "monoMz_": 722.913955828125,
                    "intensity_": 3153.6504497057213
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "0",
                    "pos_": 722.4126055781251,
                    "monoMz_": 722.4126055781251,
                    "intensity_": 2198.3903466796874
                },
                {
                    "displayLevel_": -1,
                    "peakId_": "3",
                    "pos_": 723.164524328125,
                    "monoMz_": 723.164524328125,
                    "intensity_": 1978.3853660127531
                }
            ],
            "monoMass_": 2885.621316444984,
            "charge_": 4,
            "intensity_": -1
        }
    }
];
const matchedIonMassGraph = [
    {
        "mz": 729.5172463826483,
        "intensity": 23177.394901236927,
        "text": "Y6",
        "pos": "6",
        "error": -0.0053
    },
    {
        "mz": 1238.702976363117,
        "intensity": 40636.086431906995,
        "text": "B11",
        "pos": "11",
        "error": -0.0106
    },
    {
        "mz": 1488.8412382263161,
        "intensity": 3088.80186103542,
        "text": "B13",
        "pos": "13",
        "error": -0.0153
    },
    {
        "mz": 2626.546130132484,
        "intensity": 11396.916019150827,
        "text": "B23",
        "pos": "23",
        "error": -0.0225
    },
    {
        "mz": 2794.5228584546676,
        "intensity": 35088.7868988052,
        "text": "Y25",
        "pos": "25",
        "error": -0.0239
    },
    {
        "mz": 2885.621316444984,
        "intensity": 10903.768086721599,
        "text": "B26",
        "pos": "26",
        "error": -0.0278
    },
    {
        "mz": 2893.591468970292,
        "intensity": 37288.41876379941,
        "text": "Y26",
        "pos": "26",
        "error": -0.0237
    },
    {
        "mz": 3024.628142099198,
        "intensity": 18950.959534904283,
        "text": "Y27",
        "pos": "27",
        "error": -0.0275
    },
    {
        "mz": 3153.669287794511,
        "intensity": 12432.03412791128,
        "text": "Y28",
        "pos": "28",
        "error": -0.0289
    },
    {
        "mz": 3268.69582908748,
        "intensity": 13676.481644222957,
        "text": "Y29",
        "pos": "29",
        "error": -0.0293
    },
    {
        "mz": 3355.728695229976,
        "intensity": 83840.48787809785,
        "text": "Y30",
        "pos": "30",
        "error": -0.0285
    },
    {
        "mz": 3601.8681765971637,
        "intensity": 12944.010854001666,
        "text": "Y32",
        "pos": "32",
        "error": -0.0258
    },
    {
        "mz": 3797.982915690913,
        "intensity": 46902.56783062889,
        "text": "Y34",
        "pos": "34",
        "error": -0.0323
    },
    {
        "mz": 4642.630655202467,
        "intensity": 63778.47839404788,
        "text": "B40",
        "pos": "40",
        "error": -0.0397
    },
    {
        "mz": 4739.6872838274685,
        "intensity": 14025.415012627589,
        "text": "B41",
        "pos": "41",
        "error": -0.0359
    },
    {
        "mz": 4796.490574802159,
        "intensity": 14460.739572053508,
        "text": "Y43",
        "pos": "43",
        "error": -0.0405
    },
    {
        "mz": 4853.722441483717,
        "intensity": 443089.41626211087,
        "text": "B42",
        "pos": "42",
        "error": -0.0436
    },
    {
        "mz": 4909.572536139967,
        "intensity": 5312.34078117239,
        "text": "Y44",
        "pos": "44",
        "error": -0.0426
    },
    {
        "mz": 4984.764586493402,
        "intensity": 135674.41534727524,
        "text": "B43",
        "pos": "43",
        "error": -0.042
    },
    {
        "mz": 5010.615958171218,
        "intensity": 52177.224657721046,
        "text": "Y45",
        "pos": "45",
        "error": -0.0468
    },
    {
        "mz": 5155.839566983717,
        "intensity": 4413.030174385089,
        "text": "B45",
        "pos": "45",
        "error": -0.0725
    },
    {
        "mz": 5180.716129614495,
        "intensity": 209861.937041779,
        "text": "Y47",
        "pos": "47",
        "error": -0.0522
    },
    {
        "mz": 5293.800591270746,
        "intensity": 22062.135465818203,
        "text": "Y48",
        "pos": "48",
        "error": -0.0518
    },
    {
        "mz": 6381.500961794184,
        "intensity": 4434.182316541002,
        "text": "B56",
        "pos": "56",
        "error": -0.072
    }
];
