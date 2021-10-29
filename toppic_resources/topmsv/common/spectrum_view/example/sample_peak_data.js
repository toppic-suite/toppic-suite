"use strict";
let samplePeaks = [
    {
        "displayLevel_": 7,
        "peakId_": "505",
        "pos_": 608.0986,
        "monoMz_": 608.0986,
        "intensity_": 108110
    },
    {
        "displayLevel_": 0,
        "peakId_": "502",
        "pos_": 607.9734,
        "monoMz_": 607.9734,
        "intensity_": 105560
    },
    {
        "displayLevel_": 0,
        "peakId_": "507",
        "pos_": 608.2237,
        "monoMz_": 608.2237,
        "intensity_": 77823
    },
    {
        "displayLevel_": 0,
        "peakId_": "501",
        "pos_": 607.8482,
        "monoMz_": 607.8482,
        "intensity_": 76613
    },
    {
        "displayLevel_": 2,
        "peakId_": "509",
        "pos_": 608.3484,
        "monoMz_": 608.3484,
        "intensity_": 58436
    },
    {
        "displayLevel_": 6,
        "peakId_": "330",
        "pos_": 576.9766,
        "monoMz_": 576.9766,
        "intensity_": 52395
    },
    {
        "displayLevel_": 0,
        "peakId_": "329",
        "pos_": 576.866,
        "monoMz_": 576.866,
        "intensity_": 44055
    },
    {
        "displayLevel_": 7,
        "peakId_": "914",
        "pos_": 1877.0382,
        "monoMz_": 1877.0382,
        "intensity_": 42944
    },
    {
        "displayLevel_": 0,
        "peakId_": "333",
        "pos_": 577.0887,
        "monoMz_": 577.0887,
        "intensity_": 36573
    },
    {
        "displayLevel_": 7,
        "peakId_": "209",
        "pos_": 555.2042,
        "monoMz_": 555.2042,
        "intensity_": 32616
    },
    {
        "displayLevel_": 0,
        "peakId_": "327",
        "pos_": 576.7546,
        "monoMz_": 576.7546,
        "intensity_": 32268
    },
    {
        "displayLevel_": 0,
        "peakId_": "208",
        "pos_": 555.093,
        "monoMz_": 555.093,
        "intensity_": 30600
    },
    {
        "displayLevel_": 3,
        "peakId_": "500",
        "pos_": 607.7231,
        "monoMz_": 607.7231,
        "intensity_": 29555
    },
    {
        "displayLevel_": 0,
        "peakId_": "210",
        "pos_": 555.3156,
        "monoMz_": 555.3156,
        "intensity_": 26747
    },
    {
        "displayLevel_": 0,
        "peakId_": "336",
        "pos_": 577.1998,
        "monoMz_": 577.1998,
        "intensity_": 26738
    },
    {
        "displayLevel_": 4,
        "peakId_": "241",
        "pos_": 560.4628,
        "monoMz_": 560.4628,
        "intensity_": 24586
    },
    {
        "displayLevel_": 0,
        "peakId_": "243",
        "pos_": 560.6296,
        "monoMz_": 560.6296,
        "intensity_": 23809
    },
    {
        "displayLevel_": 2,
        "peakId_": "313",
        "pos_": 574.9761,
        "monoMz_": 574.9761,
        "intensity_": 23278
    },
    {
        "displayLevel_": 0,
        "peakId_": "511",
        "pos_": 608.4741,
        "monoMz_": 608.4741,
        "intensity_": 23208
    },
    {
        "displayLevel_": 4,
        "peakId_": "548",
        "pos_": 620.3588,
        "monoMz_": 620.3588,
        "intensity_": 20672
    },
    {
        "displayLevel_": 0,
        "peakId_": "207",
        "pos_": 554.9816,
        "monoMz_": 554.9816,
        "intensity_": 19026
    },
    {
        "displayLevel_": 2,
        "peakId_": "455",
        "pos_": 603.8451,
        "monoMz_": 603.8451,
        "intensity_": 18354
    },
    {
        "displayLevel_": 0,
        "peakId_": "337",
        "pos_": 577.3104,
        "monoMz_": 577.3104,
        "intensity_": 18142
    },
    {
        "displayLevel_": 4,
        "peakId_": "365",
        "pos_": 581.5868,
        "monoMz_": 581.5868,
        "intensity_": 16971
    },
    {
        "displayLevel_": 0,
        "peakId_": "211",
        "pos_": 555.427,
        "monoMz_": 555.427,
        "intensity_": 16758
    },
    {
        "displayLevel_": 0,
        "peakId_": "312",
        "pos_": 574.8649,
        "monoMz_": 574.8649,
        "intensity_": 16727
    },
    {
        "displayLevel_": 0,
        "peakId_": "367",
        "pos_": 581.7121,
        "monoMz_": 581.7121,
        "intensity_": 16288
    },
    {
        "displayLevel_": 4,
        "peakId_": "454",
        "pos_": 603.6782,
        "monoMz_": 603.6782,
        "intensity_": 15970
    },
    {
        "displayLevel_": 6,
        "peakId_": "60",
        "pos_": 365.7659,
        "monoMz_": 365.7659,
        "intensity_": 15907
    },
    {
        "displayLevel_": 0,
        "peakId_": "245",
        "pos_": 560.7965,
        "monoMz_": 560.7965,
        "intensity_": 14955
    },
    {
        "displayLevel_": 5,
        "peakId_": "790",
        "pos_": 694.8259,
        "monoMz_": 694.8259,
        "intensity_": 14600
    },
    {
        "displayLevel_": 2,
        "peakId_": "240",
        "pos_": 560.2985,
        "monoMz_": 560.2985,
        "intensity_": 14578
    },
    {
        "displayLevel_": 5,
        "peakId_": "549",
        "pos_": 620.8601,
        "monoMz_": 620.8601,
        "intensity_": 14325
    },
    {
        "displayLevel_": 0,
        "peakId_": "310",
        "pos_": 574.7534,
        "monoMz_": 574.7534,
        "intensity_": 14165
    },
    {
        "displayLevel_": 0,
        "peakId_": "325",
        "pos_": 576.6427,
        "monoMz_": 576.6427,
        "intensity_": 13922
    },
    {
        "displayLevel_": 4,
        "peakId_": "621",
        "pos_": 634.3387,
        "monoMz_": 634.3387,
        "intensity_": 13634
    },
    {
        "displayLevel_": 0,
        "peakId_": "789",
        "pos_": 694.6822,
        "monoMz_": 694.6822,
        "intensity_": 13378
    },
    {
        "displayLevel_": 2,
        "peakId_": "357",
        "pos_": 579.9261,
        "monoMz_": 579.9261,
        "intensity_": 13130
    },
    {
        "displayLevel_": 3,
        "peakId_": "595",
        "pos_": 627.7103,
        "monoMz_": 627.7103,
        "intensity_": 13129
    },
    {
        "displayLevel_": 0,
        "peakId_": "239",
        "pos_": 560.1124,
        "monoMz_": 560.1124,
        "intensity_": 12800
    },
    {
        "displayLevel_": 0,
        "peakId_": "315",
        "pos_": 575.0873,
        "monoMz_": 575.0873,
        "intensity_": 12492
    },
    {
        "displayLevel_": 0,
        "peakId_": "513",
        "pos_": 608.5999,
        "monoMz_": 608.5999,
        "intensity_": 12411
    },
    {
        "displayLevel_": 0,
        "peakId_": "458",
        "pos_": 604.0122,
        "monoMz_": 604.0122,
        "intensity_": 12298
    },
    {
        "displayLevel_": 4,
        "peakId_": "920",
        "pos_": 1904.9218,
        "monoMz_": 1904.9218,
        "intensity_": 11877
    },
    {
        "displayLevel_": 0,
        "peakId_": "316",
        "pos_": 575.1989,
        "monoMz_": 575.1989,
        "intensity_": 11756
    },
    {
        "displayLevel_": 0,
        "peakId_": "620",
        "pos_": 634.1708,
        "monoMz_": 634.1708,
        "intensity_": 11709
    },
    {
        "displayLevel_": 3,
        "peakId_": "782",
        "pos_": 692.3787,
        "monoMz_": 692.3787,
        "intensity_": 11585
    },
    {
        "displayLevel_": 2,
        "peakId_": "596",
        "pos_": 627.8353,
        "monoMz_": 627.8353,
        "intensity_": 11517
    },
    {
        "displayLevel_": 0,
        "peakId_": "784",
        "pos_": 692.8797,
        "monoMz_": 692.8797,
        "intensity_": 11071
    },
    {
        "displayLevel_": 0,
        "peakId_": "623",
        "pos_": 634.5053,
        "monoMz_": 634.5053,
        "intensity_": 11013
    },
    {
        "displayLevel_": 0,
        "peakId_": "792",
        "pos_": 694.9695,
        "monoMz_": 694.9695,
        "intensity_": 11004
    },
    {
        "displayLevel_": 0,
        "peakId_": "788",
        "pos_": 694.5395,
        "monoMz_": 694.5395,
        "intensity_": 10831
    },
    {
        "displayLevel_": 5,
        "peakId_": "739",
        "pos_": 672.3534,
        "monoMz_": 672.3534,
        "intensity_": 10688
    },
    {
        "displayLevel_": 2,
        "peakId_": "569",
        "pos_": 624.2277,
        "monoMz_": 624.2277,
        "intensity_": 10484
    },
    {
        "displayLevel_": 0,
        "peakId_": "570",
        "pos_": 624.3529,
        "monoMz_": 624.3529,
        "intensity_": 10399
    },
    {
        "displayLevel_": 0,
        "peakId_": "594",
        "pos_": 627.5851,
        "monoMz_": 627.5851,
        "intensity_": 10108
    },
    {
        "displayLevel_": 0,
        "peakId_": "358",
        "pos_": 580.1265,
        "monoMz_": 580.1265,
        "intensity_": 10077
    },
    {
        "displayLevel_": 6,
        "peakId_": "814",
        "pos_": 713.5457,
        "monoMz_": 713.5457,
        "intensity_": 9660.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "571",
        "pos_": 624.4787,
        "monoMz_": 624.4787,
        "intensity_": 9567.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "247",
        "pos_": 560.9639,
        "monoMz_": 560.9639,
        "intensity_": 9395.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "740",
        "pos_": 672.5538,
        "monoMz_": 672.5538,
        "intensity_": 9241.6
    },
    {
        "displayLevel_": 1,
        "peakId_": "368",
        "pos_": 581.8369,
        "monoMz_": 581.8369,
        "intensity_": 9212.2
    },
    {
        "displayLevel_": 4,
        "peakId_": "672",
        "pos_": 648.848,
        "monoMz_": 648.848,
        "intensity_": 9037.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "363",
        "pos_": 581.4614,
        "monoMz_": 581.4614,
        "intensity_": 8874.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "212",
        "pos_": 555.5384,
        "monoMz_": 555.5384,
        "intensity_": 8239
    },
    {
        "displayLevel_": 0,
        "peakId_": "206",
        "pos_": 554.8701,
        "monoMz_": 554.8701,
        "intensity_": 7989.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "572",
        "pos_": 624.6035,
        "monoMz_": 624.6035,
        "intensity_": 7980.4
    },
    {
        "displayLevel_": 1,
        "peakId_": "354",
        "pos_": 579.7252,
        "monoMz_": 579.7252,
        "intensity_": 7791.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "453",
        "pos_": 603.5102,
        "monoMz_": 603.5102,
        "intensity_": 7737.8
    },
    {
        "displayLevel_": 3,
        "peakId_": "650",
        "pos_": 642.8436,
        "monoMz_": 642.8436,
        "intensity_": 7676.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "238",
        "pos_": 559.9122,
        "monoMz_": 559.9122,
        "intensity_": 7556.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "593",
        "pos_": 627.46,
        "monoMz_": 627.46,
        "intensity_": 7546.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "673",
        "pos_": 648.9736,
        "monoMz_": 648.9736,
        "intensity_": 7487.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "508",
        "pos_": 608.2525,
        "monoMz_": 608.2525,
        "intensity_": 7372.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "741",
        "pos_": 672.754,
        "monoMz_": 672.754,
        "intensity_": 7160.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "317",
        "pos_": 575.3101,
        "monoMz_": 575.3101,
        "intensity_": 7036.7
    },
    {
        "displayLevel_": 4,
        "peakId_": "827",
        "pos_": 719.7056,
        "monoMz_": 719.7056,
        "intensity_": 6987.1
    },
    {
        "displayLevel_": 1,
        "peakId_": "483",
        "pos_": 606.1335,
        "monoMz_": 606.1335,
        "intensity_": 6820.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "828",
        "pos_": 720.0391,
        "monoMz_": 720.0391,
        "intensity_": 6815.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "793",
        "pos_": 695.1127,
        "monoMz_": 695.1127,
        "intensity_": 6782.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "813",
        "pos_": 713.4034,
        "monoMz_": 713.4034,
        "intensity_": 6714.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "625",
        "pos_": 634.6719,
        "monoMz_": 634.6719,
        "intensity_": 6634
    },
    {
        "displayLevel_": 0,
        "peakId_": "242",
        "pos_": 560.5135,
        "monoMz_": 560.5135,
        "intensity_": 6622.4
    },
    {
        "displayLevel_": 5,
        "peakId_": "151",
        "pos_": 505.2794,
        "monoMz_": 505.2794,
        "intensity_": 6273.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "671",
        "pos_": 648.7228,
        "monoMz_": 648.7228,
        "intensity_": 6213.2
    },
    {
        "displayLevel_": 2,
        "peakId_": "812",
        "pos_": 713.2589,
        "monoMz_": 713.2589,
        "intensity_": 6136.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "460",
        "pos_": 604.1797,
        "monoMz_": 604.1797,
        "intensity_": 6055.6
    },
    {
        "displayLevel_": 1,
        "peakId_": "225",
        "pos_": 557.4612,
        "monoMz_": 557.4612,
        "intensity_": 6052
    },
    {
        "displayLevel_": 0,
        "peakId_": "553",
        "pos_": 621.3618,
        "monoMz_": 621.3618,
        "intensity_": 5848.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "359",
        "pos_": 580.3283,
        "monoMz_": 580.3283,
        "intensity_": 5694.2
    },
    {
        "displayLevel_": 4,
        "peakId_": "701",
        "pos_": 657.8113,
        "monoMz_": 657.8113,
        "intensity_": 5679.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "226",
        "pos_": 557.6284,
        "monoMz_": 557.6284,
        "intensity_": 5568.8
    },
    {
        "displayLevel_": 4,
        "peakId_": "126",
        "pos_": 493.0328,
        "monoMz_": 493.0328,
        "intensity_": 5419.8
    },
    {
        "displayLevel_": 5,
        "peakId_": "58",
        "pos_": 356.7606,
        "monoMz_": 356.7606,
        "intensity_": 5410.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "506",
        "pos_": 608.1639,
        "monoMz_": 608.1639,
        "intensity_": 5353.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "338",
        "pos_": 577.4226,
        "monoMz_": 577.4226,
        "intensity_": 5329.8
    },
    {
        "displayLevel_": 2,
        "peakId_": "524",
        "pos_": 611.3502,
        "monoMz_": 611.3502,
        "intensity_": 5277.6
    },
    {
        "displayLevel_": 3,
        "peakId_": "274",
        "pos_": 568.496,
        "monoMz_": 568.496,
        "intensity_": 5276.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "675",
        "pos_": 649.2239,
        "monoMz_": 649.2239,
        "intensity_": 5276
    },
    {
        "displayLevel_": 0,
        "peakId_": "698",
        "pos_": 657.5604,
        "monoMz_": 657.5604,
        "intensity_": 5246.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "738",
        "pos_": 672.153,
        "monoMz_": 672.153,
        "intensity_": 5192.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "510",
        "pos_": 608.4352,
        "monoMz_": 608.4352,
        "intensity_": 5175.5
    },
    {
        "displayLevel_": 3,
        "peakId_": "390",
        "pos_": 589.6543,
        "monoMz_": 589.6543,
        "intensity_": 4969.3
    },
    {
        "displayLevel_": 6,
        "peakId_": "70",
        "pos_": 413.9083,
        "monoMz_": 413.9083,
        "intensity_": 4929.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "597",
        "pos_": 627.9612,
        "monoMz_": 627.9612,
        "intensity_": 4907
    },
    {
        "displayLevel_": 1,
        "peakId_": "579",
        "pos_": 625.3343,
        "monoMz_": 625.3343,
        "intensity_": 4901.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "346",
        "pos_": 578.333,
        "monoMz_": 578.333,
        "intensity_": 4898.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "71",
        "pos_": 414.243,
        "monoMz_": 414.243,
        "intensity_": 4888.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "63",
        "pos_": 366.2677,
        "monoMz_": 366.2677,
        "intensity_": 4840.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "815",
        "pos_": 713.6864,
        "monoMz_": 713.6864,
        "intensity_": 4833.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "463",
        "pos_": 604.3471,
        "monoMz_": 604.3471,
        "intensity_": 4816.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "128",
        "pos_": 493.2825,
        "monoMz_": 493.2825,
        "intensity_": 4812.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "485",
        "pos_": 606.3347,
        "monoMz_": 606.3347,
        "intensity_": 4788
    },
    {
        "displayLevel_": 0,
        "peakId_": "334",
        "pos_": 577.1401,
        "monoMz_": 577.1401,
        "intensity_": 4786.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "248",
        "pos_": 561.1306,
        "monoMz_": 561.1306,
        "intensity_": 4770
    },
    {
        "displayLevel_": 2,
        "peakId_": "608",
        "pos_": 631.3362,
        "monoMz_": 631.3362,
        "intensity_": 4764.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "480",
        "pos_": 605.933,
        "monoMz_": 605.933,
        "intensity_": 4755.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "573",
        "pos_": 624.7283,
        "monoMz_": 624.7283,
        "intensity_": 4689.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "674",
        "pos_": 649.0991,
        "monoMz_": 649.0991,
        "intensity_": 4553.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "309",
        "pos_": 574.642,
        "monoMz_": 574.642,
        "intensity_": 4550.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "388",
        "pos_": 589.4302,
        "monoMz_": 589.4302,
        "intensity_": 4377.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "653",
        "pos_": 643.3462,
        "monoMz_": 643.3462,
        "intensity_": 4375
    },
    {
        "displayLevel_": 5,
        "peakId_": "83",
        "pos_": 461.9211,
        "monoMz_": 461.9211,
        "intensity_": 4356.4
    },
    {
        "displayLevel_": 1,
        "peakId_": "231",
        "pos_": 558.0532,
        "monoMz_": 558.0532,
        "intensity_": 4332.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "490",
        "pos_": 606.707,
        "monoMz_": 606.707,
        "intensity_": 4319
    },
    {
        "displayLevel_": 1,
        "peakId_": "297",
        "pos_": 572.9743,
        "monoMz_": 572.9743,
        "intensity_": 4221.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "348",
        "pos_": 578.5328,
        "monoMz_": 578.5328,
        "intensity_": 4216
    },
    {
        "displayLevel_": 0,
        "peakId_": "275",
        "pos_": 568.6635,
        "monoMz_": 568.6635,
        "intensity_": 4146.9
    },
    {
        "displayLevel_": 3,
        "peakId_": "690",
        "pos_": 655.1476,
        "monoMz_": 655.1476,
        "intensity_": 4127.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "619",
        "pos_": 634.0042,
        "monoMz_": 634.0042,
        "intensity_": 4125.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "915",
        "pos_": 1877.2362,
        "monoMz_": 1877.2362,
        "intensity_": 4123
    },
    {
        "displayLevel_": 2,
        "peakId_": "636",
        "pos_": 637.0378,
        "monoMz_": 637.0378,
        "intensity_": 4122.1
    },
    {
        "displayLevel_": 4,
        "peakId_": "182",
        "pos_": 543.6466,
        "monoMz_": 543.6466,
        "intensity_": 4106.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "369",
        "pos_": 581.9625,
        "monoMz_": 581.9625,
        "intensity_": 4049.4
    },
    {
        "displayLevel_": 4,
        "peakId_": "753",
        "pos_": 677.0254,
        "monoMz_": 677.0254,
        "intensity_": 4043.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "362",
        "pos_": 581.3359,
        "monoMz_": 581.3359,
        "intensity_": 4025.7
    },
    {
        "displayLevel_": 3,
        "peakId_": "434",
        "pos_": 601.4847,
        "monoMz_": 601.4847,
        "intensity_": 4025.4
    },
    {
        "displayLevel_": 3,
        "peakId_": "770",
        "pos_": 686.6495,
        "monoMz_": 686.6495,
        "intensity_": 4018.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "703",
        "pos_": 658.0607,
        "monoMz_": 658.0607,
        "intensity_": 4013.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "607",
        "pos_": 631.1688,
        "monoMz_": 631.1688,
        "intensity_": 4010.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "387",
        "pos_": 589.3201,
        "monoMz_": 589.3201,
        "intensity_": 4009.5
    },
    {
        "displayLevel_": 7,
        "peakId_": "892",
        "pos_": 1198.0404,
        "monoMz_": 1198.0404,
        "intensity_": 3989.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "523",
        "pos_": 611.2252,
        "monoMz_": 611.2252,
        "intensity_": 3943.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "689",
        "pos_": 654.947,
        "monoMz_": 654.947,
        "intensity_": 3867.2
    },
    {
        "displayLevel_": 3,
        "peakId_": "720",
        "pos_": 667.3771,
        "monoMz_": 667.3771,
        "intensity_": 3836.4
    },
    {
        "displayLevel_": 2,
        "peakId_": "711",
        "pos_": 664.5269,
        "monoMz_": 664.5269,
        "intensity_": 3835.4
    },
    {
        "displayLevel_": 1,
        "peakId_": "751",
        "pos_": 676.6915,
        "monoMz_": 676.6915,
        "intensity_": 3793.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "702",
        "pos_": 657.8945,
        "monoMz_": 657.8945,
        "intensity_": 3771.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "296",
        "pos_": 572.8635,
        "monoMz_": 572.8635,
        "intensity_": 3771
    },
    {
        "displayLevel_": 2,
        "peakId_": "838",
        "pos_": 722.9141,
        "monoMz_": 722.9141,
        "intensity_": 3763.4
    },
    {
        "displayLevel_": 2,
        "peakId_": "761",
        "pos_": 683.3734,
        "monoMz_": 683.3734,
        "intensity_": 3676.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "769",
        "pos_": 686.5067,
        "monoMz_": 686.5067,
        "intensity_": 3676.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "637",
        "pos_": 637.1806,
        "monoMz_": 637.1806,
        "intensity_": 3599.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "609",
        "pos_": 631.5041,
        "monoMz_": 631.5041,
        "intensity_": 3590.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "230",
        "pos_": 557.965,
        "monoMz_": 557.965,
        "intensity_": 3568.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "324",
        "pos_": 576.5242,
        "monoMz_": 576.5242,
        "intensity_": 3542.6
    },
    {
        "displayLevel_": 2,
        "peakId_": "726",
        "pos_": 668.9531,
        "monoMz_": 668.9531,
        "intensity_": 3499.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "580",
        "pos_": 625.4587,
        "monoMz_": 625.4587,
        "intensity_": 3495.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "154",
        "pos_": 505.6133,
        "monoMz_": 505.6133,
        "intensity_": 3473.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "697",
        "pos_": 657.3738,
        "monoMz_": 657.3738,
        "intensity_": 3466.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "435",
        "pos_": 601.6529,
        "monoMz_": 601.6529,
        "intensity_": 3403.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "276",
        "pos_": 568.8302,
        "monoMz_": 568.8302,
        "intensity_": 3388.9
    },
    {
        "displayLevel_": 1,
        "peakId_": "742",
        "pos_": 672.9543,
        "monoMz_": 672.9543,
        "intensity_": 3380.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "816",
        "pos_": 713.8322,
        "monoMz_": 713.8322,
        "intensity_": 3380.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "328",
        "pos_": 576.8076,
        "monoMz_": 576.8076,
        "intensity_": 3376.4
    },
    {
        "displayLevel_": 4,
        "peakId_": "107",
        "pos_": 483.2721,
        "monoMz_": 483.2721,
        "intensity_": 3332.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "512",
        "pos_": 608.5283,
        "monoMz_": 608.5283,
        "intensity_": 3307.8
    },
    {
        "displayLevel_": 3,
        "peakId_": "401",
        "pos_": 593.8442,
        "monoMz_": 593.8442,
        "intensity_": 3300.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "389",
        "pos_": 589.5416,
        "monoMz_": 589.5416,
        "intensity_": 3280
    },
    {
        "displayLevel_": 0,
        "peakId_": "220",
        "pos_": 556.7116,
        "monoMz_": 556.7116,
        "intensity_": 3279.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "108",
        "pos_": 483.4395,
        "monoMz_": 483.4395,
        "intensity_": 3275.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "574",
        "pos_": 624.8536,
        "monoMz_": 624.8536,
        "intensity_": 3273.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "829",
        "pos_": 720.3749,
        "monoMz_": 720.3749,
        "intensity_": 3238.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "180",
        "pos_": 543.5351,
        "monoMz_": 543.5351,
        "intensity_": 3230.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "592",
        "pos_": 627.3342,
        "monoMz_": 627.3342,
        "intensity_": 3226.4
    },
    {
        "displayLevel_": 6,
        "peakId_": "24",
        "pos_": 199.2472,
        "monoMz_": 199.2472,
        "intensity_": 3220.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "522",
        "pos_": 611.1008,
        "monoMz_": 611.1008,
        "intensity_": 3198.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "785",
        "pos_": 693.3811,
        "monoMz_": 693.3811,
        "intensity_": 3186
    },
    {
        "displayLevel_": 0,
        "peakId_": "497",
        "pos_": 607.3271,
        "monoMz_": 607.3271,
        "intensity_": 3184.6
    },
    {
        "displayLevel_": 1,
        "peakId_": "837",
        "pos_": 722.6997,
        "monoMz_": 722.6997,
        "intensity_": 3182.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "219",
        "pos_": 556.5094,
        "monoMz_": 556.5094,
        "intensity_": 3145.3
    },
    {
        "displayLevel_": 2,
        "peakId_": "684",
        "pos_": 652.3293,
        "monoMz_": 652.3293,
        "intensity_": 3101.6
    },
    {
        "displayLevel_": 2,
        "peakId_": "776",
        "pos_": 689.6772,
        "monoMz_": 689.6772,
        "intensity_": 3094.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "318",
        "pos_": 575.4209,
        "monoMz_": 575.4209,
        "intensity_": 3091.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "727",
        "pos_": 668.9756,
        "monoMz_": 668.9756,
        "intensity_": 3070.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "627",
        "pos_": 634.8401,
        "monoMz_": 634.8401,
        "intensity_": 3045.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "298",
        "pos_": 573.0863,
        "monoMz_": 573.0863,
        "intensity_": 3034.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "724",
        "pos_": 668.752,
        "monoMz_": 668.752,
        "intensity_": 2963.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "704",
        "pos_": 658.1457,
        "monoMz_": 658.1457,
        "intensity_": 2934
    },
    {
        "displayLevel_": 0,
        "peakId_": "598",
        "pos_": 628.0872,
        "monoMz_": 628.0872,
        "intensity_": 2931.7
    },
    {
        "displayLevel_": 1,
        "peakId_": "616",
        "pos_": 632.3466,
        "monoMz_": 632.3466,
        "intensity_": 2920.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "235",
        "pos_": 558.3041,
        "monoMz_": 558.3041,
        "intensity_": 2914.6
    },
    {
        "displayLevel_": 1,
        "peakId_": "763",
        "pos_": 683.8725,
        "monoMz_": 683.8725,
        "intensity_": 2896.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "447",
        "pos_": 602.7941,
        "monoMz_": 602.7941,
        "intensity_": 2878.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "514",
        "pos_": 608.7234,
        "monoMz_": 608.7234,
        "intensity_": 2873.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "360",
        "pos_": 580.5284,
        "monoMz_": 580.5284,
        "intensity_": 2859.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "818",
        "pos_": 714.0356,
        "monoMz_": 714.0356,
        "intensity_": 2840.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "213",
        "pos_": 555.6489,
        "monoMz_": 555.6489,
        "intensity_": 2831.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "234",
        "pos_": 558.1893,
        "monoMz_": 558.1893,
        "intensity_": 2815.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "787",
        "pos_": 694.3964,
        "monoMz_": 694.3964,
        "intensity_": 2802.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "718",
        "pos_": 667.177,
        "monoMz_": 667.177,
        "intensity_": 2798.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "736",
        "pos_": 671.0231,
        "monoMz_": 671.0231,
        "intensity_": 2764.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "638",
        "pos_": 637.3224,
        "monoMz_": 637.3224,
        "intensity_": 2753.4
    },
    {
        "displayLevel_": 3,
        "peakId_": "142",
        "pos_": 499.7845,
        "monoMz_": 499.7845,
        "intensity_": 2749.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "728",
        "pos_": 669.0801,
        "monoMz_": 669.0801,
        "intensity_": 2742
    },
    {
        "displayLevel_": 0,
        "peakId_": "224",
        "pos_": 557.2936,
        "monoMz_": 557.2936,
        "intensity_": 2738.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "699",
        "pos_": 657.6433,
        "monoMz_": 657.6433,
        "intensity_": 2720.4
    },
    {
        "displayLevel_": 1,
        "peakId_": "666",
        "pos_": 647.5658,
        "monoMz_": 647.5658,
        "intensity_": 2711.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "525",
        "pos_": 611.4763,
        "monoMz_": 611.4763,
        "intensity_": 2687
    },
    {
        "displayLevel_": 2,
        "peakId_": "754",
        "pos_": 677.359,
        "monoMz_": 677.359,
        "intensity_": 2682.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "109",
        "pos_": 483.6065,
        "monoMz_": 483.6065,
        "intensity_": 2671.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "402",
        "pos_": 593.9696,
        "monoMz_": 593.9696,
        "intensity_": 2661.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "613",
        "pos_": 631.9417,
        "monoMz_": 631.9417,
        "intensity_": 2656.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "399",
        "pos_": 593.5932,
        "monoMz_": 593.5932,
        "intensity_": 2649.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "601",
        "pos_": 628.3378,
        "monoMz_": 628.3378,
        "intensity_": 2648.3
    },
    {
        "displayLevel_": 1,
        "peakId_": "392",
        "pos_": 589.876,
        "monoMz_": 589.876,
        "intensity_": 2640.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "400",
        "pos_": 593.7189,
        "monoMz_": 593.7189,
        "intensity_": 2640.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "295",
        "pos_": 572.7524,
        "monoMz_": 572.7524,
        "intensity_": 2623.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "729",
        "pos_": 669.1556,
        "monoMz_": 669.1556,
        "intensity_": 2621.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "705",
        "pos_": 658.312,
        "monoMz_": 658.312,
        "intensity_": 2619.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "84",
        "pos_": 462.2557,
        "monoMz_": 462.2557,
        "intensity_": 2592.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "150",
        "pos_": 505.1121,
        "monoMz_": 505.1121,
        "intensity_": 2589.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "344",
        "pos_": 578.1328,
        "monoMz_": 578.1328,
        "intensity_": 2588.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "183",
        "pos_": 543.7203,
        "monoMz_": 543.7203,
        "intensity_": 2563.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "152",
        "pos_": 505.4456,
        "monoMz_": 505.4456,
        "intensity_": 2555.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "618",
        "pos_": 633.8403,
        "monoMz_": 633.8403,
        "intensity_": 2554.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "615",
        "pos_": 632.143,
        "monoMz_": 632.143,
        "intensity_": 2545
    },
    {
        "displayLevel_": 0,
        "peakId_": "665",
        "pos_": 647.4651,
        "monoMz_": 647.4651,
        "intensity_": 2531.5
    },
    {
        "displayLevel_": 5,
        "peakId_": "859",
        "pos_": 768.3751,
        "monoMz_": 768.3751,
        "intensity_": 2524.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "611",
        "pos_": 631.7428,
        "monoMz_": 631.7428,
        "intensity_": 2512.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "836",
        "pos_": 722.6633,
        "monoMz_": 722.6633,
        "intensity_": 2511.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "678",
        "pos_": 649.3669,
        "monoMz_": 649.3669,
        "intensity_": 2500.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "227",
        "pos_": 557.791,
        "monoMz_": 557.791,
        "intensity_": 2473.1
    },
    {
        "displayLevel_": 1,
        "peakId_": "187",
        "pos_": 544.1491,
        "monoMz_": 544.1491,
        "intensity_": 2467.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "526",
        "pos_": 611.6015,
        "monoMz_": 611.6015,
        "intensity_": 2453.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "214",
        "pos_": 555.7609,
        "monoMz_": 555.7609,
        "intensity_": 2452.8
    },
    {
        "displayLevel_": 2,
        "peakId_": "90",
        "pos_": 466.9282,
        "monoMz_": 466.9282,
        "intensity_": 2452.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "716",
        "pos_": 666.977,
        "monoMz_": 666.977,
        "intensity_": 2449.4
    },
    {
        "displayLevel_": 2,
        "peakId_": "373",
        "pos_": 584.6701,
        "monoMz_": 584.6701,
        "intensity_": 2439.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "817",
        "pos_": 713.9767,
        "monoMz_": 713.9767,
        "intensity_": 2436
    },
    {
        "displayLevel_": 0,
        "peakId_": "391",
        "pos_": 589.7635,
        "monoMz_": 589.7635,
        "intensity_": 2429.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "449",
        "pos_": 603.045,
        "monoMz_": 603.045,
        "intensity_": 2421.2
    },
    {
        "displayLevel_": 1,
        "peakId_": "516",
        "pos_": 608.85,
        "monoMz_": 608.85,
        "intensity_": 2418.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "442",
        "pos_": 602.4314,
        "monoMz_": 602.4314,
        "intensity_": 2386.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "244",
        "pos_": 560.7131,
        "monoMz_": 560.7131,
        "intensity_": 2366.5
    },
    {
        "displayLevel_": 3,
        "peakId_": "801",
        "pos_": 705.3867,
        "monoMz_": 705.3867,
        "intensity_": 2365.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "129",
        "pos_": 493.5333,
        "monoMz_": 493.5333,
        "intensity_": 2363.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "492",
        "pos_": 606.7988,
        "monoMz_": 606.7988,
        "intensity_": 2360.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "713",
        "pos_": 664.8138,
        "monoMz_": 664.8138,
        "intensity_": 2348
    },
    {
        "displayLevel_": 0,
        "peakId_": "835",
        "pos_": 722.4124,
        "monoMz_": 722.4124,
        "intensity_": 2338.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "719",
        "pos_": 667.2778,
        "monoMz_": 667.2778,
        "intensity_": 2324.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "141",
        "pos_": 499.6839,
        "monoMz_": 499.6839,
        "intensity_": 2276.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "730",
        "pos_": 669.28,
        "monoMz_": 669.28,
        "intensity_": 2269.2
    },
    {
        "displayLevel_": 5,
        "peakId_": "868",
        "pos_": 790.4379,
        "monoMz_": 790.4379,
        "intensity_": 2261
    },
    {
        "displayLevel_": 0,
        "peakId_": "626",
        "pos_": 634.7495,
        "monoMz_": 634.7495,
        "intensity_": 2258.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "891",
        "pos_": 1197.9624,
        "monoMz_": 1197.9624,
        "intensity_": 2251.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "602",
        "pos_": 628.5405,
        "monoMz_": 628.5405,
        "intensity_": 2243.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "624",
        "pos_": 634.6066,
        "monoMz_": 634.6066,
        "intensity_": 2243.7
    },
    {
        "displayLevel_": 1,
        "peakId_": "857",
        "pos_": 766.3801,
        "monoMz_": 766.3801,
        "intensity_": 2219.2
    },
    {
        "displayLevel_": 2,
        "peakId_": "798",
        "pos_": 698.3994,
        "monoMz_": 698.3994,
        "intensity_": 2215.4
    },
    {
        "displayLevel_": 3,
        "peakId_": "542",
        "pos_": 615.0798,
        "monoMz_": 615.0798,
        "intensity_": 2208.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "494",
        "pos_": 606.9814,
        "monoMz_": 606.9814,
        "intensity_": 2198.9
    },
    {
        "displayLevel_": 5,
        "peakId_": "80",
        "pos_": 456.2509,
        "monoMz_": 456.2509,
        "intensity_": 2197.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "610",
        "pos_": 631.6696,
        "monoMz_": 631.6696,
        "intensity_": 2193.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "228",
        "pos_": 557.8051,
        "monoMz_": 557.8051,
        "intensity_": 2186.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "488",
        "pos_": 606.5381,
        "monoMz_": 606.5381,
        "intensity_": 2167.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "339",
        "pos_": 577.4743,
        "monoMz_": 577.4743,
        "intensity_": 2164.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "322",
        "pos_": 576.3223,
        "monoMz_": 576.3223,
        "intensity_": 2161.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "794",
        "pos_": 695.2565,
        "monoMz_": 695.2565,
        "intensity_": 2161.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "775",
        "pos_": 689.3441,
        "monoMz_": 689.3441,
        "intensity_": 2126.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "185",
        "pos_": 543.863,
        "monoMz_": 543.863,
        "intensity_": 2112.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "581",
        "pos_": 625.5159,
        "monoMz_": 625.5159,
        "intensity_": 2109
    },
    {
        "displayLevel_": 0,
        "peakId_": "22",
        "pos_": 199.2377,
        "monoMz_": 199.2377,
        "intensity_": 2099.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "430",
        "pos_": 600.9452,
        "monoMz_": 600.9452,
        "intensity_": 2096.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "706",
        "pos_": 658.3939,
        "monoMz_": 658.3939,
        "intensity_": 2085.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "606",
        "pos_": 631.003,
        "monoMz_": 631.003,
        "intensity_": 2077.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "110",
        "pos_": 483.7744,
        "monoMz_": 483.7744,
        "intensity_": 2063.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "340",
        "pos_": 577.5344,
        "monoMz_": 577.5344,
        "intensity_": 2048.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "823",
        "pos_": 714.688,
        "monoMz_": 714.688,
        "intensity_": 2047.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "277",
        "pos_": 568.9987,
        "monoMz_": 568.9987,
        "intensity_": 2045.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "831",
        "pos_": 721.5798,
        "monoMz_": 721.5798,
        "intensity_": 2044.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "436",
        "pos_": 601.8209,
        "monoMz_": 601.8209,
        "intensity_": 2041.4
    },
    {
        "displayLevel_": 1,
        "peakId_": "687",
        "pos_": 654.746,
        "monoMz_": 654.746,
        "intensity_": 2038.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "186",
        "pos_": 544.0067,
        "monoMz_": 544.0067,
        "intensity_": 2031
    },
    {
        "displayLevel_": 0,
        "peakId_": "481",
        "pos_": 605.967,
        "monoMz_": 605.967,
        "intensity_": 2030.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "834",
        "pos_": 722.365,
        "monoMz_": 722.365,
        "intensity_": 2019.6
    },
    {
        "displayLevel_": 3,
        "peakId_": "192",
        "pos_": 550.3151,
        "monoMz_": 550.3151,
        "intensity_": 2019
    },
    {
        "displayLevel_": 0,
        "peakId_": "840",
        "pos_": 723.1641,
        "monoMz_": 723.1641,
        "intensity_": 2013.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "861",
        "pos_": 768.6241,
        "monoMz_": 768.6241,
        "intensity_": 2011.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "305",
        "pos_": 574.1029,
        "monoMz_": 574.1029,
        "intensity_": 2003
    },
    {
        "displayLevel_": 0,
        "peakId_": "710",
        "pos_": 664.3856,
        "monoMz_": 664.3856,
        "intensity_": 2001
    },
    {
        "displayLevel_": 2,
        "peakId_": "289",
        "pos_": 571.6438,
        "monoMz_": 571.6438,
        "intensity_": 1992.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "735",
        "pos_": 670.6874,
        "monoMz_": 670.6874,
        "intensity_": 1979.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "668",
        "pos_": 647.7647,
        "monoMz_": 647.7647,
        "intensity_": 1976.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "462",
        "pos_": 604.25,
        "monoMz_": 604.25,
        "intensity_": 1970.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "677",
        "pos_": 649.3506,
        "monoMz_": 649.3506,
        "intensity_": 1969.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "639",
        "pos_": 637.467,
        "monoMz_": 637.467,
        "intensity_": 1961.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "759",
        "pos_": 682.5781,
        "monoMz_": 682.5781,
        "intensity_": 1950.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "503",
        "pos_": 608.0096,
        "monoMz_": 608.0096,
        "intensity_": 1948.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "712",
        "pos_": 664.6698,
        "monoMz_": 664.6698,
        "intensity_": 1920.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "586",
        "pos_": 625.8499,
        "monoMz_": 625.8499,
        "intensity_": 1920
    },
    {
        "displayLevel_": 0,
        "peakId_": "284",
        "pos_": 570.6484,
        "monoMz_": 570.6484,
        "intensity_": 1917.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "57",
        "pos_": 355.7345,
        "monoMz_": 355.7345,
        "intensity_": 1916.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "474",
        "pos_": 605.6178,
        "monoMz_": 605.6178,
        "intensity_": 1894.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "869",
        "pos_": 790.938,
        "monoMz_": 790.938,
        "intensity_": 1887.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "722",
        "pos_": 668.5513,
        "monoMz_": 668.5513,
        "intensity_": 1878.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "477",
        "pos_": 605.8005,
        "monoMz_": 605.8005,
        "intensity_": 1874.4
    },
    {
        "displayLevel_": 6,
        "peakId_": "903",
        "pos_": 1717.957,
        "monoMz_": 1717.957,
        "intensity_": 1873.9
    },
    {
        "displayLevel_": 3,
        "peakId_": "842",
        "pos_": 724.4103,
        "monoMz_": 724.4103,
        "intensity_": 1871.4
    },
    {
        "displayLevel_": 1,
        "peakId_": "88",
        "pos_": 464.7726,
        "monoMz_": 464.7726,
        "intensity_": 1867.4
    },
    {
        "displayLevel_": 2,
        "peakId_": "662",
        "pos_": 646.7238,
        "monoMz_": 646.7238,
        "intensity_": 1857.8
    },
    {
        "displayLevel_": 2,
        "peakId_": "795",
        "pos_": 695.3992,
        "monoMz_": 695.3992,
        "intensity_": 1844.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "802",
        "pos_": 705.8871,
        "monoMz_": 705.8871,
        "intensity_": 1828.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "599",
        "pos_": 628.1415,
        "monoMz_": 628.1415,
        "intensity_": 1822.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "578",
        "pos_": 625.2078,
        "monoMz_": 625.2078,
        "intensity_": 1808.9
    },
    {
        "displayLevel_": 6,
        "peakId_": "876",
        "pos_": 863.9684,
        "monoMz_": 863.9684,
        "intensity_": 1804.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "17",
        "pos_": 198.8254,
        "monoMz_": 198.8254,
        "intensity_": 1800.4
    },
    {
        "displayLevel_": 3,
        "peakId_": "114",
        "pos_": 486.6802,
        "monoMz_": 486.6802,
        "intensity_": 1784.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "679",
        "pos_": 649.4708,
        "monoMz_": 649.4708,
        "intensity_": 1782
    },
    {
        "displayLevel_": 0,
        "peakId_": "764",
        "pos_": 683.9307,
        "monoMz_": 683.9307,
        "intensity_": 1774.8
    },
    {
        "displayLevel_": 4,
        "peakId_": "99",
        "pos_": 480.3974,
        "monoMz_": 480.3974,
        "intensity_": 1773.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "178",
        "pos_": 543.3124,
        "monoMz_": 543.3124,
        "intensity_": 1765.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "877",
        "pos_": 864.4744,
        "monoMz_": 864.4744,
        "intensity_": 1764
    },
    {
        "displayLevel_": 1,
        "peakId_": "655",
        "pos_": 643.8473,
        "monoMz_": 643.8473,
        "intensity_": 1760.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "584",
        "pos_": 625.7087,
        "monoMz_": 625.7087,
        "intensity_": 1745.9
    },
    {
        "displayLevel_": 5,
        "peakId_": "884",
        "pos_": 877.08,
        "monoMz_": 877.08,
        "intensity_": 1742.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "370",
        "pos_": 582.089,
        "monoMz_": 582.089,
        "intensity_": 1741
    },
    {
        "displayLevel_": 0,
        "peakId_": "714",
        "pos_": 665.3462,
        "monoMz_": 665.3462,
        "intensity_": 1735.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "700",
        "pos_": 657.7097,
        "monoMz_": 657.7097,
        "intensity_": 1723.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "143",
        "pos_": 499.8849,
        "monoMz_": 499.8849,
        "intensity_": 1716.2
    },
    {
        "displayLevel_": 1,
        "peakId_": "777",
        "pos_": 690.013,
        "monoMz_": 690.013,
        "intensity_": 1707.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "273",
        "pos_": 568.3272,
        "monoMz_": 568.3272,
        "intensity_": 1705.8
    },
    {
        "displayLevel_": 2,
        "peakId_": "862",
        "pos_": 768.8729,
        "monoMz_": 768.8729,
        "intensity_": 1702.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "144",
        "pos_": 499.9854,
        "monoMz_": 499.9854,
        "intensity_": 1699.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "429",
        "pos_": 600.8199,
        "monoMz_": 600.8199,
        "intensity_": 1698.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "766",
        "pos_": 684.0778,
        "monoMz_": 684.0778,
        "intensity_": 1680.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "279",
        "pos_": 570.0052,
        "monoMz_": 570.0052,
        "intensity_": 1677.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "371",
        "pos_": 582.3683,
        "monoMz_": 582.3683,
        "intensity_": 1676.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "641",
        "pos_": 637.6118,
        "monoMz_": 637.6118,
        "intensity_": 1673
    },
    {
        "displayLevel_": 1,
        "peakId_": "266",
        "pos_": 566.4363,
        "monoMz_": 566.4363,
        "intensity_": 1663.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "696",
        "pos_": 657.0393,
        "monoMz_": 657.0393,
        "intensity_": 1654.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "443",
        "pos_": 602.5295,
        "monoMz_": 602.5295,
        "intensity_": 1654
    },
    {
        "displayLevel_": 0,
        "peakId_": "113",
        "pos_": 486.4792,
        "monoMz_": 486.4792,
        "intensity_": 1653.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "768",
        "pos_": 686.3618,
        "monoMz_": 686.3618,
        "intensity_": 1652.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "830",
        "pos_": 720.7087,
        "monoMz_": 720.7087,
        "intensity_": 1631
    },
    {
        "displayLevel_": 4,
        "peakId_": "845",
        "pos_": 741.3104,
        "monoMz_": 741.3104,
        "intensity_": 1627
    },
    {
        "displayLevel_": 0,
        "peakId_": "303",
        "pos_": 573.8801,
        "monoMz_": 573.8801,
        "intensity_": 1618.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "467",
        "pos_": 604.6007,
        "monoMz_": 604.6007,
        "intensity_": 1601.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "550",
        "pos_": 620.9792,
        "monoMz_": 620.9792,
        "intensity_": 1599.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "105",
        "pos_": 482.2784,
        "monoMz_": 482.2784,
        "intensity_": 1598.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "786",
        "pos_": 693.8843,
        "monoMz_": 693.8843,
        "intensity_": 1587
    },
    {
        "displayLevel_": 0,
        "peakId_": "433",
        "pos_": 601.3179,
        "monoMz_": 601.3179,
        "intensity_": 1584.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "832",
        "pos_": 721.7806,
        "monoMz_": 721.7806,
        "intensity_": 1579.1
    },
    {
        "displayLevel_": 1,
        "peakId_": "530",
        "pos_": 613.3531,
        "monoMz_": 613.3531,
        "intensity_": 1578.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "155",
        "pos_": 505.7805,
        "monoMz_": 505.7805,
        "intensity_": 1577.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "321",
        "pos_": 576.1237,
        "monoMz_": 576.1237,
        "intensity_": 1575.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "159",
        "pos_": 506.5917,
        "monoMz_": 506.5917,
        "intensity_": 1575.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "635",
        "pos_": 636.8926,
        "monoMz_": 636.8926,
        "intensity_": 1568.2
    },
    {
        "displayLevel_": 2,
        "peakId_": "415",
        "pos_": 598.3206,
        "monoMz_": 598.3206,
        "intensity_": 1567.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "833",
        "pos_": 721.9831,
        "monoMz_": 721.9831,
        "intensity_": 1566.1
    },
    {
        "displayLevel_": 3,
        "peakId_": "169",
        "pos_": 532.4794,
        "monoMz_": 532.4794,
        "intensity_": 1565.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "743",
        "pos_": 673.1565,
        "monoMz_": 673.1565,
        "intensity_": 1562.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "622",
        "pos_": 634.4666,
        "monoMz_": 634.4666,
        "intensity_": 1558.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "780",
        "pos_": 692.0521,
        "monoMz_": 692.0521,
        "intensity_": 1557.9
    },
    {
        "displayLevel_": 1,
        "peakId_": "59",
        "pos_": 357.2623,
        "monoMz_": 357.2623,
        "intensity_": 1556.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "535",
        "pos_": 614.3445,
        "monoMz_": 614.3445,
        "intensity_": 1550.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "79",
        "pos_": 455.9184,
        "monoMz_": 455.9184,
        "intensity_": 1550.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "537",
        "pos_": 614.6791,
        "monoMz_": 614.6791,
        "intensity_": 1544.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "372",
        "pos_": 584.4199,
        "monoMz_": 584.4199,
        "intensity_": 1541.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "326",
        "pos_": 576.7267,
        "monoMz_": 576.7267,
        "intensity_": 1539.2
    },
    {
        "displayLevel_": 2,
        "peakId_": "163",
        "pos_": 508.1919,
        "monoMz_": 508.1919,
        "intensity_": 1536.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "858",
        "pos_": 768.1231,
        "monoMz_": 768.1231,
        "intensity_": 1528.9
    },
    {
        "displayLevel_": 2,
        "peakId_": "148",
        "pos_": 502.6122,
        "monoMz_": 502.6122,
        "intensity_": 1524.5
    },
    {
        "displayLevel_": 1,
        "peakId_": "808",
        "pos_": 710.3977,
        "monoMz_": 710.3977,
        "intensity_": 1517.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "493",
        "pos_": 606.8913,
        "monoMz_": 606.8913,
        "intensity_": 1515.6
    },
    {
        "displayLevel_": 2,
        "peakId_": "122",
        "pos_": 489.0269,
        "monoMz_": 489.0269,
        "intensity_": 1515.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "465",
        "pos_": 604.4746,
        "monoMz_": 604.4746,
        "intensity_": 1514.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "664",
        "pos_": 647.3642,
        "monoMz_": 647.3642,
        "intensity_": 1496.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "737",
        "pos_": 671.3553,
        "monoMz_": 671.3553,
        "intensity_": 1494.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "291",
        "pos_": 571.9797,
        "monoMz_": 571.9797,
        "intensity_": 1491.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "629",
        "pos_": 635.0068,
        "monoMz_": 635.0068,
        "intensity_": 1483.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "748",
        "pos_": 675.7416,
        "monoMz_": 675.7416,
        "intensity_": 1477.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "350",
        "pos_": 578.7354,
        "monoMz_": 578.7354,
        "intensity_": 1476.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "515",
        "pos_": 608.7997,
        "monoMz_": 608.7997,
        "intensity_": 1475.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "281",
        "pos_": 570.3145,
        "monoMz_": 570.3145,
        "intensity_": 1470.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "781",
        "pos_": 692.1637,
        "monoMz_": 692.1637,
        "intensity_": 1469.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "448",
        "pos_": 602.9349,
        "monoMz_": 602.9349,
        "intensity_": 1463.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "249",
        "pos_": 561.3058,
        "monoMz_": 561.3058,
        "intensity_": 1461.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "237",
        "pos_": 558.5533,
        "monoMz_": 558.5533,
        "intensity_": 1448.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "640",
        "pos_": 637.5573,
        "monoMz_": 637.5573,
        "intensity_": 1448.5
    },
    {
        "displayLevel_": 2,
        "peakId_": "133",
        "pos_": 497.2877,
        "monoMz_": 497.2877,
        "intensity_": 1434.4
    },
    {
        "displayLevel_": 4,
        "peakId_": "865",
        "pos_": 778.9199,
        "monoMz_": 778.9199,
        "intensity_": 1428.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "771",
        "pos_": 686.7937,
        "monoMz_": 686.7937,
        "intensity_": 1427.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "583",
        "pos_": 625.6892,
        "monoMz_": 625.6892,
        "intensity_": 1417.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "745",
        "pos_": 674.183,
        "monoMz_": 674.183,
        "intensity_": 1415.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "663",
        "pos_": 646.8443,
        "monoMz_": 646.8443,
        "intensity_": 1408.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "670",
        "pos_": 648.5986,
        "monoMz_": 648.5986,
        "intensity_": 1402.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "707",
        "pos_": 658.5623,
        "monoMz_": 658.5623,
        "intensity_": 1401.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "222",
        "pos_": 556.9122,
        "monoMz_": 556.9122,
        "intensity_": 1393.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "101",
        "pos_": 480.5402,
        "monoMz_": 480.5402,
        "intensity_": 1392
    },
    {
        "displayLevel_": 0,
        "peakId_": "398",
        "pos_": 593.4665,
        "monoMz_": 593.4665,
        "intensity_": 1392
    },
    {
        "displayLevel_": 0,
        "peakId_": "432",
        "pos_": 601.0698,
        "monoMz_": 601.0698,
        "intensity_": 1389
    },
    {
        "displayLevel_": 0,
        "peakId_": "476",
        "pos_": 605.7258,
        "monoMz_": 605.7258,
        "intensity_": 1386
    },
    {
        "displayLevel_": 0,
        "peakId_": "772",
        "pos_": 686.9371,
        "monoMz_": 686.9371,
        "intensity_": 1380.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "459",
        "pos_": 604.0706,
        "monoMz_": 604.0706,
        "intensity_": 1378.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "304",
        "pos_": 573.9919,
        "monoMz_": 573.9919,
        "intensity_": 1374.4
    },
    {
        "displayLevel_": 2,
        "peakId_": "173",
        "pos_": 540.6437,
        "monoMz_": 540.6437,
        "intensity_": 1371.6
    },
    {
        "displayLevel_": 1,
        "peakId_": "428",
        "pos_": 600.6854,
        "monoMz_": 600.6854,
        "intensity_": 1368.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "420",
        "pos_": 598.8234,
        "monoMz_": 598.8234,
        "intensity_": 1368.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "100",
        "pos_": 480.4391,
        "monoMz_": 480.4391,
        "intensity_": 1365.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "558",
        "pos_": 622.2248,
        "monoMz_": 622.2248,
        "intensity_": 1362.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "469",
        "pos_": 604.9775,
        "monoMz_": 604.9775,
        "intensity_": 1360.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "612",
        "pos_": 631.8373,
        "monoMz_": 631.8373,
        "intensity_": 1359.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "168",
        "pos_": 532.4577,
        "monoMz_": 532.4577,
        "intensity_": 1358.2
    },
    {
        "displayLevel_": 3,
        "peakId_": "94",
        "pos_": 470.5913,
        "monoMz_": 470.5913,
        "intensity_": 1357.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "902",
        "pos_": 1717.8245,
        "monoMz_": 1717.8245,
        "intensity_": 1356
    },
    {
        "displayLevel_": 0,
        "peakId_": "131",
        "pos_": 493.7831,
        "monoMz_": 493.7831,
        "intensity_": 1355.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "688",
        "pos_": 654.861,
        "monoMz_": 654.861,
        "intensity_": 1353.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "733",
        "pos_": 669.6807,
        "monoMz_": 669.6807,
        "intensity_": 1350.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "797",
        "pos_": 698.2573,
        "monoMz_": 698.2573,
        "intensity_": 1341.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "218",
        "pos_": 556.3103,
        "monoMz_": 556.3103,
        "intensity_": 1339.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "657",
        "pos_": 645.6153,
        "monoMz_": 645.6153,
        "intensity_": 1338.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "747",
        "pos_": 674.5829,
        "monoMz_": 674.5829,
        "intensity_": 1338.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "568",
        "pos_": 624.1005,
        "monoMz_": 624.1005,
        "intensity_": 1335.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "746",
        "pos_": 674.3831,
        "monoMz_": 674.3831,
        "intensity_": 1333.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "81",
        "pos_": 456.3269,
        "monoMz_": 456.3269,
        "intensity_": 1322.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "414",
        "pos_": 598.0076,
        "monoMz_": 598.0076,
        "intensity_": 1318.2
    },
    {
        "displayLevel_": 4,
        "peakId_": "850",
        "pos_": 757.6035,
        "monoMz_": 757.6035,
        "intensity_": 1306
    },
    {
        "displayLevel_": 0,
        "peakId_": "491",
        "pos_": 606.7341,
        "monoMz_": 606.7341,
        "intensity_": 1299.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "538",
        "pos_": 614.7031,
        "monoMz_": 614.7031,
        "intensity_": 1284.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "232",
        "pos_": 558.0778,
        "monoMz_": 558.0778,
        "intensity_": 1279.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "268",
        "pos_": 566.9866,
        "monoMz_": 566.9866,
        "intensity_": 1268.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "489",
        "pos_": 606.6157,
        "monoMz_": 606.6157,
        "intensity_": 1266.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "103",
        "pos_": 480.6831,
        "monoMz_": 480.6831,
        "intensity_": 1265.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "446",
        "pos_": 602.7358,
        "monoMz_": 602.7358,
        "intensity_": 1265.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "555",
        "pos_": 621.8627,
        "monoMz_": 621.8627,
        "intensity_": 1265.5
    },
    {
        "displayLevel_": 2,
        "peakId_": "194",
        "pos_": 552.0957,
        "monoMz_": 552.0957,
        "intensity_": 1265
    },
    {
        "displayLevel_": 0,
        "peakId_": "731",
        "pos_": 669.3515,
        "monoMz_": 669.3515,
        "intensity_": 1262.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "528",
        "pos_": 611.7744,
        "monoMz_": 611.7744,
        "intensity_": 1261.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "575",
        "pos_": 624.9777,
        "monoMz_": 624.9777,
        "intensity_": 1259.3
    },
    {
        "displayLevel_": 1,
        "peakId_": "382",
        "pos_": 587.6521,
        "monoMz_": 587.6521,
        "intensity_": 1251
    },
    {
        "displayLevel_": 0,
        "peakId_": "411",
        "pos_": 597.5184,
        "monoMz_": 597.5184,
        "intensity_": 1244
    },
    {
        "displayLevel_": 0,
        "peakId_": "773",
        "pos_": 688.7621,
        "monoMz_": 688.7621,
        "intensity_": 1243.8
    },
    {
        "displayLevel_": 1,
        "peakId_": "826",
        "pos_": 715.1931,
        "monoMz_": 715.1931,
        "intensity_": 1239.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "171",
        "pos_": 540.4219,
        "monoMz_": 540.4219,
        "intensity_": 1229.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "630",
        "pos_": 635.116,
        "monoMz_": 635.116,
        "intensity_": 1229.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "496",
        "pos_": 607.0745,
        "monoMz_": 607.0745,
        "intensity_": 1228
    },
    {
        "displayLevel_": 0,
        "peakId_": "482",
        "pos_": 606.0965,
        "monoMz_": 606.0965,
        "intensity_": 1224.7
    },
    {
        "displayLevel_": 3,
        "peakId_": "854",
        "pos_": 761.206,
        "monoMz_": 761.206,
        "intensity_": 1222.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "203",
        "pos_": 553.4164,
        "monoMz_": 553.4164,
        "intensity_": 1217.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "445",
        "pos_": 602.5869,
        "monoMz_": 602.5869,
        "intensity_": 1210.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "527",
        "pos_": 611.7267,
        "monoMz_": 611.7267,
        "intensity_": 1201.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "167",
        "pos_": 532.4306,
        "monoMz_": 532.4306,
        "intensity_": 1199.6
    },
    {
        "displayLevel_": 1,
        "peakId_": "645",
        "pos_": 638.6809,
        "monoMz_": 638.6809,
        "intensity_": 1199.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "681",
        "pos_": 649.7667,
        "monoMz_": 649.7667,
        "intensity_": 1198.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "683",
        "pos_": 651.5585,
        "monoMz_": 651.5585,
        "intensity_": 1192.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "676",
        "pos_": 649.2657,
        "monoMz_": 649.2657,
        "intensity_": 1190.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "863",
        "pos_": 769.1275,
        "monoMz_": 769.1275,
        "intensity_": 1186.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "912",
        "pos_": 1876.4398,
        "monoMz_": 1876.4398,
        "intensity_": 1186.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "236",
        "pos_": 558.4116,
        "monoMz_": 558.4116,
        "intensity_": 1183.5
    },
    {
        "displayLevel_": 2,
        "peakId_": "263",
        "pos_": 565.3245,
        "monoMz_": 565.3245,
        "intensity_": 1176
    },
    {
        "displayLevel_": 0,
        "peakId_": "374",
        "pos_": 584.8414,
        "monoMz_": 584.8414,
        "intensity_": 1175.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "587",
        "pos_": 626.0203,
        "monoMz_": 626.0203,
        "intensity_": 1171.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "532",
        "pos_": 613.6059,
        "monoMz_": 613.6059,
        "intensity_": 1170.5
    },
    {
        "displayLevel_": 2,
        "peakId_": "395",
        "pos_": 591.0774,
        "monoMz_": 591.0774,
        "intensity_": 1166.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "844",
        "pos_": 741.1949,
        "monoMz_": 741.1949,
        "intensity_": 1161.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "179",
        "pos_": 543.4237,
        "monoMz_": 543.4237,
        "intensity_": 1159.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "809",
        "pos_": 710.5068,
        "monoMz_": 710.5068,
        "intensity_": 1155
    },
    {
        "displayLevel_": 0,
        "peakId_": "19",
        "pos_": 199.1493,
        "monoMz_": 199.1493,
        "intensity_": 1153.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "461",
        "pos_": 604.2182,
        "monoMz_": 604.2182,
        "intensity_": 1153
    },
    {
        "displayLevel_": 0,
        "peakId_": "452",
        "pos_": 603.3359,
        "monoMz_": 603.3359,
        "intensity_": 1151.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "229",
        "pos_": 557.8543,
        "monoMz_": 557.8543,
        "intensity_": 1143.1
    },
    {
        "displayLevel_": 2,
        "peakId_": "807",
        "pos_": 710.2861,
        "monoMz_": 710.2861,
        "intensity_": 1142
    },
    {
        "displayLevel_": 0,
        "peakId_": "86",
        "pos_": 462.5901,
        "monoMz_": 462.5901,
        "intensity_": 1141.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "839",
        "pos_": 723.0339,
        "monoMz_": 723.0339,
        "intensity_": 1136.9
    },
    {
        "displayLevel_": 6,
        "peakId_": "888",
        "pos_": 1055.8302,
        "monoMz_": 1055.8302,
        "intensity_": 1133.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "146",
        "pos_": 500.2612,
        "monoMz_": 500.2612,
        "intensity_": 1133.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "604",
        "pos_": 628.7404,
        "monoMz_": 628.7404,
        "intensity_": 1132.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "117",
        "pos_": 486.9804,
        "monoMz_": 486.9804,
        "intensity_": 1129.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "283",
        "pos_": 570.529,
        "monoMz_": 570.529,
        "intensity_": 1125.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "682",
        "pos_": 651.3459,
        "monoMz_": 651.3459,
        "intensity_": 1124.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "557",
        "pos_": 622.1,
        "monoMz_": 622.1,
        "intensity_": 1123.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "319",
        "pos_": 575.5333,
        "monoMz_": 575.5333,
        "intensity_": 1115.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "302",
        "pos_": 573.8047,
        "monoMz_": 573.8047,
        "intensity_": 1110.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "644",
        "pos_": 638.6113,
        "monoMz_": 638.6113,
        "intensity_": 1110.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "717",
        "pos_": 667.0795,
        "monoMz_": 667.0795,
        "intensity_": 1106.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "643",
        "pos_": 638.4895,
        "monoMz_": 638.4895,
        "intensity_": 1102.8
    },
    {
        "displayLevel_": 2,
        "peakId_": "189",
        "pos_": 545.3162,
        "monoMz_": 545.3162,
        "intensity_": 1102
    },
    {
        "displayLevel_": 3,
        "peakId_": "68",
        "pos_": 404.5777,
        "monoMz_": 404.5777,
        "intensity_": 1100.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "533",
        "pos_": 613.73,
        "monoMz_": 613.73,
        "intensity_": 1095
    },
    {
        "displayLevel_": 0,
        "peakId_": "64",
        "pos_": 366.569,
        "monoMz_": 366.569,
        "intensity_": 1094.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "475",
        "pos_": 605.7034,
        "monoMz_": 605.7034,
        "intensity_": 1091.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "306",
        "pos_": 574.1382,
        "monoMz_": 574.1382,
        "intensity_": 1087.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "174",
        "pos_": 540.8663,
        "monoMz_": 540.8663,
        "intensity_": 1087.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "431",
        "pos_": 601.0073,
        "monoMz_": 601.0073,
        "intensity_": 1085.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "866",
        "pos_": 779.4207,
        "monoMz_": 779.4207,
        "intensity_": 1084.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "921",
        "pos_": 1905.2799,
        "monoMz_": 1905.2799,
        "intensity_": 1082.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "157",
        "pos_": 506.3892,
        "monoMz_": 506.3892,
        "intensity_": 1080.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "223",
        "pos_": 557.112,
        "monoMz_": 557.112,
        "intensity_": 1078.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "89",
        "pos_": 466.7597,
        "monoMz_": 466.7597,
        "intensity_": 1078
    },
    {
        "displayLevel_": 0,
        "peakId_": "565",
        "pos_": 623.0163,
        "monoMz_": 623.0163,
        "intensity_": 1075.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "693",
        "pos_": 655.5499,
        "monoMz_": 655.5499,
        "intensity_": 1074.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "686",
        "pos_": 653.5569,
        "monoMz_": 653.5569,
        "intensity_": 1073.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "351",
        "pos_": 578.9365,
        "monoMz_": 578.9365,
        "intensity_": 1071.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "134",
        "pos_": 497.6224,
        "monoMz_": 497.6224,
        "intensity_": 1070.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "132",
        "pos_": 494.2723,
        "monoMz_": 494.2723,
        "intensity_": 1064.7
    },
    {
        "displayLevel_": 4,
        "peakId_": "878",
        "pos_": 864.9756,
        "monoMz_": 864.9756,
        "intensity_": 1060
    },
    {
        "displayLevel_": 0,
        "peakId_": "139",
        "pos_": 499.5837,
        "monoMz_": 499.5837,
        "intensity_": 1058.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "694",
        "pos_": 655.7515,
        "monoMz_": 655.7515,
        "intensity_": 1057.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "413",
        "pos_": 597.7859,
        "monoMz_": 597.7859,
        "intensity_": 1057
    },
    {
        "displayLevel_": 0,
        "peakId_": "72",
        "pos_": 414.5768,
        "monoMz_": 414.5768,
        "intensity_": 1051
    },
    {
        "displayLevel_": 0,
        "peakId_": "471",
        "pos_": 605.2513,
        "monoMz_": 605.2513,
        "intensity_": 1049.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "440",
        "pos_": 602.3288,
        "monoMz_": 602.3288,
        "intensity_": 1046.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "778",
        "pos_": 691.8862,
        "monoMz_": 691.8862,
        "intensity_": 1045.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "820",
        "pos_": 714.3676,
        "monoMz_": 714.3676,
        "intensity_": 1044.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "846",
        "pos_": 741.5295,
        "monoMz_": 741.5295,
        "intensity_": 1043.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "164",
        "pos_": 508.2945,
        "monoMz_": 508.2945,
        "intensity_": 1037.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "799",
        "pos_": 698.5427,
        "monoMz_": 698.5427,
        "intensity_": 1035
    },
    {
        "displayLevel_": 0,
        "peakId_": "25",
        "pos_": 199.2534,
        "monoMz_": 199.2534,
        "intensity_": 1024.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "478",
        "pos_": 605.8506,
        "monoMz_": 605.8506,
        "intensity_": 1019.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "774",
        "pos_": 688.9649,
        "monoMz_": 688.9649,
        "intensity_": 1016.4
    },
    {
        "displayLevel_": 5,
        "peakId_": "871",
        "pos_": 815.8806,
        "monoMz_": 815.8806,
        "intensity_": 1012.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "92",
        "pos_": 467.2641,
        "monoMz_": 467.2641,
        "intensity_": 1011.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "484",
        "pos_": 606.2173,
        "monoMz_": 606.2173,
        "intensity_": 1011
    },
    {
        "displayLevel_": 0,
        "peakId_": "819",
        "pos_": 714.1196,
        "monoMz_": 714.1196,
        "intensity_": 1010.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "544",
        "pos_": 615.3292,
        "monoMz_": 615.3292,
        "intensity_": 1010.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "554",
        "pos_": 621.7155,
        "monoMz_": 621.7155,
        "intensity_": 1008
    },
    {
        "displayLevel_": 0,
        "peakId_": "825",
        "pos_": 715.0223,
        "monoMz_": 715.0223,
        "intensity_": 1006
    },
    {
        "displayLevel_": 0,
        "peakId_": "121",
        "pos_": 488.7754,
        "monoMz_": 488.7754,
        "intensity_": 1004
    },
    {
        "displayLevel_": 0,
        "peakId_": "221",
        "pos_": 556.8027,
        "monoMz_": 556.8027,
        "intensity_": 1002.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "191",
        "pos_": 549.9821,
        "monoMz_": 549.9821,
        "intensity_": 998.38
    },
    {
        "displayLevel_": 0,
        "peakId_": "181",
        "pos_": 543.5758,
        "monoMz_": 543.5758,
        "intensity_": 998.24
    },
    {
        "displayLevel_": 0,
        "peakId_": "116",
        "pos_": 486.8813,
        "monoMz_": 486.8813,
        "intensity_": 996.16
    },
    {
        "displayLevel_": 1,
        "peakId_": "806",
        "pos_": 709.5487,
        "monoMz_": 709.5487,
        "intensity_": 993.26
    },
    {
        "displayLevel_": 0,
        "peakId_": "416",
        "pos_": 598.4842,
        "monoMz_": 598.4842,
        "intensity_": 992.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "62",
        "pos_": 366.2345,
        "monoMz_": 366.2345,
        "intensity_": 991.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "762",
        "pos_": 683.7862,
        "monoMz_": 683.7862,
        "intensity_": 988.18
    },
    {
        "displayLevel_": 0,
        "peakId_": "521",
        "pos_": 610.0073,
        "monoMz_": 610.0073,
        "intensity_": 983.72
    },
    {
        "displayLevel_": 0,
        "peakId_": "320",
        "pos_": 575.6448,
        "monoMz_": 575.6448,
        "intensity_": 981.82
    },
    {
        "displayLevel_": 0,
        "peakId_": "721",
        "pos_": 667.4773,
        "monoMz_": 667.4773,
        "intensity_": 979.79
    },
    {
        "displayLevel_": 0,
        "peakId_": "341",
        "pos_": 577.6432,
        "monoMz_": 577.6432,
        "intensity_": 973.67
    },
    {
        "displayLevel_": 0,
        "peakId_": "779",
        "pos_": 691.9664,
        "monoMz_": 691.9664,
        "intensity_": 973.25
    },
    {
        "displayLevel_": 0,
        "peakId_": "202",
        "pos_": 553.1984,
        "monoMz_": 553.1984,
        "intensity_": 972.21
    },
    {
        "displayLevel_": 0,
        "peakId_": "115",
        "pos_": 486.7805,
        "monoMz_": 486.7805,
        "intensity_": 969.98
    },
    {
        "displayLevel_": 0,
        "peakId_": "856",
        "pos_": 766.044,
        "monoMz_": 766.044,
        "intensity_": 969.82
    },
    {
        "displayLevel_": 0,
        "peakId_": "796",
        "pos_": 698.1122,
        "monoMz_": 698.1122,
        "intensity_": 968.97
    },
    {
        "displayLevel_": 0,
        "peakId_": "560",
        "pos_": 622.4785,
        "monoMz_": 622.4785,
        "intensity_": 968.51
    },
    {
        "displayLevel_": 0,
        "peakId_": "559",
        "pos_": 622.3448,
        "monoMz_": 622.3448,
        "intensity_": 967.92
    },
    {
        "displayLevel_": 0,
        "peakId_": "658",
        "pos_": 645.7382,
        "monoMz_": 645.7382,
        "intensity_": 964.69
    },
    {
        "displayLevel_": 0,
        "peakId_": "233",
        "pos_": 558.1289,
        "monoMz_": 558.1289,
        "intensity_": 964.12
    },
    {
        "displayLevel_": 0,
        "peakId_": "285",
        "pos_": 570.7502,
        "monoMz_": 570.7502,
        "intensity_": 959.59
    },
    {
        "displayLevel_": 0,
        "peakId_": "418",
        "pos_": 598.6508,
        "monoMz_": 598.6508,
        "intensity_": 958.94
    },
    {
        "displayLevel_": 0,
        "peakId_": "539",
        "pos_": 614.8273,
        "monoMz_": 614.8273,
        "intensity_": 954.18
    },
    {
        "displayLevel_": 2,
        "peakId_": "880",
        "pos_": 872.4261,
        "monoMz_": 872.4261,
        "intensity_": 949.68
    },
    {
        "displayLevel_": 0,
        "peakId_": "172",
        "pos_": 540.5312,
        "monoMz_": 540.5312,
        "intensity_": 946.45
    },
    {
        "displayLevel_": 0,
        "peakId_": "661",
        "pos_": 646.4702,
        "monoMz_": 646.4702,
        "intensity_": 944.71
    },
    {
        "displayLevel_": 1,
        "peakId_": "254",
        "pos_": 563.2024,
        "monoMz_": 563.2024,
        "intensity_": 944.09
    },
    {
        "displayLevel_": 0,
        "peakId_": "841",
        "pos_": 724.2123,
        "monoMz_": 724.2123,
        "intensity_": 944.08
    },
    {
        "displayLevel_": 0,
        "peakId_": "299",
        "pos_": 573.1973,
        "monoMz_": 573.1973,
        "intensity_": 937.89
    },
    {
        "displayLevel_": 0,
        "peakId_": "883",
        "pos_": 876.7471,
        "monoMz_": 876.7471,
        "intensity_": 937.51
    },
    {
        "displayLevel_": 0,
        "peakId_": "380",
        "pos_": 587.4301,
        "monoMz_": 587.4301,
        "intensity_": 934.22
    },
    {
        "displayLevel_": 0,
        "peakId_": "800",
        "pos_": 698.826,
        "monoMz_": 698.826,
        "intensity_": 934.14
    },
    {
        "displayLevel_": 0,
        "peakId_": "540",
        "pos_": 614.8985,
        "monoMz_": 614.8985,
        "intensity_": 933.81
    },
    {
        "displayLevel_": 0,
        "peakId_": "822",
        "pos_": 714.5246,
        "monoMz_": 714.5246,
        "intensity_": 932.89
    },
    {
        "displayLevel_": 0,
        "peakId_": "651",
        "pos_": 642.9312,
        "monoMz_": 642.9312,
        "intensity_": 932.17
    },
    {
        "displayLevel_": 0,
        "peakId_": "824",
        "pos_": 714.8581,
        "monoMz_": 714.8581,
        "intensity_": 928.99
    },
    {
        "displayLevel_": 0,
        "peakId_": "659",
        "pos_": 645.8644,
        "monoMz_": 645.8644,
        "intensity_": 928.65
    },
    {
        "displayLevel_": 0,
        "peakId_": "377",
        "pos_": 586.6749,
        "monoMz_": 586.6749,
        "intensity_": 925.96
    },
    {
        "displayLevel_": 0,
        "peakId_": "464",
        "pos_": 604.436,
        "monoMz_": 604.436,
        "intensity_": 924.66
    },
    {
        "displayLevel_": 0,
        "peakId_": "695",
        "pos_": 656.7568,
        "monoMz_": 656.7568,
        "intensity_": 923.76
    },
    {
        "displayLevel_": 0,
        "peakId_": "196",
        "pos_": 552.3159,
        "monoMz_": 552.3159,
        "intensity_": 922.38
    },
    {
        "displayLevel_": 0,
        "peakId_": "161",
        "pos_": 506.7932,
        "monoMz_": 506.7932,
        "intensity_": 921.91
    },
    {
        "displayLevel_": 0,
        "peakId_": "87",
        "pos_": 464.7208,
        "monoMz_": 464.7208,
        "intensity_": 920.88
    },
    {
        "displayLevel_": 0,
        "peakId_": "744",
        "pos_": 673.9825,
        "monoMz_": 673.9825,
        "intensity_": 919.71
    },
    {
        "displayLevel_": 0,
        "peakId_": "104",
        "pos_": 480.8251,
        "monoMz_": 480.8251,
        "intensity_": 919.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "519",
        "pos_": 609.6721,
        "monoMz_": 609.6721,
        "intensity_": 916.35
    },
    {
        "displayLevel_": 0,
        "peakId_": "669",
        "pos_": 648.1682,
        "monoMz_": 648.1682,
        "intensity_": 913.32
    },
    {
        "displayLevel_": 0,
        "peakId_": "26",
        "pos_": 199.2627,
        "monoMz_": 199.2627,
        "intensity_": 912.86
    },
    {
        "displayLevel_": 5,
        "peakId_": "44",
        "pos_": 292.7136,
        "monoMz_": 292.7136,
        "intensity_": 912.05
    },
    {
        "displayLevel_": 4,
        "peakId_": "873",
        "pos_": 840.1907,
        "monoMz_": 840.1907,
        "intensity_": 909.68
    },
    {
        "displayLevel_": 0,
        "peakId_": "308",
        "pos_": 574.314,
        "monoMz_": 574.314,
        "intensity_": 908.63
    },
    {
        "displayLevel_": 0,
        "peakId_": "347",
        "pos_": 578.4617,
        "monoMz_": 578.4617,
        "intensity_": 908.48
    },
    {
        "displayLevel_": 0,
        "peakId_": "803",
        "pos_": 708.1755,
        "monoMz_": 708.1755,
        "intensity_": 903.16
    },
    {
        "displayLevel_": 2,
        "peakId_": "843",
        "pos_": 737.9929,
        "monoMz_": 737.9929,
        "intensity_": 903.04
    },
    {
        "displayLevel_": 0,
        "peakId_": "881",
        "pos_": 872.924,
        "monoMz_": 872.924,
        "intensity_": 901.68
    },
    {
        "displayLevel_": 0,
        "peakId_": "656",
        "pos_": 644.7192,
        "monoMz_": 644.7192,
        "intensity_": 900.1
    },
    {
        "displayLevel_": 1,
        "peakId_": "409",
        "pos_": 596.3328,
        "monoMz_": 596.3328,
        "intensity_": 898.26
    },
    {
        "displayLevel_": 0,
        "peakId_": "412",
        "pos_": 597.6098,
        "monoMz_": 597.6098,
        "intensity_": 893.82
    },
    {
        "displayLevel_": 0,
        "peakId_": "441",
        "pos_": 602.3518,
        "monoMz_": 602.3518,
        "intensity_": 892.97
    },
    {
        "displayLevel_": 0,
        "peakId_": "408",
        "pos_": 596.0986,
        "monoMz_": 596.0986,
        "intensity_": 892.48
    },
    {
        "displayLevel_": 2,
        "peakId_": "855",
        "pos_": 763.7203,
        "monoMz_": 763.7203,
        "intensity_": 892.4
    },
    {
        "displayLevel_": 3,
        "peakId_": "875",
        "pos_": 855.4632,
        "monoMz_": 855.4632,
        "intensity_": 890.52
    },
    {
        "displayLevel_": 0,
        "peakId_": "193",
        "pos_": 551.9835,
        "monoMz_": 551.9835,
        "intensity_": 889.15
    },
    {
        "displayLevel_": 0,
        "peakId_": "628",
        "pos_": 634.8967,
        "monoMz_": 634.8967,
        "intensity_": 886.87
    },
    {
        "displayLevel_": 0,
        "peakId_": "386",
        "pos_": 589.2052,
        "monoMz_": 589.2052,
        "intensity_": 885.97
    },
    {
        "displayLevel_": 0,
        "peakId_": "680",
        "pos_": 649.6646,
        "monoMz_": 649.6646,
        "intensity_": 884.51
    },
    {
        "displayLevel_": 0,
        "peakId_": "642",
        "pos_": 637.7598,
        "monoMz_": 637.7598,
        "intensity_": 883.21
    },
    {
        "displayLevel_": 0,
        "peakId_": "165",
        "pos_": 508.3956,
        "monoMz_": 508.3956,
        "intensity_": 882.66
    },
    {
        "displayLevel_": 1,
        "peakId_": "709",
        "pos_": 663.3633,
        "monoMz_": 663.3633,
        "intensity_": 877.55
    },
    {
        "displayLevel_": 0,
        "peakId_": "444",
        "pos_": 602.5464,
        "monoMz_": 602.5464,
        "intensity_": 875.41
    },
    {
        "displayLevel_": 0,
        "peakId_": "406",
        "pos_": 595.8488,
        "monoMz_": 595.8488,
        "intensity_": 874.97
    },
    {
        "displayLevel_": 0,
        "peakId_": "269",
        "pos_": 567.3195,
        "monoMz_": 567.3195,
        "intensity_": 874.76
    },
    {
        "displayLevel_": 0,
        "peakId_": "176",
        "pos_": 541.2914,
        "monoMz_": 541.2914,
        "intensity_": 874.72
    },
    {
        "displayLevel_": 0,
        "peakId_": "147",
        "pos_": 502.278,
        "monoMz_": 502.278,
        "intensity_": 872.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "582",
        "pos_": 625.5844,
        "monoMz_": 625.5844,
        "intensity_": 872.12
    },
    {
        "displayLevel_": 0,
        "peakId_": "913",
        "pos_": 1876.6373,
        "monoMz_": 1876.6373,
        "intensity_": 869.12
    },
    {
        "displayLevel_": 5,
        "peakId_": "899",
        "pos_": 1654.0144,
        "monoMz_": 1654.0144,
        "intensity_": 868.49
    },
    {
        "displayLevel_": 0,
        "peakId_": "28",
        "pos_": 199.345,
        "monoMz_": 199.345,
        "intensity_": 866.89
    },
    {
        "displayLevel_": 4,
        "peakId_": "900",
        "pos_": 1707.3395,
        "monoMz_": 1707.3395,
        "intensity_": 866.55
    },
    {
        "displayLevel_": 3,
        "peakId_": "97",
        "pos_": 478.1113,
        "monoMz_": 478.1113,
        "intensity_": 866.38
    },
    {
        "displayLevel_": 0,
        "peakId_": "437",
        "pos_": 601.9869,
        "monoMz_": 601.9869,
        "intensity_": 866.05
    },
    {
        "displayLevel_": 0,
        "peakId_": "451",
        "pos_": 603.2942,
        "monoMz_": 603.2942,
        "intensity_": 865.67
    },
    {
        "displayLevel_": 0,
        "peakId_": "783",
        "pos_": 692.5708,
        "monoMz_": 692.5708,
        "intensity_": 864.91
    },
    {
        "displayLevel_": 0,
        "peakId_": "882",
        "pos_": 876.4098,
        "monoMz_": 876.4098,
        "intensity_": 864.76
    },
    {
        "displayLevel_": 0,
        "peakId_": "410",
        "pos_": 597.1442,
        "monoMz_": 597.1442,
        "intensity_": 859.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "561",
        "pos_": 622.683,
        "monoMz_": 622.683,
        "intensity_": 859.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "258",
        "pos_": 564.3101,
        "monoMz_": 564.3101,
        "intensity_": 856.33
    },
    {
        "displayLevel_": 0,
        "peakId_": "257",
        "pos_": 563.8194,
        "monoMz_": 563.8194,
        "intensity_": 854.18
    },
    {
        "displayLevel_": 0,
        "peakId_": "190",
        "pos_": 549.0889,
        "monoMz_": 549.0889,
        "intensity_": 853.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "634",
        "pos_": 636.3887,
        "monoMz_": 636.3887,
        "intensity_": 852.16
    },
    {
        "displayLevel_": 0,
        "peakId_": "267",
        "pos_": 566.6237,
        "monoMz_": 566.6237,
        "intensity_": 850.41
    },
    {
        "displayLevel_": 0,
        "peakId_": "752",
        "pos_": 676.9005,
        "monoMz_": 676.9005,
        "intensity_": 848.62
    },
    {
        "displayLevel_": 0,
        "peakId_": "457",
        "pos_": 603.9695,
        "monoMz_": 603.9695,
        "intensity_": 845.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "487",
        "pos_": 606.518,
        "monoMz_": 606.518,
        "intensity_": 839.88
    },
    {
        "displayLevel_": 0,
        "peakId_": "91",
        "pos_": 467.097,
        "monoMz_": 467.097,
        "intensity_": 835.22
    },
    {
        "displayLevel_": 0,
        "peakId_": "23",
        "pos_": 199.2431,
        "monoMz_": 199.2431,
        "intensity_": 831.17
    },
    {
        "displayLevel_": 0,
        "peakId_": "136",
        "pos_": 497.9834,
        "monoMz_": 497.9834,
        "intensity_": 830.67
    },
    {
        "displayLevel_": 0,
        "peakId_": "518",
        "pos_": 609.5887,
        "monoMz_": 609.5887,
        "intensity_": 830.22
    },
    {
        "displayLevel_": 0,
        "peakId_": "270",
        "pos_": 567.7769,
        "monoMz_": 567.7769,
        "intensity_": 829.49
    },
    {
        "displayLevel_": 0,
        "peakId_": "300",
        "pos_": 573.4211,
        "monoMz_": 573.4211,
        "intensity_": 829.48
    },
    {
        "displayLevel_": 0,
        "peakId_": "546",
        "pos_": 619.3473,
        "monoMz_": 619.3473,
        "intensity_": 829.34
    },
    {
        "displayLevel_": 0,
        "peakId_": "314",
        "pos_": 575.0333,
        "monoMz_": 575.0333,
        "intensity_": 829.06
    },
    {
        "displayLevel_": 0,
        "peakId_": "685",
        "pos_": 653.3401,
        "monoMz_": 653.3401,
        "intensity_": 828.19
    },
    {
        "displayLevel_": 0,
        "peakId_": "543",
        "pos_": 615.2022,
        "monoMz_": 615.2022,
        "intensity_": 824.02
    },
    {
        "displayLevel_": 5,
        "peakId_": "48",
        "pos_": 308.8542,
        "monoMz_": 308.8542,
        "intensity_": 823.46
    },
    {
        "displayLevel_": 2,
        "peakId_": "867",
        "pos_": 781.431,
        "monoMz_": 781.431,
        "intensity_": 823.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "468",
        "pos_": 604.8505,
        "monoMz_": 604.8505,
        "intensity_": 822.23
    },
    {
        "displayLevel_": 0,
        "peakId_": "725",
        "pos_": 668.8773,
        "monoMz_": 668.8773,
        "intensity_": 821.13
    },
    {
        "displayLevel_": 0,
        "peakId_": "153",
        "pos_": 505.5292,
        "monoMz_": 505.5292,
        "intensity_": 818.19
    },
    {
        "displayLevel_": 0,
        "peakId_": "529",
        "pos_": 613.2286,
        "monoMz_": 613.2286,
        "intensity_": 816.12
    },
    {
        "displayLevel_": 0,
        "peakId_": "853",
        "pos_": 761.0085,
        "monoMz_": 761.0085,
        "intensity_": 813.64
    },
    {
        "displayLevel_": 0,
        "peakId_": "160",
        "pos_": 506.6886,
        "monoMz_": 506.6886,
        "intensity_": 811.77
    },
    {
        "displayLevel_": 0,
        "peakId_": "473",
        "pos_": 605.528,
        "monoMz_": 605.528,
        "intensity_": 811.76
    },
    {
        "displayLevel_": 2,
        "peakId_": "849",
        "pos_": 754.1694,
        "monoMz_": 754.1694,
        "intensity_": 808.26
    },
    {
        "displayLevel_": 0,
        "peakId_": "379",
        "pos_": 587.3402,
        "monoMz_": 587.3402,
        "intensity_": 807.09
    },
    {
        "displayLevel_": 0,
        "peakId_": "734",
        "pos_": 669.7576,
        "monoMz_": 669.7576,
        "intensity_": 806.43
    },
    {
        "displayLevel_": 0,
        "peakId_": "566",
        "pos_": 623.4595,
        "monoMz_": 623.4595,
        "intensity_": 804.92
    },
    {
        "displayLevel_": 0,
        "peakId_": "272",
        "pos_": 567.9895,
        "monoMz_": 567.9895,
        "intensity_": 802.96
    },
    {
        "displayLevel_": 0,
        "peakId_": "479",
        "pos_": 605.8897,
        "monoMz_": 605.8897,
        "intensity_": 802.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "667",
        "pos_": 647.6652,
        "monoMz_": 647.6652,
        "intensity_": 801.87
    },
    {
        "displayLevel_": 0,
        "peakId_": "407",
        "pos_": 596.0029,
        "monoMz_": 596.0029,
        "intensity_": 801.77
    },
    {
        "displayLevel_": 0,
        "peakId_": "184",
        "pos_": 543.7598,
        "monoMz_": 543.7598,
        "intensity_": 798.82
    },
    {
        "displayLevel_": 0,
        "peakId_": "498",
        "pos_": 607.4625,
        "monoMz_": 607.4625,
        "intensity_": 797.92
    },
    {
        "displayLevel_": 0,
        "peakId_": "158",
        "pos_": 506.4905,
        "monoMz_": 506.4905,
        "intensity_": 796.11
    },
    {
        "displayLevel_": 0,
        "peakId_": "343",
        "pos_": 578.0869,
        "monoMz_": 578.0869,
        "intensity_": 794.38
    },
    {
        "displayLevel_": 0,
        "peakId_": "758",
        "pos_": 682.1915,
        "monoMz_": 682.1915,
        "intensity_": 794.16
    },
    {
        "displayLevel_": 0,
        "peakId_": "649",
        "pos_": 642.0966,
        "monoMz_": 642.0966,
        "intensity_": 791.69
    },
    {
        "displayLevel_": 0,
        "peakId_": "394",
        "pos_": 590.425,
        "monoMz_": 590.425,
        "intensity_": 790.94
    },
    {
        "displayLevel_": 0,
        "peakId_": "422",
        "pos_": 599.3242,
        "monoMz_": 599.3242,
        "intensity_": 787.47
    },
    {
        "displayLevel_": 0,
        "peakId_": "810",
        "pos_": 711.1174,
        "monoMz_": 711.1174,
        "intensity_": 786.98
    },
    {
        "displayLevel_": 0,
        "peakId_": "760",
        "pos_": 683.0826,
        "monoMz_": 683.0826,
        "intensity_": 784.73
    },
    {
        "displayLevel_": 3,
        "peakId_": "908",
        "pos_": 1725.9979,
        "monoMz_": 1725.9979,
        "intensity_": 783.59
    },
    {
        "displayLevel_": 0,
        "peakId_": "805",
        "pos_": 708.5084,
        "monoMz_": 708.5084,
        "intensity_": 783.36
    },
    {
        "displayLevel_": 0,
        "peakId_": "631",
        "pos_": 635.8611,
        "monoMz_": 635.8611,
        "intensity_": 783.07
    },
    {
        "displayLevel_": 5,
        "peakId_": "922",
        "pos_": 1930.2637,
        "monoMz_": 1930.2637,
        "intensity_": 782.02
    },
    {
        "displayLevel_": 0,
        "peakId_": "765",
        "pos_": 683.9871,
        "monoMz_": 683.9871,
        "intensity_": 781.26
    },
    {
        "displayLevel_": 0,
        "peakId_": "118",
        "pos_": 487.0838,
        "monoMz_": 487.0838,
        "intensity_": 780.47
    },
    {
        "displayLevel_": 0,
        "peakId_": "427",
        "pos_": 600.5734,
        "monoMz_": 600.5734,
        "intensity_": 779.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "504",
        "pos_": 608.0539,
        "monoMz_": 608.0539,
        "intensity_": 778.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "767",
        "pos_": 684.3821,
        "monoMz_": 684.3821,
        "intensity_": 778.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "531",
        "pos_": 613.4779,
        "monoMz_": 613.4779,
        "intensity_": 778.26
    },
    {
        "displayLevel_": 0,
        "peakId_": "280",
        "pos_": 570.0854,
        "monoMz_": 570.0854,
        "intensity_": 778.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "65",
        "pos_": 366.7692,
        "monoMz_": 366.7692,
        "intensity_": 774.06
    },
    {
        "displayLevel_": 2,
        "peakId_": "96",
        "pos_": 474.3365,
        "monoMz_": 474.3365,
        "intensity_": 772.96
    },
    {
        "displayLevel_": 0,
        "peakId_": "98",
        "pos_": 480.2696,
        "monoMz_": 480.2696,
        "intensity_": 772.91
    },
    {
        "displayLevel_": 0,
        "peakId_": "821",
        "pos_": 714.4002,
        "monoMz_": 714.4002,
        "intensity_": 770.22
    },
    {
        "displayLevel_": 0,
        "peakId_": "917",
        "pos_": 1877.6901,
        "monoMz_": 1877.6901,
        "intensity_": 766.33
    },
    {
        "displayLevel_": 4,
        "peakId_": "42",
        "pos_": 280.1518,
        "monoMz_": 280.1518,
        "intensity_": 765.37
    },
    {
        "displayLevel_": 0,
        "peakId_": "170",
        "pos_": 540.3135,
        "monoMz_": 540.3135,
        "intensity_": 765.15
    },
    {
        "displayLevel_": 0,
        "peakId_": "124",
        "pos_": 489.7825,
        "monoMz_": 489.7825,
        "intensity_": 764.74
    },
    {
        "displayLevel_": 0,
        "peakId_": "633",
        "pos_": 636.3454,
        "monoMz_": 636.3454,
        "intensity_": 763.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "603",
        "pos_": 628.6053,
        "monoMz_": 628.6053,
        "intensity_": 763
    },
    {
        "displayLevel_": 0,
        "peakId_": "450",
        "pos_": 603.1344,
        "monoMz_": 603.1344,
        "intensity_": 762.74
    },
    {
        "displayLevel_": 0,
        "peakId_": "652",
        "pos_": 643.0551,
        "monoMz_": 643.0551,
        "intensity_": 759.41
    },
    {
        "displayLevel_": 0,
        "peakId_": "576",
        "pos_": 625.1027,
        "monoMz_": 625.1027,
        "intensity_": 756
    },
    {
        "displayLevel_": 0,
        "peakId_": "261",
        "pos_": 565.124,
        "monoMz_": 565.124,
        "intensity_": 754.74
    },
    {
        "displayLevel_": 0,
        "peakId_": "271",
        "pos_": 567.8808,
        "monoMz_": 567.8808,
        "intensity_": 754.48
    },
    {
        "displayLevel_": 0,
        "peakId_": "552",
        "pos_": 621.2301,
        "monoMz_": 621.2301,
        "intensity_": 752.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "396",
        "pos_": 591.3314,
        "monoMz_": 591.3314,
        "intensity_": 745.18
    },
    {
        "displayLevel_": 0,
        "peakId_": "127",
        "pos_": 493.2284,
        "monoMz_": 493.2284,
        "intensity_": 744.37
    },
    {
        "displayLevel_": 0,
        "peakId_": "195",
        "pos_": 552.2082,
        "monoMz_": 552.2082,
        "intensity_": 744.28
    },
    {
        "displayLevel_": 3,
        "peakId_": "55",
        "pos_": 348.7478,
        "monoMz_": 348.7478,
        "intensity_": 743.46
    },
    {
        "displayLevel_": 0,
        "peakId_": "919",
        "pos_": 1904.7081,
        "monoMz_": 1904.7081,
        "intensity_": 743.24
    },
    {
        "displayLevel_": 0,
        "peakId_": "419",
        "pos_": 598.697,
        "monoMz_": 598.697,
        "intensity_": 743.22
    },
    {
        "displayLevel_": 0,
        "peakId_": "466",
        "pos_": 604.512,
        "monoMz_": 604.512,
        "intensity_": 742.43
    },
    {
        "displayLevel_": 0,
        "peakId_": "499",
        "pos_": 607.5937,
        "monoMz_": 607.5937,
        "intensity_": 740.15
    },
    {
        "displayLevel_": 0,
        "peakId_": "495",
        "pos_": 607.0426,
        "monoMz_": 607.0426,
        "intensity_": 737.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "352",
        "pos_": 579.4562,
        "monoMz_": 579.4562,
        "intensity_": 737.32
    },
    {
        "displayLevel_": 0,
        "peakId_": "119",
        "pos_": 487.9571,
        "monoMz_": 487.9571,
        "intensity_": 736.25
    },
    {
        "displayLevel_": 0,
        "peakId_": "335",
        "pos_": 577.1766,
        "monoMz_": 577.1766,
        "intensity_": 736.09
    },
    {
        "displayLevel_": 5,
        "peakId_": "889",
        "pos_": 1079.5471,
        "monoMz_": 1079.5471,
        "intensity_": 735.39
    },
    {
        "displayLevel_": 0,
        "peakId_": "307",
        "pos_": 574.2185,
        "monoMz_": 574.2185,
        "intensity_": 735.09
    },
    {
        "displayLevel_": 0,
        "peakId_": "20",
        "pos_": 199.2238,
        "monoMz_": 199.2238,
        "intensity_": 733.24
    },
    {
        "displayLevel_": 0,
        "peakId_": "198",
        "pos_": 552.6376,
        "monoMz_": 552.6376,
        "intensity_": 732.17
    },
    {
        "displayLevel_": 0,
        "peakId_": "260",
        "pos_": 564.8198,
        "monoMz_": 564.8198,
        "intensity_": 731.78
    },
    {
        "displayLevel_": 0,
        "peakId_": "617",
        "pos_": 633.0089,
        "monoMz_": 633.0089,
        "intensity_": 731.47
    },
    {
        "displayLevel_": 0,
        "peakId_": "732",
        "pos_": 669.5528,
        "monoMz_": 669.5528,
        "intensity_": 729.89
    },
    {
        "displayLevel_": 0,
        "peakId_": "205",
        "pos_": 554.6257,
        "monoMz_": 554.6257,
        "intensity_": 729.13
    },
    {
        "displayLevel_": 5,
        "peakId_": "923",
        "pos_": 1962.7487,
        "monoMz_": 1962.7487,
        "intensity_": 728.55
    },
    {
        "displayLevel_": 0,
        "peakId_": "804",
        "pos_": 708.2803,
        "monoMz_": 708.2803,
        "intensity_": 727.81
    },
    {
        "displayLevel_": 0,
        "peakId_": "292",
        "pos_": 572.0797,
        "monoMz_": 572.0797,
        "intensity_": 727.51
    },
    {
        "displayLevel_": 0,
        "peakId_": "282",
        "pos_": 570.4193,
        "monoMz_": 570.4193,
        "intensity_": 726.61
    },
    {
        "displayLevel_": 0,
        "peakId_": "123",
        "pos_": 489.2761,
        "monoMz_": 489.2761,
        "intensity_": 724.19
    },
    {
        "displayLevel_": 0,
        "peakId_": "381",
        "pos_": 587.54,
        "monoMz_": 587.54,
        "intensity_": 723.79
    },
    {
        "displayLevel_": 0,
        "peakId_": "376",
        "pos_": 585.4215,
        "monoMz_": 585.4215,
        "intensity_": 723.09
    },
    {
        "displayLevel_": 1,
        "peakId_": "847",
        "pos_": 743.0837,
        "monoMz_": 743.0837,
        "intensity_": 722.97
    },
    {
        "displayLevel_": 0,
        "peakId_": "375",
        "pos_": 585.0032,
        "monoMz_": 585.0032,
        "intensity_": 719.99
    },
    {
        "displayLevel_": 0,
        "peakId_": "811",
        "pos_": 713.1168,
        "monoMz_": 713.1168,
        "intensity_": 718.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "755",
        "pos_": 678.3819,
        "monoMz_": 678.3819,
        "intensity_": 718.58
    },
    {
        "displayLevel_": 0,
        "peakId_": "199",
        "pos_": 552.868,
        "monoMz_": 552.868,
        "intensity_": 718.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "138",
        "pos_": 498.2923,
        "monoMz_": 498.2923,
        "intensity_": 716
    },
    {
        "displayLevel_": 0,
        "peakId_": "874",
        "pos_": 854.9665,
        "monoMz_": 854.9665,
        "intensity_": 714.75
    },
    {
        "displayLevel_": 0,
        "peakId_": "646",
        "pos_": 638.8651,
        "monoMz_": 638.8651,
        "intensity_": 713.92
    },
    {
        "displayLevel_": 0,
        "peakId_": "660",
        "pos_": 645.999,
        "monoMz_": 645.999,
        "intensity_": 713.45
    },
    {
        "displayLevel_": 1,
        "peakId_": "851",
        "pos_": 759.6461,
        "monoMz_": 759.6461,
        "intensity_": 713.09
    },
    {
        "displayLevel_": 0,
        "peakId_": "149",
        "pos_": 505.0284,
        "monoMz_": 505.0284,
        "intensity_": 711.04
    },
    {
        "displayLevel_": 0,
        "peakId_": "287",
        "pos_": 570.9845,
        "monoMz_": 570.9845,
        "intensity_": 710.28
    },
    {
        "displayLevel_": 0,
        "peakId_": "216",
        "pos_": 556.1255,
        "monoMz_": 556.1255,
        "intensity_": 709.28
    },
    {
        "displayLevel_": 0,
        "peakId_": "424",
        "pos_": 600.0004,
        "monoMz_": 600.0004,
        "intensity_": 707.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "250",
        "pos_": 561.8192,
        "monoMz_": 561.8192,
        "intensity_": 706.99
    },
    {
        "displayLevel_": 3,
        "peakId_": "870",
        "pos_": 809.3957,
        "monoMz_": 809.3957,
        "intensity_": 706.06
    },
    {
        "displayLevel_": 0,
        "peakId_": "647",
        "pos_": 639.4575,
        "monoMz_": 639.4575,
        "intensity_": 705.6
    },
    {
        "displayLevel_": 5,
        "peakId_": "37",
        "pos_": 245.0747,
        "monoMz_": 245.0747,
        "intensity_": 704.66
    },
    {
        "displayLevel_": 0,
        "peakId_": "879",
        "pos_": 872.0254,
        "monoMz_": 872.0254,
        "intensity_": 702.92
    },
    {
        "displayLevel_": 1,
        "peakId_": "95",
        "pos_": 470.9298,
        "monoMz_": 470.9298,
        "intensity_": 702.17
    },
    {
        "displayLevel_": 0,
        "peakId_": "423",
        "pos_": 599.4301,
        "monoMz_": 599.4301,
        "intensity_": 701.08
    },
    {
        "displayLevel_": 0,
        "peakId_": "262",
        "pos_": 565.251,
        "monoMz_": 565.251,
        "intensity_": 700.28
    },
    {
        "displayLevel_": 0,
        "peakId_": "691",
        "pos_": 655.3506,
        "monoMz_": 655.3506,
        "intensity_": 699.37
    },
    {
        "displayLevel_": 0,
        "peakId_": "197",
        "pos_": 552.4269,
        "monoMz_": 552.4269,
        "intensity_": 696.26
    },
    {
        "displayLevel_": 0,
        "peakId_": "470",
        "pos_": 605.1583,
        "monoMz_": 605.1583,
        "intensity_": 696.16
    },
    {
        "displayLevel_": 0,
        "peakId_": "361",
        "pos_": 580.6504,
        "monoMz_": 580.6504,
        "intensity_": 695.88
    },
    {
        "displayLevel_": 0,
        "peakId_": "288",
        "pos_": 571.1985,
        "monoMz_": 571.1985,
        "intensity_": 694.59
    },
    {
        "displayLevel_": 0,
        "peakId_": "600",
        "pos_": 628.2092,
        "monoMz_": 628.2092,
        "intensity_": 693.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "692",
        "pos_": 655.3642,
        "monoMz_": 655.3642,
        "intensity_": 693.42
    },
    {
        "displayLevel_": 0,
        "peakId_": "342",
        "pos_": 577.8378,
        "monoMz_": 577.8378,
        "intensity_": 693.25
    },
    {
        "displayLevel_": 1,
        "peakId_": "111",
        "pos_": 484.1091,
        "monoMz_": 484.1091,
        "intensity_": 692.68
    },
    {
        "displayLevel_": 0,
        "peakId_": "791",
        "pos_": 694.9265,
        "monoMz_": 694.9265,
        "intensity_": 691.81
    },
    {
        "displayLevel_": 0,
        "peakId_": "286",
        "pos_": 570.8611,
        "monoMz_": 570.8611,
        "intensity_": 689.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "188",
        "pos_": 544.2986,
        "monoMz_": 544.2986,
        "intensity_": 686.11
    },
    {
        "displayLevel_": 6,
        "peakId_": "894",
        "pos_": 1306.8431,
        "monoMz_": 1306.8431,
        "intensity_": 685.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "356",
        "pos_": 579.8748,
        "monoMz_": 579.8748,
        "intensity_": 684.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "860",
        "pos_": 768.5773,
        "monoMz_": 768.5773,
        "intensity_": 683.94
    },
    {
        "displayLevel_": 0,
        "peakId_": "301",
        "pos_": 573.6492,
        "monoMz_": 573.6492,
        "intensity_": 683.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "852",
        "pos_": 760.8002,
        "monoMz_": 760.8002,
        "intensity_": 682.35
    },
    {
        "displayLevel_": 0,
        "peakId_": "562",
        "pos_": 622.7296,
        "monoMz_": 622.7296,
        "intensity_": 681.5
    },
    {
        "displayLevel_": 3,
        "peakId_": "848",
        "pos_": 745.4351,
        "monoMz_": 745.4351,
        "intensity_": 680.47
    },
    {
        "displayLevel_": 0,
        "peakId_": "403",
        "pos_": 595.5778,
        "monoMz_": 595.5778,
        "intensity_": 679.04
    },
    {
        "displayLevel_": 0,
        "peakId_": "404",
        "pos_": 595.6721,
        "monoMz_": 595.6721,
        "intensity_": 669.48
    },
    {
        "displayLevel_": 0,
        "peakId_": "145",
        "pos_": 500.0371,
        "monoMz_": 500.0371,
        "intensity_": 669.16
    },
    {
        "displayLevel_": 0,
        "peakId_": "907",
        "pos_": 1725.8854,
        "monoMz_": 1725.8854,
        "intensity_": 668.05
    },
    {
        "displayLevel_": 0,
        "peakId_": "708",
        "pos_": 662.6848,
        "monoMz_": 662.6848,
        "intensity_": 667.09
    },
    {
        "displayLevel_": 0,
        "peakId_": "654",
        "pos_": 643.5372,
        "monoMz_": 643.5372,
        "intensity_": 663.81
    },
    {
        "displayLevel_": 0,
        "peakId_": "536",
        "pos_": 614.5922,
        "monoMz_": 614.5922,
        "intensity_": 661.88
    },
    {
        "displayLevel_": 0,
        "peakId_": "69",
        "pos_": 404.9106,
        "monoMz_": 404.9106,
        "intensity_": 661.24
    },
    {
        "displayLevel_": 0,
        "peakId_": "632",
        "pos_": 636.0557,
        "monoMz_": 636.0557,
        "intensity_": 660.9
    },
    {
        "displayLevel_": 5,
        "peakId_": "910",
        "pos_": 1842.212,
        "monoMz_": 1842.212,
        "intensity_": 660.23
    },
    {
        "displayLevel_": 0,
        "peakId_": "102",
        "pos_": 480.6033,
        "monoMz_": 480.6033,
        "intensity_": 658.62
    },
    {
        "displayLevel_": 0,
        "peakId_": "162",
        "pos_": 507.7725,
        "monoMz_": 507.7725,
        "intensity_": 657.45
    },
    {
        "displayLevel_": 0,
        "peakId_": "567",
        "pos_": 623.9819,
        "monoMz_": 623.9819,
        "intensity_": 657.45
    },
    {
        "displayLevel_": 0,
        "peakId_": "384",
        "pos_": 588.6943,
        "monoMz_": 588.6943,
        "intensity_": 657.25
    },
    {
        "displayLevel_": 0,
        "peakId_": "364",
        "pos_": 581.5333,
        "monoMz_": 581.5333,
        "intensity_": 657.09
    },
    {
        "displayLevel_": 0,
        "peakId_": "140",
        "pos_": 499.6154,
        "monoMz_": 499.6154,
        "intensity_": 655.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "27",
        "pos_": 199.3355,
        "monoMz_": 199.3355,
        "intensity_": 654.01
    },
    {
        "displayLevel_": 0,
        "peakId_": "120",
        "pos_": 488.2706,
        "monoMz_": 488.2706,
        "intensity_": 653.74
    },
    {
        "displayLevel_": 0,
        "peakId_": "112",
        "pos_": 486.3555,
        "monoMz_": 486.3555,
        "intensity_": 646.37
    },
    {
        "displayLevel_": 0,
        "peakId_": "605",
        "pos_": 630.3274,
        "monoMz_": 630.3274,
        "intensity_": 645.83
    },
    {
        "displayLevel_": 0,
        "peakId_": "21",
        "pos_": 199.2327,
        "monoMz_": 199.2327,
        "intensity_": 645.74
    },
    {
        "displayLevel_": 0,
        "peakId_": "425",
        "pos_": 600.3434,
        "monoMz_": 600.3434,
        "intensity_": 645.67
    },
    {
        "displayLevel_": 2,
        "peakId_": "906",
        "pos_": 1725.7805,
        "monoMz_": 1725.7805,
        "intensity_": 644.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "294",
        "pos_": 572.6415,
        "monoMz_": 572.6415,
        "intensity_": 643.94
    },
    {
        "displayLevel_": 0,
        "peakId_": "366",
        "pos_": 581.6207,
        "monoMz_": 581.6207,
        "intensity_": 643.17
    },
    {
        "displayLevel_": 4,
        "peakId_": "75",
        "pos_": 433.3825,
        "monoMz_": 433.3825,
        "intensity_": 641.92
    },
    {
        "displayLevel_": 0,
        "peakId_": "456",
        "pos_": 603.8902,
        "monoMz_": 603.8902,
        "intensity_": 641.86
    },
    {
        "displayLevel_": 4,
        "peakId_": "895",
        "pos_": 1326.2926,
        "monoMz_": 1326.2926,
        "intensity_": 641.09
    },
    {
        "displayLevel_": 0,
        "peakId_": "916",
        "pos_": 1877.561,
        "monoMz_": 1877.561,
        "intensity_": 638.71
    },
    {
        "displayLevel_": 0,
        "peakId_": "345",
        "pos_": 578.2129,
        "monoMz_": 578.2129,
        "intensity_": 637.53
    },
    {
        "displayLevel_": 0,
        "peakId_": "293",
        "pos_": 572.3082,
        "monoMz_": 572.3082,
        "intensity_": 634.07
    },
    {
        "displayLevel_": 0,
        "peakId_": "175",
        "pos_": 541.15,
        "monoMz_": 541.15,
        "intensity_": 632.49
    },
    {
        "displayLevel_": 0,
        "peakId_": "749",
        "pos_": 675.7811,
        "monoMz_": 675.7811,
        "intensity_": 632.02
    },
    {
        "displayLevel_": 0,
        "peakId_": "541",
        "pos_": 614.9539,
        "monoMz_": 614.9539,
        "intensity_": 631.77
    },
    {
        "displayLevel_": 0,
        "peakId_": "278",
        "pos_": 569.5243,
        "monoMz_": 569.5243,
        "intensity_": 629.14
    },
    {
        "displayLevel_": 4,
        "peakId_": "67",
        "pos_": 381.3012,
        "monoMz_": 381.3012,
        "intensity_": 628.6
    },
    {
        "displayLevel_": 0,
        "peakId_": "750",
        "pos_": 676.1672,
        "monoMz_": 676.1672,
        "intensity_": 628.04
    },
    {
        "displayLevel_": 5,
        "peakId_": "890",
        "pos_": 1112.9398,
        "monoMz_": 1112.9398,
        "intensity_": 627.63
    },
    {
        "displayLevel_": 0,
        "peakId_": "918",
        "pos_": 1904.5476,
        "monoMz_": 1904.5476,
        "intensity_": 626.04
    },
    {
        "displayLevel_": 0,
        "peakId_": "614",
        "pos_": 632.0109,
        "monoMz_": 632.0109,
        "intensity_": 624.4
    },
    {
        "displayLevel_": 3,
        "peakId_": "46",
        "pos_": 301.6928,
        "monoMz_": 301.6928,
        "intensity_": 624.3
    },
    {
        "displayLevel_": 1,
        "peakId_": "904",
        "pos_": 1723.5221,
        "monoMz_": 1723.5221,
        "intensity_": 623.18
    },
    {
        "displayLevel_": 0,
        "peakId_": "486",
        "pos_": 606.465,
        "monoMz_": 606.465,
        "intensity_": 622.25
    },
    {
        "displayLevel_": 0,
        "peakId_": "18",
        "pos_": 199.1392,
        "monoMz_": 199.1392,
        "intensity_": 620.91
    },
    {
        "displayLevel_": 4,
        "peakId_": "51",
        "pos_": 328.2308,
        "monoMz_": 328.2308,
        "intensity_": 619.26
    },
    {
        "displayLevel_": 0,
        "peakId_": "82",
        "pos_": 457.5914,
        "monoMz_": 457.5914,
        "intensity_": 616.4
    },
    {
        "displayLevel_": 0,
        "peakId_": "137",
        "pos_": 498.0892,
        "monoMz_": 498.0892,
        "intensity_": 613.25
    },
    {
        "displayLevel_": 0,
        "peakId_": "405",
        "pos_": 595.7214,
        "monoMz_": 595.7214,
        "intensity_": 612.91
    },
    {
        "displayLevel_": 0,
        "peakId_": "204",
        "pos_": 553.5104,
        "monoMz_": 553.5104,
        "intensity_": 612.8
    },
    {
        "displayLevel_": 0,
        "peakId_": "520",
        "pos_": 609.8367,
        "monoMz_": 609.8367,
        "intensity_": 612.69
    },
    {
        "displayLevel_": 0,
        "peakId_": "517",
        "pos_": 609.4613,
        "monoMz_": 609.4613,
        "intensity_": 612.44
    },
    {
        "displayLevel_": 0,
        "peakId_": "589",
        "pos_": 626.1285,
        "monoMz_": 626.1285,
        "intensity_": 608.71
    },
    {
        "displayLevel_": 0,
        "peakId_": "93",
        "pos_": 469.7137,
        "monoMz_": 469.7137,
        "intensity_": 606.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "252",
        "pos_": 562.9604,
        "monoMz_": 562.9604,
        "intensity_": 606.06
    },
    {
        "displayLevel_": 0,
        "peakId_": "551",
        "pos_": 621.105,
        "monoMz_": 621.105,
        "intensity_": 605.57
    },
    {
        "displayLevel_": 0,
        "peakId_": "201",
        "pos_": 553.1035,
        "monoMz_": 553.1035,
        "intensity_": 604.28
    },
    {
        "displayLevel_": 0,
        "peakId_": "563",
        "pos_": 622.8469,
        "monoMz_": 622.8469,
        "intensity_": 602.01
    },
    {
        "displayLevel_": 0,
        "peakId_": "911",
        "pos_": 1876.2551,
        "monoMz_": 1876.2551,
        "intensity_": 601.89
    },
    {
        "displayLevel_": 5,
        "peakId_": "29",
        "pos_": 213.7705,
        "monoMz_": 213.7705,
        "intensity_": 601.82
    },
    {
        "displayLevel_": 6,
        "peakId_": "898",
        "pos_": 1525.2131,
        "monoMz_": 1525.2131,
        "intensity_": 600.6
    },
    {
        "displayLevel_": 3,
        "peakId_": "76",
        "pos_": 434.172,
        "monoMz_": 434.172,
        "intensity_": 599.47
    },
    {
        "displayLevel_": 0,
        "peakId_": "265",
        "pos_": 566.006,
        "monoMz_": 566.006,
        "intensity_": 597.96
    },
    {
        "displayLevel_": 0,
        "peakId_": "393",
        "pos_": 590.0958,
        "monoMz_": 590.0958,
        "intensity_": 597.29
    },
    {
        "displayLevel_": 1,
        "peakId_": "757",
        "pos_": 679.5858,
        "monoMz_": 679.5858,
        "intensity_": 597.18
    },
    {
        "displayLevel_": 0,
        "peakId_": "378",
        "pos_": 587.3216,
        "monoMz_": 587.3216,
        "intensity_": 597.1
    },
    {
        "displayLevel_": 0,
        "peakId_": "648",
        "pos_": 642.0369,
        "monoMz_": 642.0369,
        "intensity_": 596.59
    },
    {
        "displayLevel_": 3,
        "peakId_": "924",
        "pos_": 1984.3383,
        "monoMz_": 1984.3383,
        "intensity_": 595.39
    },
    {
        "displayLevel_": 0,
        "peakId_": "715",
        "pos_": 665.6821,
        "monoMz_": 665.6821,
        "intensity_": 594.2
    },
    {
        "displayLevel_": 5,
        "peakId_": "7",
        "pos_": 157.582,
        "monoMz_": 157.582,
        "intensity_": 593.59
    },
    {
        "displayLevel_": 0,
        "peakId_": "534",
        "pos_": 613.8014,
        "monoMz_": 613.8014,
        "intensity_": 593.57
    },
    {
        "displayLevel_": 2,
        "peakId_": "78",
        "pos_": 452.6884,
        "monoMz_": 452.6884,
        "intensity_": 590.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "591",
        "pos_": 626.3349,
        "monoMz_": 626.3349,
        "intensity_": 589.32
    },
    {
        "displayLevel_": 0,
        "peakId_": "585",
        "pos_": 625.7938,
        "monoMz_": 625.7938,
        "intensity_": 588.87
    },
    {
        "displayLevel_": 0,
        "peakId_": "246",
        "pos_": 560.9136,
        "monoMz_": 560.9136,
        "intensity_": 588.81
    },
    {
        "displayLevel_": 0,
        "peakId_": "438",
        "pos_": 602.1768,
        "monoMz_": 602.1768,
        "intensity_": 588.14
    },
    {
        "displayLevel_": 0,
        "peakId_": "545",
        "pos_": 616.0377,
        "monoMz_": 616.0377,
        "intensity_": 587.86
    },
    {
        "displayLevel_": 0,
        "peakId_": "74",
        "pos_": 433.2406,
        "monoMz_": 433.2406,
        "intensity_": 586.92
    },
    {
        "displayLevel_": 0,
        "peakId_": "864",
        "pos_": 769.6194,
        "monoMz_": 769.6194,
        "intensity_": 586.26
    },
    {
        "displayLevel_": 4,
        "peakId_": "32",
        "pos_": 224.6413,
        "monoMz_": 224.6413,
        "intensity_": 581.64
    },
    {
        "displayLevel_": 0,
        "peakId_": "472",
        "pos_": 605.4672,
        "monoMz_": 605.4672,
        "intensity_": 577.86
    },
    {
        "displayLevel_": 0,
        "peakId_": "756",
        "pos_": 679.3436,
        "monoMz_": 679.3436,
        "intensity_": 577.29
    },
    {
        "displayLevel_": 0,
        "peakId_": "723",
        "pos_": 668.576,
        "monoMz_": 668.576,
        "intensity_": 576.89
    },
    {
        "displayLevel_": 0,
        "peakId_": "290",
        "pos_": 571.8633,
        "monoMz_": 571.8633,
        "intensity_": 576.73
    },
    {
        "displayLevel_": 2,
        "peakId_": "125",
        "pos_": 490.0344,
        "monoMz_": 490.0344,
        "intensity_": 574.2
    },
    {
        "displayLevel_": 4,
        "peakId_": "872",
        "pos_": 821.3143,
        "monoMz_": 821.3143,
        "intensity_": 572.54
    },
    {
        "displayLevel_": 0,
        "peakId_": "421",
        "pos_": 598.9524,
        "monoMz_": 598.9524,
        "intensity_": 572.5
    },
    {
        "displayLevel_": 0,
        "peakId_": "47",
        "pos_": 301.7193,
        "monoMz_": 301.7193,
        "intensity_": 571.44
    },
    {
        "displayLevel_": 3,
        "peakId_": "901",
        "pos_": 1713.4988,
        "monoMz_": 1713.4988,
        "intensity_": 571.29
    },
    {
        "displayLevel_": 6,
        "peakId_": "886",
        "pos_": 990.0797,
        "monoMz_": 990.0797,
        "intensity_": 570.15
    },
    {
        "displayLevel_": 0,
        "peakId_": "251",
        "pos_": 562.1769,
        "monoMz_": 562.1769,
        "intensity_": 569.7
    },
    {
        "displayLevel_": 0,
        "peakId_": "547",
        "pos_": 620.2146,
        "monoMz_": 620.2146,
        "intensity_": 569.38
    },
    {
        "displayLevel_": 1,
        "peakId_": "8",
        "pos_": 159.4742,
        "monoMz_": 159.4742,
        "intensity_": 569.31
    },
    {
        "displayLevel_": 1,
        "peakId_": "909",
        "pos_": 1727.8619,
        "monoMz_": 1727.8619,
        "intensity_": 568.35
    },
    {
        "displayLevel_": 0,
        "peakId_": "256",
        "pos_": 563.3289,
        "monoMz_": 563.3289,
        "intensity_": 568.13
    },
    {
        "displayLevel_": 0,
        "peakId_": "426",
        "pos_": 600.4678,
        "monoMz_": 600.4678,
        "intensity_": 566.3
    },
    {
        "displayLevel_": 0,
        "peakId_": "255",
        "pos_": 563.282,
        "monoMz_": 563.282,
        "intensity_": 564.97
    },
    {
        "displayLevel_": 0,
        "peakId_": "61",
        "pos_": 365.9381,
        "monoMz_": 365.9381,
        "intensity_": 564.57
    },
    {
        "displayLevel_": 0,
        "peakId_": "259",
        "pos_": 564.442,
        "monoMz_": 564.442,
        "intensity_": 562.79
    },
    {
        "displayLevel_": 4,
        "peakId_": "13",
        "pos_": 175.7913,
        "monoMz_": 175.7913,
        "intensity_": 562.2
    },
    {
        "displayLevel_": 0,
        "peakId_": "85",
        "pos_": 462.3262,
        "monoMz_": 462.3262,
        "intensity_": 562.03
    },
    {
        "displayLevel_": 5,
        "peakId_": "897",
        "pos_": 1383.515,
        "monoMz_": 1383.515,
        "intensity_": 561.97
    },
    {
        "displayLevel_": 0,
        "peakId_": "353",
        "pos_": 579.6936,
        "monoMz_": 579.6936,
        "intensity_": 561.74
    },
    {
        "displayLevel_": 0,
        "peakId_": "135",
        "pos_": 497.7816,
        "monoMz_": 497.7816,
        "intensity_": 559.54
    },
    {
        "displayLevel_": 0,
        "peakId_": "311",
        "pos_": 574.7876,
        "monoMz_": 574.7876,
        "intensity_": 559.36
    },
    {
        "displayLevel_": 0,
        "peakId_": "30",
        "pos_": 214.4859,
        "monoMz_": 214.4859,
        "intensity_": 558.5
    },
    {
        "displayLevel_": 3,
        "peakId_": "10",
        "pos_": 163.9326,
        "monoMz_": 163.9326,
        "intensity_": 557.71
    },
    {
        "displayLevel_": 0,
        "peakId_": "577",
        "pos_": 625.1451,
        "monoMz_": 625.1451,
        "intensity_": 555.47
    },
    {
        "displayLevel_": 4,
        "peakId_": "53",
        "pos_": 335.4095,
        "monoMz_": 335.4095,
        "intensity_": 555
    },
    {
        "displayLevel_": 0,
        "peakId_": "439",
        "pos_": 602.2054,
        "monoMz_": 602.2054,
        "intensity_": 553.17
    },
    {
        "displayLevel_": 0,
        "peakId_": "217",
        "pos_": 556.1461,
        "monoMz_": 556.1461,
        "intensity_": 547.18
    },
    {
        "displayLevel_": 0,
        "peakId_": "264",
        "pos_": 565.9903,
        "monoMz_": 565.9903,
        "intensity_": 546.96
    },
    {
        "displayLevel_": 7,
        "peakId_": "5",
        "pos_": 150.8699,
        "monoMz_": 150.8699,
        "intensity_": 543.76
    },
    {
        "displayLevel_": 4,
        "peakId_": "887",
        "pos_": 1040.325,
        "monoMz_": 1040.325,
        "intensity_": 542.24
    },
    {
        "displayLevel_": 0,
        "peakId_": "200",
        "pos_": 552.9795,
        "monoMz_": 552.9795,
        "intensity_": 542.2
    },
    {
        "displayLevel_": 1,
        "peakId_": "15",
        "pos_": 177.7424,
        "monoMz_": 177.7424,
        "intensity_": 542.07
    },
    {
        "displayLevel_": 0,
        "peakId_": "588",
        "pos_": 626.0443,
        "monoMz_": 626.0443,
        "intensity_": 540.54
    },
    {
        "displayLevel_": 3,
        "peakId_": "66",
        "pos_": 374.5431,
        "monoMz_": 374.5431,
        "intensity_": 539.73
    },
    {
        "displayLevel_": 0,
        "peakId_": "41",
        "pos_": 279.769,
        "monoMz_": 279.769,
        "intensity_": 539.25
    },
    {
        "displayLevel_": 3,
        "peakId_": "166",
        "pos_": 514.2649,
        "monoMz_": 514.2649,
        "intensity_": 538.89
    },
    {
        "displayLevel_": 0,
        "peakId_": "14",
        "pos_": 176.4013,
        "monoMz_": 176.4013,
        "intensity_": 538.67
    },
    {
        "displayLevel_": 0,
        "peakId_": "52",
        "pos_": 328.4895,
        "monoMz_": 328.4895,
        "intensity_": 538.24
    },
    {
        "displayLevel_": 0,
        "peakId_": "564",
        "pos_": 622.9882,
        "monoMz_": 622.9882,
        "intensity_": 537.86
    },
    {
        "displayLevel_": 0,
        "peakId_": "16",
        "pos_": 198.7277,
        "monoMz_": 198.7277,
        "intensity_": 537.13
    },
    {
        "displayLevel_": 0,
        "peakId_": "349",
        "pos_": 578.6416,
        "monoMz_": 578.6416,
        "intensity_": 535.72
    },
    {
        "displayLevel_": 0,
        "peakId_": "590",
        "pos_": 626.2078,
        "monoMz_": 626.2078,
        "intensity_": 535.18
    },
    {
        "displayLevel_": 0,
        "peakId_": "905",
        "pos_": 1725.7031,
        "monoMz_": 1725.7031,
        "intensity_": 534.85
    },
    {
        "displayLevel_": 0,
        "peakId_": "130",
        "pos_": 493.5882,
        "monoMz_": 493.5882,
        "intensity_": 534.06
    },
    {
        "displayLevel_": 5,
        "peakId_": "893",
        "pos_": 1233.8242,
        "monoMz_": 1233.8242,
        "intensity_": 534.04
    },
    {
        "displayLevel_": 0,
        "peakId_": "323",
        "pos_": 576.4167,
        "monoMz_": 576.4167,
        "intensity_": 533.05
    },
    {
        "displayLevel_": 0,
        "peakId_": "896",
        "pos_": 1381.8534,
        "monoMz_": 1381.8534,
        "intensity_": 532.89
    },
    {
        "displayLevel_": 0,
        "peakId_": "156",
        "pos_": 506.3084,
        "monoMz_": 506.3084,
        "intensity_": 532.52
    },
    {
        "displayLevel_": 1,
        "peakId_": "6",
        "pos_": 153.0671,
        "monoMz_": 153.0671,
        "intensity_": 532.32
    },
    {
        "displayLevel_": 0,
        "peakId_": "385",
        "pos_": 589.0657,
        "monoMz_": 589.0657,
        "intensity_": 532.23
    },
    {
        "displayLevel_": 1,
        "peakId_": "45",
        "pos_": 294.4229,
        "monoMz_": 294.4229,
        "intensity_": 531.9
    },
    {
        "displayLevel_": 0,
        "peakId_": "177",
        "pos_": 542.0922,
        "monoMz_": 542.0922,
        "intensity_": 531.2
    },
    {
        "displayLevel_": 3,
        "peakId_": "49",
        "pos_": 320.9097,
        "monoMz_": 320.9097,
        "intensity_": 528.88
    },
    {
        "displayLevel_": 0,
        "peakId_": "40",
        "pos_": 278.9547,
        "monoMz_": 278.9547,
        "intensity_": 528.65
    },
    {
        "displayLevel_": 6,
        "peakId_": "2",
        "pos_": 138.251,
        "monoMz_": 138.251,
        "intensity_": 528.44
    },
    {
        "displayLevel_": 0,
        "peakId_": "355",
        "pos_": 579.7732,
        "monoMz_": 579.7732,
        "intensity_": 528.23
    },
    {
        "displayLevel_": 3,
        "peakId_": "0",
        "pos_": 135.8057,
        "monoMz_": 135.8057,
        "intensity_": 527.11
    },
    {
        "displayLevel_": 0,
        "peakId_": "11",
        "pos_": 165.6486,
        "monoMz_": 165.6486,
        "intensity_": 522.91
    },
    {
        "displayLevel_": 0,
        "peakId_": "106",
        "pos_": 482.5889,
        "monoMz_": 482.5889,
        "intensity_": 521.87
    },
    {
        "displayLevel_": 2,
        "peakId_": "73",
        "pos_": 429.1618,
        "monoMz_": 429.1618,
        "intensity_": 520.16
    },
    {
        "displayLevel_": 2,
        "peakId_": "12",
        "pos_": 166.3798,
        "monoMz_": 166.3798,
        "intensity_": 519.94
    },
    {
        "displayLevel_": 0,
        "peakId_": "925",
        "pos_": 1993.4447,
        "monoMz_": 1993.4447,
        "intensity_": 519.87
    },
    {
        "displayLevel_": 0,
        "peakId_": "383",
        "pos_": 588.4418,
        "monoMz_": 588.4418,
        "intensity_": 518.23
    },
    {
        "displayLevel_": 1,
        "peakId_": "56",
        "pos_": 349.8962,
        "monoMz_": 349.8962,
        "intensity_": 516.99
    },
    {
        "displayLevel_": 2,
        "peakId_": "36",
        "pos_": 240.5779,
        "monoMz_": 240.5779,
        "intensity_": 514.2
    },
    {
        "displayLevel_": 4,
        "peakId_": "38",
        "pos_": 263.5172,
        "monoMz_": 263.5172,
        "intensity_": 513.07
    },
    {
        "displayLevel_": 2,
        "peakId_": "50",
        "pos_": 325.5329,
        "monoMz_": 325.5329,
        "intensity_": 509.41
    },
    {
        "displayLevel_": 2,
        "peakId_": "31",
        "pos_": 216.6362,
        "monoMz_": 216.6362,
        "intensity_": 506.63
    },
    {
        "displayLevel_": 2,
        "peakId_": "9",
        "pos_": 160.0527,
        "monoMz_": 160.0527,
        "intensity_": 506.48
    },
    {
        "displayLevel_": 0,
        "peakId_": "397",
        "pos_": 593.0668,
        "monoMz_": 593.0668,
        "intensity_": 503.88
    },
    {
        "displayLevel_": 0,
        "peakId_": "331",
        "pos_": 577.0149,
        "monoMz_": 577.0149,
        "intensity_": 502.12
    },
    {
        "displayLevel_": 1,
        "peakId_": "54",
        "pos_": 338.1656,
        "monoMz_": 338.1656,
        "intensity_": 499.18
    },
    {
        "displayLevel_": 3,
        "peakId_": "33",
        "pos_": 235.4364,
        "monoMz_": 235.4364,
        "intensity_": 495.05
    },
    {
        "displayLevel_": 0,
        "peakId_": "215",
        "pos_": 555.8754,
        "monoMz_": 555.8754,
        "intensity_": 494.04
    },
    {
        "displayLevel_": 1,
        "peakId_": "34",
        "pos_": 238.286,
        "monoMz_": 238.286,
        "intensity_": 490.65
    },
    {
        "displayLevel_": 5,
        "peakId_": "885",
        "pos_": 924.889,
        "monoMz_": 924.889,
        "intensity_": 489.57
    },
    {
        "displayLevel_": 0,
        "peakId_": "417",
        "pos_": 598.5823,
        "monoMz_": 598.5823,
        "intensity_": 487.8
    },
    {
        "displayLevel_": 2,
        "peakId_": "3",
        "pos_": 143.5742,
        "monoMz_": 143.5742,
        "intensity_": 487.29
    },
    {
        "displayLevel_": 3,
        "peakId_": "43",
        "pos_": 283.8446,
        "monoMz_": 283.8446,
        "intensity_": 478.34
    },
    {
        "displayLevel_": 0,
        "peakId_": "1",
        "pos_": 135.9892,
        "monoMz_": 135.9892,
        "intensity_": 475.06
    },
    {
        "displayLevel_": 3,
        "peakId_": "4",
        "pos_": 145.6313,
        "monoMz_": 145.6313,
        "intensity_": 475.02
    },
    {
        "displayLevel_": 3,
        "peakId_": "39",
        "pos_": 271.0363,
        "monoMz_": 271.0363,
        "intensity_": 472.84
    },
    {
        "displayLevel_": 3,
        "peakId_": "77",
        "pos_": 446.9713,
        "monoMz_": 446.9713,
        "intensity_": 465.42
    },
    {
        "displayLevel_": 0,
        "peakId_": "35",
        "pos_": 240.1386,
        "monoMz_": 240.1386,
        "intensity_": 463.56
    },
    {
        "displayLevel_": 0,
        "peakId_": "556",
        "pos_": 622.0079,
        "monoMz_": 622.0079,
        "intensity_": 458.82
    },
    {
        "displayLevel_": 0,
        "peakId_": "332",
        "pos_": 577.022,
        "monoMz_": 577.022,
        "intensity_": 454.77
    },
    {
        "displayLevel_": 0,
        "peakId_": "253",
        "pos_": 563.1798,
        "monoMz_": 563.1798,
        "intensity_": 454.43
    }
];
