from pytest import approx

from astropc.moon import apparent, lunar_node, true_position


class TestTruePosition:
    delta = 1e-4  # result precision
    cases = (
        {
            "djd": -1.000050e04,
            "coords": [253.85478, -0.35884, 0.98681],
            "delta": 0.002475,
            "motion": 14.073505,
        },
        {
            "djd": -7.000500e03,
            "coords": [183.03298, -5.10613, 0.96482],
            "delta": 0.0025318451878263725,
            "motion": 13.614904285991807,
        },
        {
            "djd": -4.000500e03,
            "coords": [114.49714, 0.29899, 0.91786],
            "delta": 0.002661387458557927,
            "motion": 12.284203442108854,
        },
        {
            "djd": -1.000500e03,
            "coords": [46.33258, 5.03904, 0.89971],
            "delta": 0.0027150753763781643,
            "motion": 11.86016463804351,
        },
        {
            "djd": 1.999500e03,
            "coords": [340.74811, -0.76686, 0.91660],
            "delta": 0.002665042330735118,
            "motion": 12.137096046101872,
        },
        {
            "djd": 4.999500e03,
            "coords": [273.11888, -5.22297, 0.93431],
            "delta": 0.0026145243283283597,
            "motion": 12.706509343283184,
        },
        {
            "djd": 7.999500e03,
            "coords": [198.76809, 0.13467, 0.97476],
            "delta": 0.0025060386483159646,
            "motion": 13.79049510733593,
        },
        {
            "djd": 1.099950e04,
            "coords": [123.17331, 5.01217, 1.02067],
            "delta": 0.002393311553917354,
            "motion": 15.25014893599962,
        },
        {
            "djd": 1.399950e04,
            "coords": [50.40519, 0.59539, 1.00077],
            "delta": 0.002440903741394359,
            "motion": 14.567332957243783,
        },
        {
            "djd": 1.699950e04,
            "coords": [336.88148, -5.04905, 0.94329],
            "delta": 0.0025896311772431097,
            "motion": 13.015006558327384,
        },
        {
            "djd": 1.999950e04,
            "coords": [266.43192, -1.18331, 0.91398],
            "delta": 0.0026726946153555506,
            "motion": 12.05705112860313,
        },
        {
            "djd": 2.299950e04,
            "coords": [200.91657, 5.13843, 0.90354],
            "delta": 0.00270357511434672,
            "motion": 11.883519914105939,
        },
        {
            "djd": 2.599950e04,
            "coords": [134.05765, 0.87204, 0.90670],
            "delta": 0.0026941433274419316,
            "motion": 11.945823078908266,
        },
        {
            "djd": 2.899950e04,
            "coords": [64.16216, -4.94147, 0.94934],
            "delta": 0.0025731373392409293,
            "motion": 13.2314091077357,
        },
        {
            "djd": 3.199950e04,
            "coords": [354.53313, -0.77311, 0.99650],
            "delta": 0.0024513561792419898,
            "motion": 14.398538661582212,
        },
        {
            "djd": 3.499950e04,
            "coords": [280.10165, 5.06817, 0.99501],
            "delta": 0.002455022531559789,
            "motion": 14.431229034360273,
        },
        {
            "djd": 3.799950e04,
            "coords": [201.62149, 2.25573, 0.97435],
            "delta": 0.0025070947174279036,
            "motion": 13.731560363493482,
        },
        {
            "djd": 4.099950e04,
            "coords": [128.41649, -4.51661, 0.95591],
            "delta": 0.0025554365866768364,
            "motion": 13.279315343224834,
        },
        {
            "djd": 4.399950e04,
            "coords": [61.54198, -2.45092, 0.92162],
            "delta": 0.0026505164857345337,
            "motion": 12.374443595332336,
        },
        {
            "djd": 4.699950e04,
            "coords": [353.93133, 4.49791, 0.89930],
            "delta": 0.00271630014162966,
            "motion": 11.857387239871063,
        },
    )

    def test_lambda(self):
        for case in self.cases:
            pos = true_position(case["djd"])
            assert approx(pos.lmbda, abs=self.delta) == case["coords"][0]

    def test_beta(self):
        for case in self.cases:
            pos = true_position(case["djd"])
            assert approx(pos.beta, abs=self.delta) == case["coords"][1]

    def test_parallax(self):
        for case in self.cases:
            pos = true_position(case["djd"])
            assert approx(pos.parallax, abs=self.delta) == case["coords"][2]

    def test_delta(self):
        for case in self.cases:
            pos = true_position(case["djd"])
            assert approx(pos.delta, abs=self.delta) == case["delta"]


class TestApparent:
    delta = 1e-4  # result precision
    cases = (
        {"djd": 23772.99027777778, "lng": 310.19998902960941},
        {"djd": 30735.5, "lng": 260.7128333333333},  # 1984-2-25.0
        {"djd": 16773.8121, "lng": 246.94925},  # 1945-12-4.3121,
    )

    def test_longitude(self):
        for c in self.cases:
            pos = apparent(c["djd"])
            assert approx(pos.lmbda, abs=self.delta) == c["lng"]


class TestLunarNode:
    delta = 1e-4  # result precision
    djd = 23772.99027777778

    def test_mean_node(self):
        assert (
            approx(lunar_node(self.djd, true_node=False), abs=self.delta)
            == 80.31173473979322
        )

    def test_true_node(self):
        assert (
            approx(lunar_node(self.djd, true_node=True), abs=self.delta)
            == 81.86652882901491
        )
