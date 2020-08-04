class IonMassShift {
    ionType;

    constructor(ionType) {
        this.ionType = ionType;
    }

    getIonTypeMass(){
		let ionTypeMassList={
			"A":-27.9949,
			"A-H2O":-46.0149,
			"A-NH3":-45.02542,
			"B":0,
			"B-H2O":-18.02,
			"B-NH3":-17.03052,
			"C":17.0266,
			"C-H2O":-0.9934,
			"C-NH3":-0.00392,
			"X":43.99,
			"X-H2O":25.97,
			"X-NH3":26.95948,
			"Y":18.0106,
			"Y-H2O":0,
			"Y-NH3":0.98008,
			"Z":0.984,
			"Z-H2O":-17.026,
			"Z-NH3":-16.04652,
			"Z0":1.9919,
			"Z0-H2O":-16.018664,
			"Z0-NH3":-15.03862
        };
        if (ionTypeMassList.hasOwnProperty(this.ionType.toUpperCase())) {
            return ionTypeMassList[this.ionType.toUpperCase()];
        } else {
            return undefined;
        }
	}
}