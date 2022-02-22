//class for formatting string or number
class FormatUtil {
  constructor(){};
  static formatFloat(number: string | number, type:string): string {
    let formatInfo: {[key: string]: number} = {
      "protMass": 3,
      "massDiff": 3,
      "massShift": 3,
      "precMz": 3,
      "precMass": 3,
      "ppmError": 2,
      "dataTable":4
    };

    if (typeof(number) == "string") {
        number = parseFloat(number);
        if (isNaN(number)) {
          console.error("ERROR: invalid number", number, " passed for formating");
          return "-1";
        }
    }
    return number.toFixed(formatInfo[type]);
  }
}