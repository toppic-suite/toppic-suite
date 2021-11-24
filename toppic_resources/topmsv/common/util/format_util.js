"use strict";
//class for formatting string or number
class FormatUtil {
    constructor() { }
    ;
    static formatFloat(number, digit) {
        if (typeof (number) == "string") {
            number = parseFloat(number);
            if (isNaN(number)) {
                console.error("ERROR: invalid number passed for formating");
                return "-1";
            }
        }
        return number.toFixed(digit);
    }
}
