var CSQ = function (info_csq, csq_names, numeric_fields) {
    var csq = this;
    csq._values = info_csq.split('|');
    csq._numeric_fields = numeric_fields;
    // create an object of key => index from vcf_csq
    csq._names = {};
    csq_names.forEach(function (name, index) { csq._names[name] = index; });

    return new Proxy(csq, {
        get: function (target, name) {
            if (!(name in target._names)) {
                throw "unknown CSQ field: " + name;
            }

            var result = target._values[target._names[name]];
            if (target._numeric_fields[name]) {
                if(result.length == 0 || isNaN(result)) {
                    result = NaN
                } else {
                    try {
                        result = parseFloat(result);
                    } catch(e) { result = NaN }
                }
            }
            return result
        }
    });
};

function CSQs(csq_string, vcf_csq, numeric_fields) {
    // if numeric fields is an array, turn it into an object
    if (Array.isArray(numeric_fields)) {
        var numeric_fields_obj = {};
        numeric_fields.forEach(function (field) { numeric_fields_obj[field] = true; });
        numeric_fields = numeric_fields_obj;
    }

    csqs = csq_string.split(',');
    return csqs.map(function (csq) { return new CSQ(csq, vcf_csq, numeric_fields); });
}
