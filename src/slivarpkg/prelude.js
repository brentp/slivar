function hasSample(obj, key, sample) {
    // test if a sample has been added to obj[key] which is a comma-delimited string
    if(!(key in obj)) { return false; }
    if(obj[key].indexOf(sample) == -1){ return false; }
    if(obj[key].indexOf(',' + sample + ',') != -1){ return true; }
    if(obj[key] == sample){ return true; }
    var idx = obj[key].indexOf(',' + sample);

    // test if the sample is the last one in the string
    if(idx != -1 && idx + sample.length + 1 == obj[key].length) {
      return true;
    }
    return (obj[key].indexOf(sample + ',') == 0);

}
