function trio_denovo(kid, dad, mom) {
    // alts are 0:hom_ref, 1:het, 2:hom_alt, -1:unknown
    if(!(kid.alts == 1 && mom.alts == 0 && dad.alts == 0)){ return false}
    // sufficient depth in all
    if(kid.DP < 7 || mom.DP < 7 || dad.DP < 7) { return false; }
    // no evidence for alternate in the parents
    if((mom.AD[1] + dad.AD[1]) > 0) { return false; }
    // check the kid's allele balance.
    if(kid.AB < 0.2 || kid.AB > 0.8) { return false; }
    return true
}

function hq(sample) {
	// this function checks that the genotype (alts) is consistent with the information
	// and that the depth, GQ and allele balance are good.
	if(sample.alts == -1) { return false; }
	if(sample.DP < 8) { return false; }
	if(sample.GQ < 10) { return false; }
	if(sample.alts == 0) {
		// if there is more than 1 piece of evidence for the alt allele, it's not HQ
		if(sample.DP > 20 && sample.AB > 0.02) { return false; }
		if(sample.DP <= 20 && sample.AD[1] > 1) { return false; }
		return true
	}
	if(sample.alts == 1) {
        if(sample.AB < 0.2 || sample.AB > 0.8) { return false; }
		return true
	}
	if(sample.alts == 2) {
		if(sample.DP > 20 && sample.AB < 0.98) { return false; }
		if(sample.DP <= 20 && sample.AD[0] > 1) { return false; }
		return true
	}
}

function hiqual(kid, dad, mom) {
	return hq(kid) && hq(dad) && hq(mom)
}


function trio_hets(kid, dad, mom) {
    if(kid.DP < 7 || mom.DP < 7 || dad.DP < 7) { return false; }
    if(kid.GQ < 10 || mom.GQ < 10 || dad.GQ < 10) { return false; }
	if(kid.alts == 0){ return false; }
	if(kid.alts == 1 && mom.alts == 1 && dad.alts == 1){ return false; }
	// explicitly allow this
	if(kid.alts == 2 && mom.alts == 1 && dad.alts == 1){ return true; }
	return true
}

function trio_autosomal_dominant(kid, dad, mom) {
    // affected samples must be het.
    if(!(kid.affected == (kid.alts == 1) && mom.affected == (mom.alts == 1) && dad.affected == (dad.alts == 1))) { return false; }
    return kid.affected && hiqual(kid, dad, mom)
}

function trio_autosomal_recessive(kid, dad, mom) {
	return kid.affected && kid.alts == 2 && mom.alts == 1 && dad.alts == 1 && hiqual(kid, dad, mom)
}

function trio_x_linked_recessive(kid, dad, mom) {
  if(!hiqual(kid, dad, mom)) { return false; }
  if(kid.alts == 0 || kid.alts == -1) { return false; }
  if((dad.alts != 0) != dad.affected) { return false; }
  if(mom.alts != 1) { return false; }
  return kid.affected
}

function trio_x_linked_denovo(kid, dad, mom) {
  if(!hiqual(kid, dad, mom)) { return false; }
  if(kid.sex == "unknown"){ return false; }
  if(!(mom.alts == 0 && dad.alts == 0)){ return false}
  if(kid.sex == "male") {
    return kid.alts == 1 || kid.alts == 2;
  }
  // female
  return kid.alts == 1;
}
